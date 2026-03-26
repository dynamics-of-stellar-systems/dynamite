import logging
import math
import numpy as np
import os
from scipy.integrate import quad
from scipy import special
from astropy import table
from dynamite import data
from dynamite import constants

class MGE(data.Data):
    """Multi Gaussian Expansions"""

    def __init__(self, config, **kwargs):
        super().__init__(**kwargs)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.config = config
        self.validate_q_values()

    def validate_q_values(self):
        """Validates the mge's q values

        Any q 'too close to 1' will be set to q=NINES for numerical stability.
        If any changes are made, a warning message will be logged. Note that
        the 'closeness' to 1 might be machine dependent - for us 'four nines',
        0.9999, worked...

        Returns
        -------
        None. Any changes are applied to ``self.data``.

        """
        NINES = 0.9999
        new_mge = False
        for r in self.data:
            if r['q'] > NINES:
                self.logger.warning(f'changing q={r["q"]} to q={NINES} for '
                                    'numerical stability.')
                r['q'] = NINES
                new_mge = True
        if new_mge:
            self.logger.warning(f'New mge:\n{self.data}')

    def read_file_old_format(self, filename):
        """read the MGE data from a text file

        old format = text file, 4 columns (I, sigma, q, PA_twist), with one
        header line

        Parameters
        ----------
        filename : string
            name of file

        Returns
        -------
        ndarray
            MGE data in structrued numpy array

        """
        with open(filename) as fp:
            n_cmp_mge = fp.readline()
        n_cmp_mge = int(n_cmp_mge)
        dat = np.genfromtxt(filename,
                            skip_header=1,
                            max_rows=n_cmp_mge,
                            names=['I', 'sigma', 'q', 'PA_twist'])
        if np.isnan(dat).any():
            txt = f'Input file {filename} has nans'
            self.logger.error(f'{txt} at: {np.argwhere(np.isnan(dat))}.')
            raise ValueError(txt)
        return dat

    def convert_file_from_old_format(self,
                                     filename_old_format,
                                     filename_new_format):
        """convert old mge file to ECSV file

        Parameters
        ----------
        filename_old_format : string
            old filename
        filename_new_format : string
            new filename

        Returns
        -------
        None
            saves file with name ``filename_new_format``

        """
        data = self.read_file_old_format(filename_old_format)
        data = table.Table(data)
        data.write(filename_new_format, format='ascii.ecsv')

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__dict__})'

    def get_projected_masses(self, use_cache=True):
        """Calculate the mass of the mge in observed 2D apertures.

        Calculate the mass of the mge in observed 2D apertures using observed
        quantities only. Hence, this calculation is independent of the model's
        hyperparameters. The projected masses are written to mass_aper.ecsv
        in the output directory. If called again, this method will read
        previously calculated projected masses from the visible component's
        mass_aper attribute or from file unless use_cache=False.

        Parameters
        ----------
        use_cache : bool, optional
            If False, the projected masses will be recalculated.
            If True, check for projected masses already existing in the visible
            component's mass_aper attribute or on disk and if yes, return those
            projected masses. The mass_aper attribute will override data on
            disk. If data is read from disk, it will be used to update the
            visible component's mass_aper attribute. The default is True.

        Returns
        -------
        numpy array, shape=(n_spatial_bins,)
            Aperture masses of the MGE
        """
        c = self.config
        if c.system.is_bar_disk_system():
            vis_comp = c.system.get_unique_bar_component()
        else:
            vis_comp = c.system.get_unique_triaxial_visible_component()
        if use_cache and (vis_comp.mass_aper is not None):
            self.logger.info('Projected masses grabbed from vis. component.')
            return vis_comp.mass_aper  ########################################
        p_mass_fname = c.settings.io_settings['output_directory'] + \
                          'mass_aper.ecsv'
        if use_cache and os.path.isfile(p_mass_fname):
            proj_mass = table.Table.read(p_mass_fname, format='ascii')
            self.logger.info(f'Projected masses read from {p_mass_fname}.')
            vis_comp.mass_aper = np.array(proj_mass['proj_mass'])
            return vis_comp.mass_aper  ########################################
        kinematics = vis_comp.kinematic_data
        distMPc = c.system.distMPc
        arcsec_to_km = constants.ARC_KM(distMPc)
        eps_abs = eps_rel = 1.49e-08
        for i_kin, kin in enumerate(kinematics):
            self.logger.info(f'Calculating projected masses for {kin.name}...')
            min_x = kin.dp_args['min_x'] * arcsec_to_km
            min_y = kin.dp_args['min_y'] * arcsec_to_km
            x_size = kin.dp_args['x_size'] * arcsec_to_km
            y_size = kin.dp_args['y_size'] * arcsec_to_km
            n_x = kin.dp_args['n_x']
            n_y = kin.dp_args['n_y']
            angle = kin.dp_args['angle'] * math.pi / 180

            psf_weight = kin.PSF['weight']
            psf_width = [s * arcsec_to_km for s in kin.PSF['sigma']]

            # Boundaries in x and y direction
            bx = np.linspace(min_x, min_x + x_size, num=n_x + 1, endpoint=True)
            by = np.linspace(min_y, min_y + y_size, num=n_y + 1, endpoint=True)

            # Dispersion in km
            sigobs_km = self.data['sigma'] * arcsec_to_km
            # Surface brightness in L_sun/km^2 (we don't multiply by ml here)
            surf_km = self.data['I'] / constants.PARSEC_KM ** 2
            # Observed flattening
            qobs = self.data['q']
            # Offset psi viewing angle in radians
            psi_obs = self.data['PA_twist'] * math.pi / 180
            # psi_view = psi * math.pi / 180
            # psi_obs = psi_obs + psi_view
            isotwist = psi_obs  # Fortran: = psi_obs + psi_view - psi_view

            grid = np.zeros((n_x, n_y))
            # Consider parallelizing in i, j below
            for i, j in [(i, j) for i in range(n_x) for j in range(n_y)]:
                for k in range(len(self.data)):
                    for psf_idx in range(len(psf_width)):
                        res = quad(self._integrand,
                                    bx[i],
                                    bx[i + 1],
                                    args=(by[j],
                                          by[j + 1],
                                          angle,
                                          qobs[k],
                                          sigobs_km[k],
                                          surf_km[k],
                                          isotwist[k],
                                          psf_width[psf_idx],
                                          psf_weight[psf_idx]),
                                    epsabs=eps_abs,
                                    epsrel=eps_rel,
                                    limit=50)
                        if res[1] > 2 * max(eps_abs, eps_rel * abs(res[0])):
                            txt = f'Projected masses integral problem ' \
                                  f'i_x={i}, i_y={j}, MGE comp. {k}, ' \
                                  f'PSF {psf_idx}: ' \
                                  f'err={res[1]}, but should be <= ' \
                                  f'{2 * max(eps_abs, eps_rel * abs(res[0]))}.'
                            self.logger.warning(txt)
                        grid[i, j] += res[0]
            # Note that the normalization by totalmass doesn't have a ml
            # factor (would be in surf_km). Even if it had, it would cancel
            # out (cf. surf_km->surcor in _integrand vs. surf_km->totalmass).
            totalmass = 2 * math.pi * np.sum(surf_km * qobs * sigobs_km ** 2)
            grid /= totalmass
            binfile = c.settings.io_settings['input_directory'] + kin.binfile
            if i_kin == 0:
                mass = self._bin_grid(np.ravel(grid, order='F'), binfile)
            else:
                mass = np.concatenate((mass,
                                       self._bin_grid(np.ravel(grid,order='F'),
                                                      binfile)))
        proj_mass = table.Table([mass], names=('proj_mass',))
        proj_mass.write(p_mass_fname, format='ascii.ecsv', overwrite=True)
        self.logger.info(f'Projected masses written to {p_mass_fname}.')
        return mass

    def _bin_grid(self, grid, binfile):
        """Read the binning file and assign the projected mass to the spatial
        bins.

        This is a helper method for get_projected_masses().

        Parameters
        ----------
        grid : numpy array, shape=(n_pixels,)
            The flattened mass contribution in the pixels.
        binfile : str
            The file name of the bins file.

        Returns
        -------
        numpy array, shape=(n_spatial_bins,)
            The projected mass for each spatial bin.

        Raises
        ------
        ValueError
            If there is a conflict in the bin and/or pixel numbers.
        """
        bininfo = [line.rstrip('\n').split()
                   for line in open(binfile) if line.lstrip(' ')[0] != '#']
        n_pixels = int(bininfo.pop(0)[0])
        bins = np.array([int(b) for c in bininfo for b in c])
        if len(bins) != n_pixels:
            txt = f'{binfile} says it has {n_pixels}, but it has {len(bins)}!'
            self.logger.error(txt)
            raise ValueError(txt)
        if max(bins) > grid.size:
            txt = f'Cannot store {max(bins)} bins in a grid sized {grid.size}!'
            self.logger.error(txt)
            raise ValueError(txt)
        n_bins = max(bins)
        txt = f'{binfile}: {n_pixels} pixels -> {n_bins} Voronoi bins.'
        self.logger.debug(txt)
        binned = np.zeros(n_bins + 1)  # Account for unused pixels
        for i in range(grid.size):
            binned[bins[i]] += grid[i]
        apermass = binned[1:]
        return apermass

    def _integrand(self, x, y0, y1, angle, qobs, sigobs_km,
                   surf_km, isotwist, psf_width, psf_weight):
        """The projected mass calulation's integrand at an independent variable
        value x

        This is a helper method for get_projected_masses().

        Parameters
        ----------
        x : float
            The integrand evaluation's independent variable value
        y0 , y1, angle, qobs, sigobs_km, surf_km, isotwist, psf_width,
        psf_weight : all floats
            Extra 'user parameters' for the integration

        Returns
        -------
        float
            The integrand evaluated at x
        """
        # compute new gaussian size after convolving with the psf:
        sb = math.sqrt(sigobs_km**2 + psf_width**2)
        qb = math.sqrt((sigobs_km*sigobs_km*qobs*qobs + psf_width*psf_width) /
                       (sigobs_km*sigobs_km + psf_width*psf_width))
        surcor = surf_km*qobs/qb*(sigobs_km**2) / (sb**2)

        # Angle is the angle form the PA to the X-axis ccw
        # To rotate the grid to the PA rotate of -angle ccw
        # Isotwist is the rotation from the PA to the major axis of the gaussian ccw
        # so this can just be added.
        alpha = -angle + isotwist

        # Glenn formula
        # k= sqrt( 1.0_dp+qb*qb+(1.0_dp-qb*qb)*cos(2.0_dp*alpha))
        # f0= 1.0_dp/(2.0_dp*k*qb*sb) * ( k*k*y0 - (1.0_dp-qb*qb)*sin(2.0_dp*alpha)*x)
        # f1= 1.0_dp/(2.0_dp*k*qb*sb) * ( k*k*y1 - (1.0_dp-qb*qb)*sin(2.0_dp*alpha)*x)
        # res = -psfweight * surcor * qb * sb * sqrt(pi_d) &
        #       / k * (derf(f0)-derf(f1)) * exp ( - (x*x)/(k*k*sb*sb))

        # Cappellari formula from ic1459 paper (2002) Appendix B3, formula (B6 and B7)
        p = math.sqrt(1 + qb*qb + (1.0 - qb*qb)*math.cos(2*alpha))
        f0 = 1.0/(2*p*qb*sb) * ((1.0 - qb*qb)*x*math.sin(2*alpha) - p*p*y0)
        f1 = 1.0/(2*p*qb*sb) * ((1.0 - qb*qb)*x*math.sin(2*alpha) - p*p*y1)

        return psf_weight*surcor*qb*sb*math.sqrt(math.pi)/p * \
               (special.erf(f0) - special.erf(f1))*math.exp(-(x/(p*sb))**2)

    def get_projected_masses_from_file(self, directory_noml):
        """read mge projected masses from ``mass_aper.dat``

        Parameters
        ----------
        directory_noml : string
            name of model directory exclusing the ``ml/`` extension

        Returns
        -------
        array
            array of aperture masses of the MGE

        """
        fname = f'{directory_noml}datfil/mass_aper.dat'
        aperture_masses = np.loadtxt(fname, skiprows=1)
        # remove first column (aperture index)
        aperture_masses = aperture_masses[:,1]
        return aperture_masses

    def get_intrinsic_masses(self, parset, grid):
        # TODO: reimplement
        # calculate the mass of the mge in observed 3D grid given the
        # parameter set containing intrinsic axis ratios (p, q, u)
        # for now, use legacy implementation below which reads from file
        pass

    def get_intrinsic_masses_from_file(self, directory_noml):
        """read mge intrinsic masses from ``mass_qgrid.dat``

        Parameters
        ----------
        directory_noml : string
            name of model directory exclusing the ``ml/`` extension

        Returns
        -------
        array
            3D intrinsic_masses masses of the MGE in a polar grid with sizes
            (n_r, n_theta, n_phi) which are defined in the config file.
            Their defaults are (6, 6, 10)

        """
        fname = f'{directory_noml}datfil/mass_qgrid.dat'
        shape = np.loadtxt(fname, max_rows=1, dtype=int)
        intrinsic_masses = np.loadtxt(fname, skiprows=1)
        intrinsic_masses = np.reshape(intrinsic_masses, shape)
        return intrinsic_masses

    def __add__(self,other):
        """Concatenate two MGEs, preserving row order.

        The input_directory and filename attributes are inherited from the
        first MGE.

        Parameters
        ----------
        other : Object of type MGE
            the MGE to add to this one

        Returns
        -------
        new_mge : Object of type MGE

        """
        mge1_data = self.data
        mge2_data = other.data
        mge1_data['row_merge_ID'] = list(range(1,len(mge1_data)+1))
        mge2_data['row_merge_ID'] = list(range(len(mge1_data)+1,
                                               len(mge1_data)+len(mge2_data)+1))

        new_data = table.join(mge1_data, mge2_data, join_type='outer')
        new_data.sort('row_merge_ID')
        new_data.remove_columns('row_merge_ID')

        new_mge = MGE(input_directory=self.input_directory,
                      datafile=self.datafile,
                      config=self.config)
        new_mge.data = new_data

        return new_mge



# end
