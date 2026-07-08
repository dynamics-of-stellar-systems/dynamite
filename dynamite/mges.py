import logging
import math
import numpy as np
import os
from scipy import integrate
from scipy import special
from astropy import table
from pathos.multiprocessing import Pool
from dynamite import data
from dynamite import constants
from dynamite.physical_system import TriaxialVisibleComponent as vis_comp

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

    def get_projected_masses(self,
                             use_cache=True,
                             nocalc=False,
                             parallel=True):
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
        nocalc : bool, optional
            If True, ignore the values of use_cache and parallel and return
            the projected masses from the visible component's mass_aper
            attribute or from disk. If the data does not exist, raise an
            exception instead of attempting to calculate. The default is False.
        parallel : bool, optional
            If True, then the mass integration will be done in `ncpus`
            parallel processes where `ncpus` is taken from the configuration's
            `multiprocessing_settings`. If False, the integration will not use
            multiprocessing. False is recommended if called from within a
            parallel process. The default is True.

        Returns
        -------
        numpy array, shape=(n_spatial_bins,)
            Aperture masses of the MGE

        Raises
        ------
        FileNotFoundError
            If nocalc=True, the visible component's mass_aper attribute is
            None, and the projected mass file doesn't exist.
        """
        c = self.config
        vis_comp = c.system.get_unique_triaxial_visible_component()
        if (use_cache or nocalc) and (vis_comp.mass_aper is not None):
            self.logger.info('Projected masses grabbed from vis. component.')
            return vis_comp.mass_aper  ########################################
        p_mass_fname = c.settings.io_settings['output_directory'] + \
                       constants.p_masses_file
        if (use_cache or nocalc) and os.path.isfile(p_mass_fname):
            proj_mass = table.Table.read(p_mass_fname, format='ascii')
            self.logger.info(f'Projected masses read from {p_mass_fname}.')
            vis_comp.mass_aper = np.array(proj_mass['proj_mass'])
            return vis_comp.mass_aper  ########################################
        if nocalc:
            txt = "Unexpected: projected masses neither in visible " \
                  f"component's mass_aper attribute nor in folder {dir}."
            self.logger.error(txt)
            raise FileNotFoundError(txt)

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

            def _integrate(ij):
                i, j = ij
                out = 0.
                for k in range(len(self.data)):
                    for psf_idx in range(len(psf_width)):
                        res = integrate.quad(self._integrand,
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
                        out += res[0]
                return out

            grid = np.zeros((n_x, n_y))
            ij = [(i, j) for i in range(n_x) for j in range(n_y)]
            if parallel:
                with Pool(c.settings.multiprocessing_settings['ncpus']) as p:
                    output = p.map(_integrate, ij)
            else:
                output = [_integrate(ij_tuple) for ij_tuple in ij]
            for ij_idx, ij_tuple in enumerate(ij):
                grid[ij_tuple] = output[ij_idx]
            # Note that the normalization by totalmass doesn't have an ml
            # factor (would be in surf_km). In Fortran it has, but it cancels
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
        sb = math.sqrt(sigobs_km ** 2 + psf_width ** 2)
        qb = math.sqrt((sigobs_km*sigobs_km*qobs*qobs + psf_width*psf_width) /
                       (sigobs_km*sigobs_km + psf_width*psf_width))
        surcor = surf_km * qobs / qb * sigobs_km ** 2 / (sb ** 2)

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
        p = math.sqrt(1 + qb*qb + (1 - qb*qb)*math.cos(2*alpha))
        f0 = 1/(2*p*qb*sb) * ((1.0 - qb*qb)*x*math.sin(2*alpha) - p*p*y0)
        f1 = 1/(2*p*qb*sb) * ((1.0 - qb*qb)*x*math.sin(2*alpha) - p*p*y1)

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

    def get_intrinsic_masses(self,
                             model,
                             len_mge_bulge=None,
                             use_cache=True,
                             nocalc=False,
                             parallel=True):
        """Calculate the mass of the mge in the intrinsic grid.

        Calculate the mass of the mge in the intrinsic grid using observation
        data and the viewing angles. The intrinsic masses in the intrinsic grid
        and radius breaks are written to mass_qgrid.ecsv and mass_radmass.ecsv
        in the model's datfil/ directory, respectively. If called again, this
        method will read previously calculated intrinsic masses from the
        respective files unless use_cache=False.

        Parameters
        ----------
        model : a ``dynamite.model.Model`` object
        len_mge_bulge : int or None
            Bar systems have two mges which are concatenated in the .data
            attribute. The first len_mge_bulge records correspond to the bulge
            mge and the remaining to the disk mge. This parameter is only used
            for bar systems and must be None for non-bar systems.
        use_cache : bool, optional
            If False, the intrinsic masses will be recalculated.
            If True, check for intrinsic masses already existing on disk and
            if yes, read and return that data. The default is True.
        nocalc : bool, optional
            If True, ignore the values of use_cache and parallel and return
            the intrinsic masses already existing on disk. If the files do not
            exist, raise an exception instead of attempting to calculate.
            The default is False.
        parallel : bool, optional
            If True, then the mass integration will be done in `ncpus`
            parallel processes where `ncpus` is taken from the configuration's
            `multiprocessing_settings`. If False, the integration will not use
            multiprocessing. False is recommended if called from within a
            parallel process. The default is True.

        Returns
        -------
        Tuple (radmass, quad_grid)
            radmass : np.array, shape=(n_r,)
                mass inside the radial shells
            quad_grid : np.array, shape=(quad_nph, quad_nth, quad_nr)
                3D intrinsic masses of the MGE in a polar grid of size
                (quad_nph, quad_nth, quad_nr) as defined in the config file

        Raises
        ------
        FileNotFoundError
            If nocalc=True and one or both of the intrinsic mass file(s) don't
            exist.
        """
        settings = self.config.settings.orblib_settings
        quad_nr = settings['quad_nr']
        quad_nth = settings['quad_nth']
        quad_nph = settings['quad_nph']
        dir = model.directory_noml + 'datfil/'
        if use_cache or nocalc:
            mr_file = dir + 'mass_radmass.ecsv'
            mq_file = dir + 'mass_qgrid.ecsv'
            mr_exists = os.path.isfile(mr_file)
            mq_exists = os.path.isfile(mq_file)
            if mr_exists and mq_exists:
                radmass = table.Table.read(mr_file,
                                           format='ascii')['mass_radmass'].data
                quad_grid = table.Table.read(mq_file,
                                             format='ascii')['mass_qgrid'].data
                self.logger.info(f'Intrinsic masses read from {mr_file} '
                                 f'and {mq_file}, respectively.')
                return radmass, \
                       np.reshape(quad_grid, (quad_nph, quad_nth, quad_nr))
                ###############################################################
            elif mr_exists or mq_exists:
                self.logger.info(f'Strange, only one of {mr_file} and '
                    f'{mq_file} exists, will pretend neither exists.')
            if nocalc:
                txt=f'Unexpected: intrinsic mass file(s) inexistent in {dir}.'
                self.logger.error(txt)
                raise FileNotFoundError(txt)
        self.logger.info('Calculating intrinsic masses...')

        def potin(n, x, y, z):
            x2 = x*x
            y2 = y*y
            z2 = z*z

            A12 = -(A1[n] - A2[n])/(1 - pintr[n]**2)
            A23 = -(A2[n] - A3[n])/(pintr[n]**2 - qintr[n]**2)
            A31 = -(A3[n] - A1[n])/(qintr[n]**2 - 1)

            A11 = 1/3 * (2 - A12 - A31)
            A22 = 1/3 * (2/(pintr[n]**2) - A23 - A12)
            A33 = 1/3 * (2/(qintr[n]**2) - A31 - A23)

            O1 = -1/(2*sigintr_km[n]**2)*(A1[n]*x2 + A2[n]*y2 + A3[n]*z2)
            O2 = 1/(8*sigintr_km[n]**4)* \
                (A11*x2*x2 + A22*y2*y2 + A33*z2*z2 +
                2*A12*x2*y2 + 2*A23*y2*z2 + 2*A31*z2*x2)

            return V0[n]/math.sqrt(1 - qintr[n]**2)*(F[n] + O1 + O2)

        def potmid(n, x, y, z):

            def potfunc(t, gx, gy, gz, gn):
                d = 1 - ((1 - pintr[gn]*pintr[gn])*t*t)
                e = 1 - ((1 - qintr[gn]*qintr[gn])*t*t)
                # Integral part of formula 12 of Cappellari 2002.
                a = math.exp(-t*t/(2*sigintr_km[gn]**2)*
                            (gx*gx + (gy*gy)/d + (gz*gz)/e))
                return a/math.sqrt(d*e)

            epsabs = 1.49e-08
            epsrel = 1.49e-08
            limit = 50 # 20
            res = integrate.quad(potfunc,
                                 0,
                                 1,
                                 args=(x, y, z, n),
                                 epsabs=epsabs,
                                 epsrel=epsrel,
                                 limit=limit)
            if res[1] > 2 * max(epsabs, epsrel * abs(res[0])):
                txt = f'Intrinsic masses potmid integral problem: ' \
                      f'err={res[1]}, but should be <= ' \
                      f'{2 * max(epsabs, epsrel * abs(res[0]))}.'
                self.logger.warning(txt)

            return V0[n] * res[0]

        def accin(n, x, y, z):
            x2 = x*x
            y2 = y*y
            z2 = z*z

            A12 = -(A1[n] - A2[n])/(1 - pintr[n]**2)
            A23 = -(A2[n] - A3[n])/(pintr[n]**2 - qintr[n]**2)
            A31 = -(A3[n] - A1[n])/(qintr[n]**2 - 1)

            A11 = 1/3 * (2 - A12 - A31)
            A22 = 1/3 * (2/(pintr[n]**2) - A23 - A12)
            A33 = 1/3 * (2/(qintr[n]**2) - A31 - A23)

            vx = -V0[n]/math.sqrt(1 - qintr[n]**2)*x/sigintr_km[n]**2* \
                (A1[n] - 1/(2*sigintr_km[n]**2)*(A11*x2 + A12*y2 + A31*z2))
            vy = -V0[n]/math.sqrt(1 - qintr[n]**2)*y/sigintr_km[n]**2* \
                (A2[n] - 1/(2*sigintr_km[n]**2)*(A12*x2 + A22*y2 + A23*z2))
            vz = -V0[n]/math.sqrt(1 - qintr[n]**2)*z/sigintr_km[n]**2* \
                (A3[n] - 1/(2*sigintr_km[n]**2)*(A31*x2 + A23*y2 + A33*z2))

            return vx, vy, vz

        def accmid(n, x, y, z):

            def axfunc(t, gx, gy, gz, gn):
                d = 1 - (1 - pintr[gn]*pintr[gn])*t*t
                e = 1 - (1 - qintr[gn]*qintr[gn])*t*t
                # Integral part of formula 12 of Cappellari 2002.
                a = math.exp(-t*t/(2*sigintr_km[gn]**2)*
                             (gx*gx + (gy*gy)/d + (gz*gz)/e))
                # dF/dx
                return -gx/sigintr_km[gn]**2*t*t*a/math.sqrt(d*e)

            def ayfunc(t, gx, gy, gz, gn):
                d = 1 - (1 - pintr[gn]*pintr[gn])*t*t
                e = 1 - (1 - qintr[gn]*qintr[gn])*t*t

                # Integral part of formula 12 of Cappellari 2002.
                a = math.exp(-t*t/(2*sigintr_km[gn]**2)*
                             (gx*gx + (gy*gy)/d + (gz*gz)/e))
                # dF/dy
                return -gy/sigintr_km[gn]**2*t*t/d*a/math.sqrt(d*e)

            def azfunc(t, gx, gy, gz, gn):
                # integral part of formula 12 of Cappellari 2002.
                d = 1 - (1 - pintr[gn]*pintr[gn])*t*t
                e = 1 - (1 - qintr[gn]*qintr[gn])*t*t

                # Integral part of formula 12 of Cappellari 2002.
                a = math.exp(-t*t/(2*sigintr_km[gn]**2)*
                             (gx*gx + (gy*gy)/d + (gz*gz)/e))
                # dF/dz
                return -gz/sigintr_km[gn]**2*t*t/e*a/math.sqrt(d*e)

            epsabs = 1.49e-08
            epsrel = 1.49e-08
            limit = 50
            res = integrate.quad(axfunc,
                                 0.,
                                 1.,
                                 args=(x, y, z, n),
                                 epsabs=epsabs,
                                 epsrel=epsrel,
                                 limit=limit)
            if res[1] > 2 * max(epsabs, epsrel * abs(res[0])):
                txt = f'Intrinsic masses vx integral problem: ' \
                      f'err={res[1]}, but should be <= ' \
                      f'{2 * max(epsabs, epsrel * abs(res[0]))}.'
                self.logger.warning(txt)
            vx = V0[n] * res[0]
            res = integrate.quad(ayfunc,
                                 0.,
                                 1.,
                                 args=(x, y, z, n),
                                 epsabs=epsabs,
                                 epsrel=epsrel,
                                 limit=limit)
            if res[1] > 2 * max(epsabs, epsrel * abs(res[0])):
                txt = f'Intrinsic masses vy integral problem: ' \
                      f'err={res[1]}, but should be <= ' \
                      f'{2 * max(epsabs, epsrel * abs(res[0]))}.'
                self.logger.warning(txt)
            vy = V0[n] * res[0]
            res = integrate.quad(azfunc,
                                 0.,
                                 1.,
                                 args=(x, y, z, n),
                                 epsabs=epsabs,
                                 epsrel=epsrel,
                                 limit=limit)
            if res[1] > 2 * max(epsabs, epsrel * abs(res[0])):
                txt = f'Intrinsic masses vy integral problem: ' \
                      f'err={res[1]}, but should be <= ' \
                      f'{2 * max(epsabs, epsrel * abs(res[0]))}.'
                self.logger.warning(txt)
            vz = V0[n] * res[0]

            return vx, vy, vz

        if self.config.system.is_bar_disk_system():
            if type(len_mge_bulge) is not int:
                txt = 'len_mge_bulge must be an integer for bar systems.'
                self.logger.error(txt)
                raise ValueError(txt)
            # DISK
            arcsec_to_km = constants.ARC_KM(self.config.system.distMPc)
            theta_view, _, _ = self._get_viewing_angles_deg(model)
            theta_view = np.deg2rad(theta_view)
            disk_data = self.data[len_mge_bulge:]
            # Dispersion in km
            sigobs_km_d = disk_data['sigma'].data * arcsec_to_km
            # Surface brightness in L_sun/km^2 (we don't multiply by ml here)
            surf_km_d = disk_data['I'].data / constants.PARSEC_KM ** 2
            # Observed flattening
            qobs_d = disk_data['q'].data
            # # Offset psi viewing angle in radians
            # psi_obs_d = self.data['PA_twist'].data
            # psi_obs_d +=  90  # psi_view_disk = pi/2
            # psi_obs_d = np.deg2rad(psi_obs_d)
            qintr_d = (qobs_d ** 2 -
                       math.cos(theta_view) ** 2) / math.sin(theta_view) ** 2
            pintr_d = 0.9999999999 * np.ones_like(qintr_d)
            # triaxpar_d = np.zeros_like(qintr_d)
            sigintr_km_d = sigobs_km_d
            if any(qintr_d < 0):
                txt = 'qintr^2 is below 0 (in disk).'
                self.logger.error(txt)
                raise ValueError(txt)
            qintr_d = np.sqrt(qintr_d)
            if any(qintr_d > pintr_d):
                txt = 'qintr > pintr in disk.'
                self.logger.error(txt)
                raise ValueError(txt)
            # BULGE
            intr_props = self._get_intrinsic(model,
                                             self.data[:len_mge_bulge],
                                             logtxt=' (bulge)')
            pintr_b, qintr_b, sigintr_km_b, surf_km_b, qobs_b, sigobs_km_b = \
                intr_props
            # Combine all the gaussians together (disk + bulge). Note that the
            # sequence is now disk first, then bulge like in LegacyFortran.
            surf_km = np.concatenate((surf_km_d, surf_km_b))
            qobs = np.concatenate((qobs_d, qobs_b))
            sigobs_km = np.concatenate((sigobs_km_d, sigobs_km_b))
            # psi_obs = np.concatenate((psi_obs_d, psi_obs_b))
            qintr = np.concatenate((qintr_d, qintr_b))
            pintr = np.concatenate((pintr_d, pintr_b))
            sigintr_km = np.concatenate((sigintr_km_d, sigintr_km_b))
            # triaxpar = np.concatenate((triaxpar_d, triaxpar_b))
        else:
            if len_mge_bulge is not None:
                txt = 'len_mge_bulge must be None for non-bar systems.'
                self.logger.error(txt)
                raise ValueError(txt)
            intr_props = self._get_intrinsic(model,
                                             self.data)
            pintr, qintr, sigintr_km, surf_km, qobs, sigobs_km = intr_props
        self.logger.debug('Unitless length of the projected major axis: u = '
                          f'{sigobs_km / sigintr_km}')
        # density factor
        dens = surf_km * qobs * sigobs_km ** 2 / \
                (math.sqrt(math.tau) * pintr * qintr * sigintr_km ** 3)
        # Integration Constant
        V0 = 4 * math.pi * constants.GRAV_CONST_KM * \
              sigintr_km ** 2 * pintr * qintr * dens
        # Total mass of the galaxy (except for ml i.e., assuming ml=1)
        total_mass = math.tau * sum(surf_km * qobs * sigobs_km ** 2)
        # Compute the constants for the inner approximation
        k = np.sqrt((1 - pintr ** 2) / (1 - qintr ** 2))
        F = special.ellipkinc(np.arccos(qintr), k ** 2)
        E = special.ellipeinc(np.arccos(qintr), k ** 2)
        A1 = (F - E) / (1 - pintr ** 2)
        A2 = ((1 - qintr ** 2) * E - (pintr ** 2 - qintr ** 2) * F -
              qintr / pintr * ( 1 - pintr ** 2) * np.sqrt(1 - qintr ** 2)) \
             / ((1 - pintr ** 2) * (pintr ** 2 - qintr ** 2))
        A3 = (pintr / qintr * np.sqrt(1 - qintr ** 2) - E) \
             / (pintr ** 2 - qintr ** 2)
        if any(abs(np.sqrt(1 - qintr ** 2) / (pintr * qintr) - A1 - A2 - A3) \
           > 1e-6):
            txt = "Failure to properly compute F, A1, A2 and/or A3. "
            txt += f"{F=}, {E=}, {A1=}, {A2=}, {A3=}, "
            txt += f"{abs((np.sqrt(1-qintr**2)/(pintr*qintr))-A1-A2-A3)=}."
            self.logger.error(txt)
            raise ValueError(txt)

        self.logger.info('Testing the accuracy of the approximation regimes.')
        inner_approx = 0.0001
        outer_approx = 300
        for i in range(len(self.data)):
            ax = potin(i, inner_approx * sigintr_km[i], 1, 1)
            ix = potmid(i, inner_approx * sigintr_km[i], 1, 1)
            if abs((ix - ax) / ix) > 1.0e-4:
                 txt = f"Failed test 1: {ix} != {ax}"
                 self.logger.error(txt)
                 raise ValueError(txt)
            ax = potin(i, 1, 1, inner_approx * sigintr_km[i])
            ix = potmid(i, 1, 1, inner_approx * sigintr_km[i])
            if abs((ix - ax) / ix) > 1.0e-4:
                 txt = f"Failed test 2: {ix} != {ax}"
                 self.logger.error(txt)
                 raise ValueError(txt)
            ax = math.sqrt(math.pi / 2.) * V0[i] / outer_approx
            ix = potmid(i, outer_approx * sigintr_km[i], 0, 0)
            if abs((ix - ax) / ix) > 1.0e-4:
                 txt = f"Failed test 3: {ix} != {ax}"
                 self.logger.error(txt)
                 raise ValueError(txt)
            ax = math.sqrt(math.pi / 2.) * V0[i] / outer_approx
            ix = potmid(i, 1, 1, outer_approx * sigintr_km[i])
            if abs((ix - ax) / ix) > 1.0e-4:
                 txt = f"Failed test 4: {ix} != {ax}"
                 self.logger.error(txt)
                 raise ValueError(txt)
            ax, ay, az = accin(i,
                               inner_approx * sigintr_km[i] * 0.95,
                               0.2 * inner_approx * sigintr_km[i],
                               0.2 * inner_approx * sigintr_km[i])
            ix, iy, iz = accmid(i,
                                inner_approx * sigintr_km[i]*0.95,
                                0.2 * inner_approx * sigintr_km[i],
                                0.2 * inner_approx * sigintr_km[i])
            if math.sqrt((ix - ax) ** 2 + (iy - ay) ** 2 + (iz - az) ** 2) / \
               math.sqrt(ix ** 2 + iy ** 2 + iz ** 2) > 1.0e-3:
                txt = "Failed test 5: " \
                      f"{math.sqrt((ix-ax)**2 + (iy-ay)**2 + (iz-az)**2)} " \
                      f"!= {math.sqrt(ix**2+iy**2+iz**2)}"
                self.logger.error(txt)
                raise ValueError(txt)
            ax, ay, az = accin(i,
                               0.2 * inner_approx * sigintr_km[i],
                               0.2 * inner_approx * sigintr_km[i],
                               0.95 * inner_approx * sigintr_km[i])
            ix, iy, iz = accmid(i,
                                0.2 * inner_approx * sigintr_km[i],
                                0.2 * inner_approx * sigintr_km[i],
                                0.95 * inner_approx * sigintr_km[i])
            if math.sqrt((ix - ax) ** 2 + (iy - ay) ** 2 + (iz - az) ** 2) / \
               math.sqrt(ix ** 2 + iy ** 2 + iz ** 2) > 1.0e-2:
                txt = "Failed test 6: " \
                      f"{math.sqrt((ix-ax)**2 + (iy-ay)**2 + (iz-az)**2)} " \
                      f"!= {math.sqrt(ix**2+iy**2+iz**2)}"
                self.logger.error(txt)
                raise ValueError(txt)
        self.logger.info('Integrating intrinsic masses: radmass...')
        radmass = self._intrin_radii(total_mass=total_mass,
                                     pintr=pintr,
                                     qintr=qintr,
                                     sigintr_km=sigintr_km,
                                     dens=dens,
                                     dir=dir)
        self.logger.info('...and qgrid...')
        quad_grid = self._intrin_spher(total_mass=total_mass,
                                      pintr=pintr,
                                      qintr=qintr,
                                      sigintr_km=sigintr_km,
                                      dens=dens,
                                      dir=dir,
                                      parallel=parallel)
        return radmass, np.reshape(quad_grid, (quad_nph, quad_nth, quad_nr))

    def _get_intrinsic(self, model, mge_data, logtxt=''):
        """Helper function: calculate intrinsic properties

        Parameters
        ----------
        model : a ``dynamite.model.Model`` object
        mge_data : astropy table of length len_mge
            The mge data delivering dispersion, surface brightness, observed
            flattening, and the offset psi viewing angle for each mge component
        logtxt : str, optional
            The logtxt is added to logged messages and is intended to help
            keeping calls to this method from different places apart. The
            default is '' (the empty string).

        Returns
        -------
        tuple (pintr, qintr, sigintr_km, surf_km, qobs, sigobs_km)
            pintr and qintr are the intrinsic axial ratios, sigintr_km the
            intrinsic dispersion, surf_km the observed surface brightness, qobs
            the observed flattening, and sigobs_km the observed dispersion.
            Each tuple element is a numpy array of shape (len_mge,).

        Raises
        ------
        ValueError
            If pintr^2 < 0 or qintr^2 < 0, qintr > pintr, pintr > 1, or the
            triaxiality parameter is less than 0 or greater than 1
        """
        arcsec_to_km = constants.ARC_KM(self.config.system.distMPc)
        theta_view, psi_view, phi_view = self._get_viewing_angles_deg(model)
        # Dispersion in km
        sigobs_km = mge_data['sigma'].data * arcsec_to_km
        # Surface brightness in L_sun/km^2 (we don't multiply by ml here)
        surf_km = mge_data['I'].data / constants.PARSEC_KM ** 2
        # Observed flattening
        qobs = mge_data['q'].data
        # Offset psi viewing angle in degrees
        psi_obs = mge_data['PA_twist'].data
        pintr, qintr, _ = vis_comp.triax_tpp2pqu(theta=theta_view,
                                                 phi=phi_view,
                                                 psi=psi_view,
                                                 qobs_pot=qobs,
                                                 psi_off=psi_obs)
        self.logger.debug(f'Middle axis ratio{logtxt} {pintr=}, '
                          f'minor axis ratio{logtxt} {qintr=}.')
        # intrinsic sigma (Cappellari 2002 eq 9.)
        theta_view = np.deg2rad(theta_view)
        phi_view = np.deg2rad(phi_view)
        sigintr_km = sigobs_km * \
            np.sqrt(qobs / np.sqrt((pintr * math.cos(theta_view)) ** 2 +
                (qintr * math.sin(theta_view)) ** 2 *
                ((pintr * math.cos(phi_view)) ** 2 + math.sin(phi_view) ** 2)))
        return pintr, qintr, sigintr_km, surf_km, qobs, sigobs_km

    def _get_viewing_angles_deg(self, model):
        """Return the model's visible component's viewing angles.

        Parameters
        ----------
        model : a ``dyn.model.Model`` object

        Returns
        -------
        Tuple (theta_view, psi_view, phi_view)
            The viewing angels are returned in degrees.
        """
        stars = self.config.system.get_unique_triaxial_visible_component()
        # used to derive the viewing angles
        parset = model.parset
        if self.config.system.is_bar_disk_system():
            theta_view = parset[f'theta-{stars.name}']
            phi_view = parset[f'phi-{stars.name}']
            psi_view = parset[f'psi-{stars.name}']
        else:
            q = parset[f'q-{stars.name}']
            p = parset[f'p-{stars.name}']
            u = parset[f'u-{stars.name}']
            theta_view, psi_view, phi_view = stars.triax_pqu2tpp(p, q, u)
        return theta_view, psi_view, phi_view

    def _intrin_spher_grid(self, low, up, P, Q, sigma, r0, r1, rho0):
        """Calculate the intrinsic mass multi-dimensional integrals"""
        epsabs = 5e-7  # 1.49e-8
        epsrel = 5e-7  # 1.49e-8
        limit = 150
        # check integration limits
        intgsign = 1
        low, up = np.array(low), np.array(up)
        if any(low == up):
            self.logger.debug("lower and upper limit are equal")
            res = 0.
        else:
            for ii in range(2):
                if low[ii] > up[ii]:
                    self.logger.debug("changing order of integration")
                    low[ii], up[ii] = up[ii], low[ii]
                    intgsign = -intgsign
            res, abserr = integrate.nquad(self._intrin_spher_grid_func,
                                          ranges=[[low[0], up[0]],
                                                  [low[1], up[1]]],
                                          args=(P, Q, sigma, r0, r1, rho0),
                                          opts={'epsabs' : epsabs,
                                                'epsrel' : epsrel,
                                                'limit' : limit})
            # Be somewhat generous with the error...
            max_err = 2 * max(epsabs, epsrel * abs(res))
            if abserr > max_err:
                txt = f'Intrinsic masses integral problem: err={abserr}, ' \
                      f'but should be <= {max_err}.'
                self.logger.warning(txt)
            res = intgsign * res
        return res

    def _intrin_spher_grid_func(self, x, y, P, Q, sigma, r0, r1, rho0):
        """Integrand of the intrinsic mass multi-dimensional integrals"""
        sth = math.sin(x)
        cth = math.cos(x)
        sph = math.sin(y)
        cph = math.cos(y)
        c = math.sqrt(2.) * sigma / \
            math.sqrt((cph**2+(sph/P)**2)*sth*sth + (cth/Q)**2)
        term1 = r0*math.exp(-(r0/c)**2) - \
                0.5*c*math.sqrt(math.pi)*special.erf(r0/c)
        term2 = r1*math.exp(-(r1/c)**2) - \
                0.5*c*math.sqrt(math.pi)*special.erf(r1/c)
        answ = 0.5*rho0*c*c*(term1 - term2)*sth
        # the orbit library folds the 8 symmetries into one octant:
        return answ * 8

    def _intrin_radii(self, total_mass, pintr, qintr, sigintr_km, dens, dir):
        """Integrate the intrinsic mass inside the radial shells

        Parameters
        ----------
        total_mass : float
            Total mass of the system
        pintr : numpy array of shape (len_mge,)
            Intrinsic axial ratio p
        qintr : numpy array of shape (len_mge,)
            Intrinsic axial ratio q
        sigintr_km : numpy array of shape (len_mge,)
            Intrinsic dispersion
        dens : numpy array of shape (len_mge,)
            Density factor
        dir : str
            Destination directory for the output file

        Returns
        -------
        np.array, shape=(n_r,)
            Intrinsic mass inside the radial shells
        """
        settings = self.config.settings.orblib_settings
        arcsec_to_km = constants.ARC_KM(self.config.system.distMPc)
        nr = settings['nE']
        rlogmin = settings['logrmin'] + math.log10(arcsec_to_km)
        rlogmax = settings['logrmax'] + math.log10(arcsec_to_km)

        radmass = np.zeros(nr)
        quad_lr = np.zeros(nr + 1)
        for i in range(nr):
            quad_lr[i] = \
                10. ** (rlogmin + (rlogmax - rlogmin) * (i - 1) / (nr - 1))
        quad_lr[0] = 0.0
        quad_lr[nr] = 10 ** rlogmax * 100.0
        quad_lth = np.array([0, math.pi / 2])
        quad_lph = np.array([0, math.pi / 2])

        low = (quad_lth[0], quad_lph[0])
        up = (quad_lth[1], quad_lph[1])
        for i in range(nr):
            for i_gauss in range(len(self.data)):
                res = self._intrin_spher_grid(low,
                                              up,
                                              P=pintr[i_gauss],
                                              Q=qintr[i_gauss],
                                              sigma=sigintr_km[i_gauss],
                                              r0=quad_lr[i],
                                              r1=quad_lr[i + 1],
                                              rho0=dens[i_gauss])
                radmass[i] += res
        self.logger.debug('Percent of the mass inside the radial shells: '
                          f'{sum(radmass) / total_mass * 100}')
        radmass /= total_mass
        self.logger.debug(f'radmass: {radmass}')
        mass_radmass_filename = dir + 'mass_radmass.ecsv'
        radmass_table = table.Table([radmass], names=('mass_radmass',))
        radmass_table.write(mass_radmass_filename,
                            format='ascii.ecsv',
                            overwrite=True)
        self.logger.info('Intrinsic masses mass_radmass written to '
                         f'{mass_radmass_filename}.')
        return radmass

    def _intrin_spher(self, total_mass, pintr, qintr, sigintr_km, dens, dir,
                      parallel=True):
        """Integrate the intrinsic mass inside the polar grid

        Calculate the 3D intrinsic masses of the MGE in a polar grid of size
        (quad_nph, quad_nth, quad_nr) as defined in the config file.

        Parameters
        ----------
        total_mass : float
            Total mass of the system
        pintr : numpy array of shape (len_mge,)
            Intrinsic axial ratio p
        qintr : numpy array of shape (len_mge,)
            Intrinsic axial ratio q
        sigintr_km : numpy array of shape (len_mge,)
            Intrinsic dispersion
        dens : numpy array of shape (len_mge,)
            Density factor
        dir : str
            Destination directory for the output file
        parallel : bool, optional
            If True, then the mass integration will be done in `ncpus`
            parallel processes where `ncpus` is taken from the configuration's
            `multiprocessing_settings`. If False, the integration will not use
            multiprocessing. False is recommended if called from within a
            parallel process. The default is True.

        Returns
        -------
        quad_grid : 1d np.array, shape=(quad_nph * quad_nth * quad_nr,)
            3D intrinsic masses of the MGE in a polar grid of size
            (quad_nph, quad_nth, quad_nr).
        """
        # Calculate like in qgrid_setup (instead of reading the orblib-file):
        # quad_nr, quad_nth, quad_nph, quad_lr, quad_lth, quad_lph
        c = self.config
        settings = c.settings.orblib_settings
        arcsec_to_km = constants.ARC_KM(self.config.system.distMPc)
        sigobs_km = self.data['sigma'] * arcsec_to_km
        rlogmin = settings['logrmin'] + math.log10(arcsec_to_km)
        rlogmax = settings['logrmax'] + math.log10(arcsec_to_km)
        quad_nr = settings['quad_nr']
        quad_nth = settings['quad_nth']
        quad_nph = settings['quad_nph']
        self.logger.debug(f'Grid dimension: {quad_nr = }, '
                          f'{quad_nth = }, {quad_nph = }.')

        # Define a grid such that the boundaries define all possible bins.
        # This also means that there are N+1 boundaries for N bins.
        quad_lr = np.zeros(quad_nr + 1)
        quad_lr[0] = 0.
        for i in range(1, quad_nr):
            quad_lr[i] = 10**(rlogmin +
                (rlogmax - rlogmin + math.log10(0.5)) * i / quad_nr)
        quad_lr[quad_nr] = max(10 ** rlogmax * 100, max(sigobs_km) * 10)
        # Define the angular bins
        quad_lth = np.zeros(quad_nth + 1)
        quad_lth[0] = 0
        for i in range(1, quad_nth):
            quad_lth[i] = math.pi / 2 * i / quad_nth
        quad_lth[quad_nth] = math.pi / 2
        # Define the angular bins
        quad_lph = np.zeros(quad_nph + 1)
        for i in range(1, quad_nph):
            quad_lph[i] = math.pi / 2 * i / quad_nph
        quad_lph[quad_nph] = math.pi / 2

        def _integrate(ijk):
            i, j, k = ijk
            out = 0
            for i_gauss in range(len(self.data)):
                low = (quad_lth[j], quad_lph[k])
                up = (quad_lth[j + 1], quad_lph[k + 1])
                res = self._intrin_spher_grid(low,
                                              up,
                                              P=pintr[i_gauss],
                                              Q=qintr[i_gauss],
                                              sigma=sigintr_km[i_gauss],
                                              r0=quad_lr[i],
                                              r1=quad_lr[i + 1],
                                              rho0=dens[i_gauss])
                out += res
            return out

        quad_grid = np.zeros((quad_nph, quad_nth, quad_nr))
        ijk = [(i, j, k) for i in range(quad_nr)
                         for j in range(quad_nth)
                         for k in range(quad_nph)]
        if parallel:
            with Pool(c.settings.multiprocessing_settings['ncpus']) as p:
                output = p.map(_integrate, ijk)
        else:
            output = [_integrate(ijk_tuple) for ijk_tuple in ijk]
        for ijk_idx, ijk_tuple in enumerate(ijk):
            quad_grid[ijk_tuple[::-1]] = output[ijk_idx]

        self.logger.debug('Percent of the Mass inside the projected grid: '
                          f'{np.sum(quad_grid) / total_mass * 100}.')
        quad_grid /= total_mass
        quad_grid = np.ravel(quad_grid, order='F')  # retain legacy format

        mass_qgrid_filename = dir + 'mass_qgrid.ecsv'
        qgrid_table = table.Table([quad_grid], names=('mass_qgrid',))
        qgrid_table.meta = \
            {'quad_nph': quad_nph,
             'quad_nth': quad_nth,
             'quad_nr': quad_nr,
             'note': '3D array flattened in Fortran index order.'}
        qgrid_table.write(mass_qgrid_filename,
                          format='ascii.ecsv',
                          overwrite=True)
        self.logger.info('Intrinsic masses mass_qgrid written to '
                         f'{mass_qgrid_filename}.')
        return quad_grid

    def get_intrinsic_masses_from_file(self, directory_noml):
        """read mge intrinsic masses from ``mass_qgrid.dat``

        Parameters
        ----------
        directory_noml : string
            name of model directory excluding the ``ml/`` extension

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
