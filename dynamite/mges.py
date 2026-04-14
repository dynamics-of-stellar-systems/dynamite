import logging
import math
import numpy as np
from scipy import integrate
from scipy import special
from astropy import table
from pathos.multiprocessing import Pool
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

    def get_projected_masses(self, parset, apertures):
        # TODO:
        # calculate the mass of the mge in observed 2D apertures given the
        # parameter set containing intrinsic axis ratios (p, q, u)
        # for now, use legacy implementation below which reads from file
        pass

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
                             parallel=True):
        """Calculate the mass of the mge in the intrinsic grid.

        Calculate the mass of the mge in the intrinsic grid using observed
        quantities and the viewing angles. The intrinsic masses in the
        intrinsic grid and radius breaks are written to mass_radmass.ecsv and
        mass_qgrid.ecsv in the model's datfil/ directory, respectively.
        If called again, this method will read previously calculated intrinsic
        masses from the respective files unless use_cache=False.

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
                3D intrinsic_masses masses of the MGE in a polar grid of size
                (quad_nph, quad_nth, quad_nr) as defined in the config file

        """
        dir = model.directory_noml + 'datfil/'
        settings = self.config.settings.orblib_settings
        quad_nr = settings['quad_nr']
        quad_nth = settings['quad_nth']
        quad_nph = settings['quad_nph']

        # calculate the mass of the mge in observed 3D grid given the
        # parameter set containing intrinsic axis ratios (p, q, u)
        # EXPERIMENTAL!
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
                a = math.exp(-t*t/(2*sigintr_km[gn]**2)*(gx*gx + (gy*gy)/d + (gz*gz)/e))
                # dF/dx
                return -gx/sigintr_km[gn]**2*t*t*a/math.sqrt(d*e)

            def ayfunc(t, gx, gy, gz, gn):
                d = 1 - (1 - pintr[gn]*pintr[gn])*t*t
                e = 1 - (1 - qintr[gn]*qintr[gn])*t*t

                # Integral part of formula 12 of Cappellari 2002.
                a = math.exp(-t*t/(2*sigintr_km[gn]**2)*(gx*gx + (gy*gy)/d + (gz*gz)/e))
                # dF/dy
                return -gy/sigintr_km[gn]**2*t*t/d*a/math.sqrt(d*e)

            def azfunc(t, gx, gy, gz, gn):
                # integral part of formula 12 of Cappellari 2002.
                d = 1 - (1 - pintr[gn]*pintr[gn])*t*t
                e = 1 - (1 - qintr[gn]*qintr[gn])*t*t

                # Integral part of formula 12 of Cappellari 2002.
                a = math.exp(-t*t/(2*sigintr_km[gn]**2)*(gx*gx + (gy*gy)/d + (gz*gz)/e))
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
            theta_view, _, _ = self._get_viewing_angles_rad(model)
            disk_data = self.data[len_mge_bulge:]
            # Dispersion in km
            sigobs_km_d = disk_data['sigma'].data * arcsec_to_km
            # Surface brightness in L_sun/km^2 (we don't multiply by ml here)
            surf_km_d = disk_data['I'].data / constants.PARSEC_KM ** 2
            # Observed flattening
            qobs_d = disk_data['q'].data
            # Offset psi viewing angle in radians
            psi_obs_d = self.data['PA_twist'].data * math.pi / 180
            psi_obs_d +=  90  # psi_view_disk = pi/2
            qintr_d = (qobs_d ** 2 -
                       math.cos(theta_view) ** 2) / math.sin(theta_view) ** 2
            pintr_d = 0.9999999999 * np.ones_like(qintr_d)
            # triaxpar_d = np.zeros_like(qintr_d)
            sigintr_km_d = sigobs_km_d
            if any(qintr_d < 0):
                txt = 'q^2 is below 0 (in disk).'
                self.logger.error(txt)
                raise ValueError(txt)
            qintr_d = np.sqrt(qintr_d)
            if any(qintr_d > pintr_d):
                txt = 'q>p in disk.'
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
        # Total mass of the galaxy (except for ml parameter)
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
        return radmass, np.reshape(quad_grid,
                                   newshape=(quad_nph, quad_nth, quad_nr))

    def _get_intrinsic(self, model, mge_data, logtxt=''):
        arcsec_to_km = constants.ARC_KM(self.config.system.distMPc)
        theta_view, psi_view, phi_view = self._get_viewing_angles_rad(model)
        secth = 1 / math.cos(theta_view)
        cotph = 1 / math.tan(phi_view)
        # Dispersion in km
        sigobs_km = mge_data['sigma'].data * arcsec_to_km
        # Surface brightness in L_sun/km^2 (we don't multiply by ml here)
        surf_km = mge_data['I'].data / constants.PARSEC_KM ** 2
        # Observed flattening
        qobs = mge_data['q'].data
        # Offset psi viewing angle in radians
        psi_obs = mge_data['PA_twist'].data * math.pi / 180
        psi_obs += psi_view
        delp = 1 - qobs ** 2
        nom1minq2 = delp * \
            (2 * np.cos(2 * psi_obs) +
             np.sin(2 * psi_obs) *
                (secth * cotph - math.cos(theta_view) * math.tan(phi_view)))
        nomp2minq2 = \
            delp * (2 * np.cos(2 * psi_obs) +
                    np.sin(2 * psi_obs) * (math.cos(theta_view) * cotph -
                                           secth * math.tan(phi_view)))
        denom = \
            2 * math.sin(theta_view) ** 2 * (delp * np.cos(psi_obs) *
                (np.cos(psi_obs) + secth * cotph * np.sin(psi_obs)) - 1)
        # These are temporary values of the squared intrinsic axial
        # ratios p^2 and q^2
        qintr = (1 - nom1minq2 / denom)
        pintr = (qintr + nomp2minq2 / denom)
        if any(qintr < 0) or any(pintr < 0):
            txt = f"p^2 or q^2 is below 0{logtxt}."
            self.logger.error(txt)
            raise ValueError(txt)
        # intrinsic axial ratios p and q
        qintr = np.sqrt(qintr)
        pintr = np.sqrt(pintr)
        self.logger.debug(f'Middle axis ratio{logtxt} p={pintr}, '
                          f'minor axis ratio{logtxt} q={qintr}.')
        if any(qintr > pintr):
            txt = f"q > p{logtxt}."
            self.logger.error(txt)
            raise ValueError(txt)
        if any(pintr > 1):
            txt = f"p > 1{logtxt}."
            self.logger.error(txt)
            raise ValueError(txt)
        # intrinsic sigma (Cappellari 2002 eq 9.)
        sigintr_km = sigobs_km * \
            np.sqrt(qobs / np.sqrt((pintr * math.cos(theta_view)) ** 2 +
                (qintr * math.sin(theta_view)) ** 2 *
                ((pintr * math.cos(phi_view)) ** 2 + math.sin(phi_view) ** 2)))
        # triaxiality parameter T = (1-p^2)/(1-q^2)
        triaxpar = (1 - pintr ** 2) / (1 - qintr ** 2)
        if any(triaxpar < 0) or any(triaxpar > 1):
            txt = f'No triaxial deprojection possible{logtxt}!'
            self.logger.error(txt)
            raise ValueError(txt)
        self.logger.debug(f'Triaxiality parameters{logtxt}: {triaxpar}')
        return pintr, qintr, sigintr_km, surf_km, qobs, sigobs_km

    def _get_viewing_angles_rad(self, model):
        """Return the model's visible component's viewing angles.

        Parameters
        ----------
        model : a ``dyn.model.Model`` object

        Returns
        -------
        Tuple (theta_view, psi_view, phi_view)
            The viewing angels are returned in radians.
        """
        if self.config.system.is_bar_disk_system():
            stars = self.config.system.get_unique_bar_component()
        else:
            stars = self.config.system.get_unique_triaxial_visible_component()
        # used to derive the viewing angles
        parset = model.parset
        if self.config.system.is_bar_disk_system():
            if self.config.system.is_bar_disk_system_with_angles():
                theta_view = parset[f'theta-{stars.name}']
                phi_view = parset[f'phi-{stars.name}']
                psi_view = parset[f'psi-{stars.name}']
            else:
                q = parset[f'q-{stars.name}']
                p = parset[f'p-{stars.name}']
                u = parset[f'u-{stars.name}']
                qdisk = parset[f'qdisk-{stars.name}']
                theta_view, psi_view, phi_view = \
                    stars.triax_pqu2tpp_bar(p, q, u, qdisk)
                phi_view = -phi_view ## FIX ME
        else:
            q = parset[f'q-{stars.name}']
            p = parset[f'p-{stars.name}']
            u = parset[f'u-{stars.name}']
            theta_view, psi_view, phi_view = stars.triax_pqu2tpp(p, q, u)
        return theta_view * math.pi / 180, psi_view * math.pi / 180, \
               phi_view * math.pi / 180

    def _intrin_spher_grid(self, low, up, P, Q, sigma, r0, r1, rho0):
        ''' calculate multi-dimensional integrals'''
        epsabs = 1.49e-8
        epsrel = 1.49e-8
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
            # Be generous with the error...
            max_err = 20 * max(epsabs, epsrel * abs(res))
            if abserr > max_err:
                txt = f'Intrinsic masses integral problem: err={abserr}, ' \
                      f'but should be <= {max_err}.'
                self.logger.warning(txt)
            res = intgsign * res
        return res

    def _intrin_spher_grid_func(self, x, y, P, Q, sigma, r0, r1, rho0):
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
                self.logger.debug(f'{i=}, {i_gauss=}, {res/total_mass*100=}')
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
        # Calculate like in qgrid_setup (instead of reading the orblib-file):
        # quad_nr, quad_nth, quad_nph, quad_lr, quad_lth, quad_lph
        #START
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
        #END

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
                          f'{sum(quad_grid) / total_mass * 100}.')
        quad_grid /= total_mass
        quad_grid = np.ravel(quad_grid, order='F')  # retain legacy format

        self.logger.debug(f'quad_grid: {quad_grid}')
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
            name of model directory exclusing the ``ml/`` extension

        Returns
        -------
        array
            3D intrinsic_masses masses of the MGE in a polar grid with sizes
            (n_r, n_theta, n_phi) which are defined in the config file.
            Their defaults are (6,6,10)

        """
        fname = f'{directory_noml}datfil/mass_qgrid.dat'
        shape = np.loadtxt(fname, max_rows=1, dtype=int)
        intrinsic_masses = np.loadtxt(fname, skiprows=1)
        intrinsic_masses = np.reshape(intrinsic_masses, shape)
        return intrinsic_masses

    def __add__(self,other):
        """Concatenate two MGEs, preserving row order.

        The input_directory and filename attributes are inherited from the first MGE.

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
        mge2_data['row_merge_ID'] = list(range(len(mge1_data)+1,len(mge1_data)+len(mge2_data)+1))

        new_data = table.join(mge1_data, mge2_data, join_type='outer')
        new_data.sort('row_merge_ID')
        new_data.remove_columns('row_merge_ID')

        new_mge = MGE(input_directory=self.input_directory,
                      datafile=self.datafile,
                      config=self.config)
        new_mge.data = new_data

        return new_mge



# end
