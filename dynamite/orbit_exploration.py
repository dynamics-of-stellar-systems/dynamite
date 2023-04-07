import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from plotbin.display_pixels import display_pixels
import astropy
import dynamite as dyn

class Decomposition:
    """
    Class for decomposition.

    Upon instatiating, the orbits are decomposed by the method
    ``decompose_orbits`` and the results stored in astropy table
    ``self.decomp``. The components' flux and moments (currently mean velocity
    and velocity dispersion only) are plotted by calling ``self.plot_decomp``
    which also writes the plotted data into the model directory.

    The methodology in this class has been contributed by Giulia Santucci.
    Please cite Santucci+22, ApJ 930, 2, 153 if used.

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object, mandatory
    model : a ``dyn.model.Model`` object, optional
        Determines which model is used.
        If model = None, the model corresponding to the minimum
        chisquare (so far) is used; the setting in the configuration
        file's parameter settings is used to determine which chisquare
        to consider. The default is None.
    kin_set : int, optional
        Determines which kinematic set to use.
        The value of this parameter is the index of the data
        set (e.g. kin_set=0 , kin_set=1). The default is 0.

    Raises
    ------
    ValueError
        if no config object is given or the kin_set does not exist.

    """
    def __init__(self, config=None, model=None, kin_set=0):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None:
            text = f'{__class__.__name__} needs configuration object, ' \
                   'None provided.'
            self.logger.error(text)
            raise ValueError(text)
        self.config = config
        if model is None:
            # Select the best model for decomposition
            best_model_idx = config.all_models.get_best_n_models_idx(n=1)[0]
            self.model = config.all_models.get_model_from_row(best_model_idx)
        stars = \
          config.system.get_component_from_class(
                                  dyn.physical_system.TriaxialVisibleComponent)
        n_kin = len(stars.kinematic_data)
        if kin_set >= n_kin:
            text = f'kin_set must be < {n_kin}, but it is {kin_set}'
            self.logger.error(text)
            raise ValueError(text)
        self.kin_set = kin_set
        self.logger.info(f'Performing decomposition for kin_set no {kin_set}: '
                         f'{stars.kinematic_data[kin_set].name}')
        # Get losvd_histograms and projected_masses
        orblib = self.model.get_orblib()
        orblib.read_losvd_histograms()
        self.losvd_histograms = orblib.losvd_histograms[self.kin_set]
        self.proj_mass = orblib.projected_masses[self.kin_set]
        self.logger.debug(f'{self.losvd_histograms.y.shape=}, '
                          f'{self.proj_mass.shape=}.')
        # Get orbit weights and store them in self.model.weights
        _ = self.model.get_weights(orblib)
        # Do the decomposition
        self.decomp = self.decompose_orbits()
        # self.losvd_histograms, self.proj_mass, self.decomp = self.run_dec()
        self.logger.info('Orbits read and velocity histogram created.')

    def plot_decomp(self, xlim, ylim, v_sigma_option='fit'):
        """ Generate decomposition plots.

        Parameters
        ----------
        xlim : float
            restricts plot x-coordinates to abs(x) <= xlim.
        ylim : float
            restricts plot y-coordinates to abs(y) <= ylim.
        v_sigma_option : str, optional
            If 'fit', v_mean and v_sigma are calculated based on fitting
            Gaussians, if 'moments', v_mean and v_sigma are calculated
            directly from the model's losvd histograms. The default is 'fit'.

        Returns
        -------
        None.

        """
        comp_kinem_moments = self.comps_aphist(v_sigma_option)
        self.logger.info('Component data done.')
        self.plot_comps_giu(xlim=xlim,
                            ylim=ylim,
                            comp_kinem_moments=comp_kinem_moments)
        self.logger.info('Plots done.')

    def comps_aphist(self, v_sigma_option='fit'):
        """Calculate components' flux, mean velocity, and velocity dispersion.


        Parameters
        ----------
        v_sigma_option : str, optional
            If 'fit', v_mean and v_sigma are calculated based on fitting
            Gaussians, if 'moments', v_mean and v_sigma are calculated
            directly from the model's losvd histograms. The default is 'fit'.

        Raises
        ------
        ValueError
            if v_sigma_option is neither 'moments' nor 'fit'.

        Returns
        -------
        comp_flux_v_sigma : astropy table
            The table columns are: aperture index (starting with 0), followed
            by three columns per component holding the flux, mean velocity,
            and velocity dispersion.
            The chosen v_sigma_option is in the table meta data.

        """
        v_sigma_options = ['moments', 'fit']
        if v_sigma_option not in v_sigma_options:
            text = f'Unknown v_sigma_option {v_sigma_option}, ' \
                   f'must be one of {v_sigma_options}.'
            self.logger.error(text)
            raise ValueError(text)
        self.logger.info('Calculating flux, v, and sigma for components '
                         f'{self.decomp.meta["comps"]}, {v_sigma_option=}.')
        comp_flux_v_sigma = astropy.table.Table(
                            {'ap_id':range(self.losvd_histograms.y.shape[-1])},
                            dtype=[int],
                            meta={'v_sigma_option':v_sigma_option})
        for comp in self.decomp.meta['comps']:
            # calculate flux and losvd histograms for component
            orb_sel = np.array([comp in s for s in self.decomp['component']],
                               dtype=bool)
            flux=np.dot(self.proj_mass[orb_sel].T, self.model.weights[orb_sel])
            losvd = np.dot(self.losvd_histograms.y[orb_sel,:,:].T,
                           self.model.weights[orb_sel]).T
            losvd = losvd[np.newaxis]
            self.logger.debug(f'{comp}: {np.count_nonzero(orb_sel)} orbits, '
                              f'{flux.shape=}, {losvd.shape=}.')
            losvd_hist = dyn.kinematics.Histogram(self.losvd_histograms.xedg,
                                                  y=losvd,
                                                  normalise=False)
            if v_sigma_option == 'moments':
                v_mean = np.squeeze(losvd_hist.get_mean())
                v_sigma = np.squeeze(losvd_hist.get_sigma())
            elif v_sigma_option == 'fit':
                v_mean, v_sigma = losvd_hist.get_mean_sigma_gaussfit()
                v_mean = np.squeeze(v_mean)
                v_sigma = np.squeeze(v_sigma)
            else:
                pass
            comp_flux_v_sigma.add_columns([flux, v_mean, v_sigma],
                                           names=[f'{comp}_lsb',
                                                  f'{comp}_v',
                                                  f'{comp}_sig'])
        return comp_flux_v_sigma

    def decompose_orbits(self, ocut=None):
        """Decompose orbits based on lambda_z.


        Parameters
        ----------
        ocut : list of floats, optional
            The cuts in lambda_z. The default is None, which translates to
            ocut=[0.8, 0.25, -0.25], the selection in lambda_z
            following Santucci+22.

        Returns
        -------
        decomp : astropy table
            The table has two columns, ``id`` and ``component``. The former is
            the orbit id (starting with 0), ``component`` is a string
            describing the component(s) an orbit belongs to. Note that an
            orbit can belong to multiple components. In that case, the
            component strings are concatenated. For easier parsing later, the
            component descriptors are surrounded by pipe symbols ``|``.
            The table's meta data ``comps`` holds a list of all components.

        """
        if ocut is None:
            ocut = [0.8, 0.25, -0.25]
        self.logger.debug(f'Cut lines are: {ocut}.')
        file2 = self.model.directory_noml + 'datfil/orblib.dat_orbclass.out'  #orbitlibraries
        file3 = self.model.directory_noml + 'datfil/orblibbox.dat_orbclass.out'
        file3_test = os.path.isfile(file3)
        if not file3_test:
            file3 = '%s' % file2

        n_orb = self.config.settings.orblib_settings['nE'] * \
                self.config.settings.orblib_settings['nI2'] * \
                self.config.settings.orblib_settings['nI3']
        n_dither = self.config.settings.orblib_settings['dithering']
        conversion_factor=self.config.all_models.system.distMPc*1.0e6*1.49598e8

#unused        nrow = norb
        ncol = n_dither ** 3
        orbclass1 = np.genfromtxt(file2).T
        orbclass1 = orbclass1.reshape((5,ncol,n_orb), order='F')
        orbclass2 = np.genfromtxt(file3).T
        orbclass2 = orbclass1.reshape((5,ncol,n_orb), order='F')

        orbw = self.model.weights
        n_orbs = len(orbw)
        self.logger.debug(f'{n_orb=}, {n_orbs=}.')

        orbclass = np.dstack((orbclass1, orbclass1, orbclass2))
        self.logger.debug(f'{len(orbclass) = }.')
        orbclass1a = np.copy(orbclass1)
        orbclass1a[0:3, :, :] *= -1  # the reverse rotating orbits of orbclass

        for i in range(n_orb):
            orbclass[:, :, i * 2] = orbclass1[:, :, i]
            orbclass[:, :, i * 2 + 1] = orbclass1a[:, :, i]

        ## define circularity of each orbit [nditcher^3, n_orb]
        lz = (orbclass[2, :, :] / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))  # lambda_z = lz/(r * Vrms)
#unused        lx = (orbclass[0, :, :] / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))  # lambda_x = lx/(r * Vrms)
#unused        l = (np.sqrt(np.sum(orbclass[0:3, :, :] ** 2, axis=0)) / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))
#unused        r = (orbclass[3, :, :] / conversion_factor)  # from km to kpc

        # average values for the orbits in the same bundle (n_dither^3).
        # Only include the orbits within Rmax_arcs

        rm = np.sum(orbclass[3, :, :]/conversion_factor, axis=0) / n_dither**3
#unused        lzm = np.sum(np.abs(lz), axis=0) / n_dither ** 3
#unused        lxm=np.sum(lx, axis=0) / n_dither ** 3
        #print("check 1", lzm, lxm)
        #s = np.ravel(np.where((rm > xrange[0]) & (rm < xrange[1])))

        # flip the sign of lz to confirm total(lz) > 0
        t = np.ravel(np.argsort(rm))

        yy = np.max(np.ravel(np.where(np.cumsum(orbw[t]) <= 0.5)))
        k = t[0:yy]
        if np.sum(np.sum(lz[:, k], axis=0) / (n_dither ** 3) * orbw[k]) < 0:
            lz *= -1.0

        lzm_sign= np.sum(lz, axis=0) / n_dither ** 3
#unused        lxm_sign= np.sum(lx, axis=0) / n_dither ** 3
        #print("check 2 - sign", lzm_sign, lxm_sign)

        comps=['disk', 'thin_d', 'warm_d', 'bulge', 'all']
        self.logger.info(f'Decomposing {n_orbs} orbits into {comps=}...')
        decomp = astropy.table.Table({'id':range(n_orbs),
                                      'component':['']*n_orbs},
                                     dtype=[int, 'U256'],
                                     meta={'comps':comps})
        # map components
        comp_map = np.zeros(n_orbs, dtype=int)
        # cold component
        comp_map[np.ravel(np.where(lzm_sign >= ocut[0]))] += \
            2**comps.index('thin_d')
        # warm component
        comp_map[np.ravel(np.where((lzm_sign > ocut[1])
                                 & (lzm_sign < ocut[0])))] += \
            2**comps.index('warm_d')
        # hot component
        comp_map[np.ravel(np.where((lzm_sign > ocut[2])
                                 & (lzm_sign <= ocut[1])))] += \
            2**comps.index('bulge') # was lzm_sign<ocut[1]
        # disk component
        comp_map[np.ravel(np.where(lzm_sign > ocut[1]))] += \
            2**comps.index('disk')
        # whole component
        comp_map += 2**comps.index('all')
        for i in np.ravel(np.where(comp_map > 0)):
            for k, comp in enumerate(comps):
                if comp_map[i] & (1 << k):
                    decomp['component'][i] += f'|{comp}|'
        return decomp

    def plot_comps_giu(self,
                       # savedata=True,
                       xlim=None,
                       ylim=None,
                       # Re=None,
                       v_sigma_option=None,
                       comp_kinem_moments=None,
                       figtype='.png'):

        v_sigma_option = comp_kinem_moments.meta['v_sigma_option'] \
                         if 'v_sigma_option' in comp_kinem_moments.meta.keys()\
                         else ''
        self.logger.info(f'Plotting decomposition for {v_sigma_option=}.')

        # read kinematic data and weights
        weights = self.model.weights
        ## COLD COMPONENT
        flux_thin = comp_kinem_moments['thin_d_lsb']
        vel_thin = comp_kinem_moments['thin_d_v']
        sig_thin = comp_kinem_moments['thin_d_sig']
        wthin = weights[['thin_d' in s for s in self.decomp['component']]]

        ## WARM COMPONENT
        flux_thick = comp_kinem_moments['warm_d_lsb']
        vel_thick = comp_kinem_moments['warm_d_v']
        sig_thick = comp_kinem_moments['warm_d_sig']
        wthick = weights[['warm_d' in s for s in self.decomp['component']]]

        ## CC COMPONENT
        flux_disk = comp_kinem_moments['disk_lsb']
        vel_disk = comp_kinem_moments['disk_v']
        sig_disk = comp_kinem_moments['disk_sig']
        wdisk = weights[['disk' in s for s in self.decomp['component']]]

        ## HOT_cr COMPONENT
        flux_bulge = comp_kinem_moments['bulge_lsb']
        vel_bulge = comp_kinem_moments['bulge_v']
        sig_bulge = comp_kinem_moments['bulge_sig']
        wbulge = weights[['bulge' in s for s in self.decomp['component']]]

        ###WHOLE component
        flux_all = comp_kinem_moments['all_lsb']
        vel_all = comp_kinem_moments['all_v']
        sig_all = comp_kinem_moments['all_sig']
        wall = weights[['all' in s for s in self.decomp['component']]]

        # read the pixel grid
        stars = \
        self.config.system.get_component_from_class(
                                dyn.physical_system.TriaxialVisibleComponent)
        dp_args = stars.kinematic_data[self.kin_set].dp_args
        xi = dp_args['x']
        yi = dp_args['y']
        dx = dp_args['dx']
        grid = dp_args['idx_bin_to_pix']
        # The angle that is saved in this file is measured counter clock-wise
        # from the galaxy major axis to the X-axis of the input data.
        angle_deg = dp_args['angle']
        self.logger.debug(f'Pixel grid dimension is {dx=}, {len(xi)=}, '
                          f'{len(yi)=}, {grid.shape}, {angle_deg=}.')

        s = np.ravel(np.where((grid >= 0) & (np.abs(xi) <= xlim)
                              & (np.abs(yi) <= ylim)))
        s_wide = np.ravel(np.where(grid >= 0))
        #normalise flux
        fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_thin))
        flux_thin = flux_thin / fhist
        fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_thick))
        flux_thick = flux_thick / fhist

        fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_disk))
        flux_disk = flux_disk / fhist

        fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_bulge))
        flux_bulge = flux_bulge / fhist

        fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_all))
        flux_all= flux_all / fhist

        tthin  = flux_thin[grid]
        tthick = flux_thick[grid]

        tdisk  = flux_disk[grid]
        tbulge =flux_bulge[grid]
        tall =flux_all[grid]

        #print('before normalization th  tw  tz tc thc tcw are:', np.sum(th),
        #      np.sum(tw), np.sum(tz), np.sum(tc), np.sum(thc),np.sum(tcw))
        #print('before normalization tthin, tthick,tdisk,tbulge, all:',
        #      np.sum(tthin),  np.sum(tthick),  np.sum(tdisk), np.sum(tbulge),
        #      np.sum(tall))

        tthin =tthin *np.sum(wthin)/np.sum(tthin)
        tthick=tthick*np.sum(wthick)/np.sum(tthick)

        tdisk =tdisk *np.sum(wdisk)/np.sum(tdisk)
        tbulge=tbulge*np.sum(wbulge)/np.sum(tbulge)
        tall =tall *np.sum(wall)/np.sum(tall)

        #print('after normalization th  tw  tz tc thc tcware:', np.sum(th),
        #      np.sum(tw), np.sum(tz), np.sum(tc), np.sum(thc), np.sum(tcw))
        #print('after normalization tthin, tthick,tdisk,tbulge:',
        #      np.sum(tthin),  np.sum(tthick),  np.sum(tdisk), np.sum(tbulge),
        #      np.sum(tall))

        #totalf = np.sum(th) + np.sum(tw) + np.sum(tz)+ np.sum(tc)

        ###if you want to check that the sum of the flux of the components is
        #the same as the total, uncomment the two extra totalf and print them
        totalf = np.sum(tthin) + np.sum(tthick) + np.sum(tbulge)
        #totalf2 = np.sum(tdisk)+ np.sum(tbulge)
        #totalf3=np.sum (tall)
        #print('total fluxes with thin thick bulge', totalf,
        #      ' and with disk bulge',totalf2, ' and tot', totalf3 )
        tthin =tthin /totalf
        tthick=tthick/totalf

        tdisk =tdisk /totalf
        tbulge=tbulge/totalf
        tall=tall/totalf
        #tcw=tcw/totalf
        #flux = th +  tw + tz + tc
        flux = tthin +  tthick + tbulge
        #flux2= tdisk+tbulge
        #print("luminosity fractions f_thin, f_thick, f_disk, f_bulge, tot1, "
        #      "tot2 (disk+bulge), all")
        #print(np.sum(th), np.sum(tw), np.sum(tz), np.sum(tc), np.sum(thc),
        #      np.sum(tcw))
        #print(np.sum(tthin),  np.sum(tthick),np.sum(tdisk), np.sum(tbulge),
        #      np.sum(flux), np.sum(flux2), np.sum(tall))

        ### SAVE DATA TO A FILE
       # if savedata:
       #     with open(figdir + 'SB_hmz.dat', 'w') as outfile:
       #          outfile.write('x/arcs,  y/arcs,  SB thin disk,  SB warm disk,'
       #                        '  SB bulge, SB counter_rot' + '\n')
       #         for j in range(0, len(s)):
       #             outfile.write(("%10.4f" % xi[s[j]]) + ("%10.4f" % yi[s[j]])
       #                           +  ("%12.3e" % th[s[j]]) +
       #                           ("%12.3e" % tw[s[j]]) + ("%12.3e" % tz[s[j]])
       #                           + ("%12.3e" % tc[s[j]]) + '\n')

        ### EVALUATE VMAX and SMAX and SMIN

        vmax = np.nanmax([vel_thin,vel_thick, vel_disk, vel_bulge, vel_all])
        sig_t = np.array((sig_thin,sig_thick,sig_disk,sig_bulge, sig_all))
        #vmax = 79


        smax = np.nanmax(sig_t[sig_t > 0])
        smin = np.nanmin(sig_t[sig_t > 0])
        #smax =346
        #smin =161

        minf=min(-2.5 * np.log10(flux))
        maxf=max(-2.5 * np.log10(flux[flux !=0]))
        xi_t= (xi[s])
        yi_t=(yi[s])

        comps_kin = astropy.table.Table({'x/arcs':xi_t,
                                         'y/arcs':yi_t,
                                         'SB_thin_disk':tthin[s],
                                         'SB_thick_disk':tthick[s],
                                         'SB_disk':tdisk[s],
                                         'SB_bulge':tbulge[s],
                                         'SB_whole':tall[s],
                                         'vel_thin_disk':vel_thin[grid[s]],
                                         'vel_thick_disk':vel_thick[grid[s]],
                                         'vel_disk':vel_disk[grid[s]],
                                         'vel_bulge':vel_bulge[grid[s]],
                                         'vel_whole':vel_all[grid[s]],
                                         'sig_thin_disk':sig_thin[grid[s]],
                                         'sig_thick_disk':sig_thick[grid[s]],
                                         'sig_disk':sig_disk[grid[s]],
                                         'sig_bulge':sig_bulge[grid[s]],
                                         'sig_whole':sig_all[grid[s]]})

        kin_name = stars.kinematic_data[self.kin_set].name
        file_name = f'comps_kin_test_s22_{v_sigma_option}_{kin_name}'
        table_file_name = self.model.directory + file_name + '.ecsv'
        plot_file_name = self.config.settings.io_settings['plot_directory'] \
                         + file_name \
                         + figtype
        comps_kin.write(f'{table_file_name}',
                        format='ascii.ecsv',
                        overwrite=True)
        self.logger.info('Component grid kinematics written to '
                         f'{table_file_name}.')

        self.logger.debug(f'{v_sigma_option}: {vmax=}, {smax=}, {smin=}.')
        # print(np.max(th),np.min(th),np.max(tw),np.min(tw),np.max(tz),
        #       np.min(tz),np.max(tc),np.min(tc))

        ### PLOT THE RESULTS
        # Plot settings
        plt.figure(figsize=(12, 18))
        #plt.subplots_adjust(hspace=0.7, wspace=0.01, left=0.01, bottom=0.05,
        #                    top=0.99, right=0.99)
        plt.subplots_adjust(hspace=0.4, wspace=0.02, left=0.01, bottom=0.05,
                            top=0.99, right=0.99)

        ### PLOT THE COMPONENTS
        ## COLD
        ax1=plt.subplot(5, 3, 1)
        display_pixels(xi[s], yi[s], -2.5 * np.log10(tthin[s]) , pixelsize=dx,
                        colorbar=True, nticks=7, cmap='YlOrRd_r',
                        label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
        #ellipse = Ellipse((0, 0),width=Re * 2,height=semi_min * 2,
        #                  facecolor='none', edgecolor = 'black')
       # ax1.add_patch(ellipse)
        ax2=plt.subplot(5, 3, 2)
        plt.title("THIN DISK COMPONENT")
        display_pixels(xi[s], yi[s], vel_thin[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='RdYlBu_r',
                       vmin=-1.0 * vmax, vmax=vmax, label='Velocity')

        #ellipse = Ellipse((0, 0),width=Re * 2,height=semi_min * 2,
        #                  facecolor='none', edgecolor = 'black')
        #ax2.add_patch(ellipse)

        ax3=plt.subplot(5, 3, 3)
        display_pixels(xi[s], yi[s], sig_thin[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='YlOrRd',
                       vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')
        #ellipse = Ellipse((0, 0),width=Re * 2,height=semi_min * 2,
        #                  facecolor='none', edgecolor = 'black')
        #ax3.add_patch(ellipse)
        ## WARM
        plt.subplot(5, 3, 4)
        display_pixels(xi[s], yi[s], -2.5 * np.log10(tthick[s]) , pixelsize=dx,
                       colorbar=True, nticks=7, cmap='YlOrRd_r',
                       label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
        plt.subplot(5, 3, 5)
        plt.title("THICK DISK COMPONENT")
        display_pixels(xi[s], yi[s], vel_thick[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='RdYlBu_r',
                       vmin=-1.0 * vmax, vmax=vmax, label='Velocity')
        plt.subplot(5, 3, 6)
        display_pixels(xi[s], yi[s], sig_thick[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='YlOrRd',
                       vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')

          #HOT+CR
        plt.subplot(5, 3, 7)
        display_pixels(xi[s], yi[s], -2.5 * np.log10(tdisk[s]) , pixelsize=dx,
                       colorbar=True, nticks=7, cmap='YlOrRd_r',
                       label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
        plt.subplot(5, 3, 8)
        #plt.title("HOT + C.R. COMPONENT")
        plt.title("DISK COMPONENT")
        display_pixels(xi[s], yi[s], vel_disk[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='RdYlBu_r',
                       vmin=-1.0 * vmax, vmax=vmax, label='Velocity')
        plt.subplot(5, 3, 9)
        display_pixels(xi[s], yi[s], sig_disk[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='YlOrRd',
                       vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')


        plt.subplot(5, 3, 10)
        display_pixels(xi[s], yi[s], -2.5 * np.log10(tbulge[s]) , pixelsize=dx,
                       colorbar=True, nticks=7, cmap='YlOrRd_r',
                       label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
        plt.subplot(5, 3, 11)
        #plt.title("HOT + C.R. COMPONENT")
        plt.title("BULGE COMPONENT")
        display_pixels(xi[s], yi[s], vel_bulge[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='RdYlBu_r',
                       vmin=-1.0 * vmax, vmax=vmax, label='Velocity')
        plt.subplot(5, 3, 12)
        display_pixels(xi[s], yi[s], sig_bulge[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='YlOrRd',
                       vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')
        #########all
        plt.subplot(5, 3, 13)
        display_pixels(xi[s], yi[s], -2.5 * np.log10(tall[s]) , pixelsize=dx,
                       colorbar=True, nticks=7, cmap='YlOrRd_r',
                       label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
        plt.subplot(5, 3, 14)
        #plt.title("HOT + C.R. COMPONENT")
        plt.title("WHOLE COMPONENT")
        display_pixels(xi[s], yi[s], vel_all[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='RdYlBu_r',
                       vmin=-1.0 * vmax, vmax=vmax, label='Velocity')
        plt.subplot(5, 3, 15)
        display_pixels(xi[s], yi[s], sig_all[grid[s]], pixelsize=dx,
                       colorbar=True, nticks=7, cmap='YlOrRd',
                       vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')

        plt.tight_layout()
        plt.savefig(plot_file_name)
        self.logger.info(f'Component plots written to {plot_file_name}.')
        plt.close()
