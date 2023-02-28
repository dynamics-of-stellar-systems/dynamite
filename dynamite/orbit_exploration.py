import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from plotbin.display_pixels import display_pixels
import lmfit
import astropy
import dynamite as dyn
from dynamite import pyfort_GaussHerm

class Decomposition:
    conversions = ['gh_expand_around_losvd_mean_and_std_deviation',
                        # 'gh_fit_with_free_v_sigma_params',
                        'gh_fit_with_free_v_sigma_params_fortran',
                        'moments']

    def __init__(self,
                 config=None,
                 model=None,
                 kin_set=0,
                 comps=['disk', 'thin_d', 'warm_d', 'bulge', 'all']):
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
        self.comps = comps
        self.logger.info(f'Performing decomposition into {comps} '
                         f'for kin_set no {kin_set}: '
                         f'{stars.kinematic_data[kin_set].name}')
        self.plotdir = config.settings.io_settings['plot_directory']
        self.losvd_histograms, self.proj_mass, self.decomp = self.run_dec()
        self.logger.info('Orbits read and velocity histogram created.')

    def plot_decomp(self,xlim, ylim, conversion):
        comp_kinem_moments = self.comps_aphist(conversion)
        self.logger.info('Components done')
        self.plot_comps_giu(xlim=xlim,
                            ylim=ylim,
                            comp_kinem_moments=comp_kinem_moments)
        self.logger.info('Plots done')

    def gaussfunc_gh(self, paramsin,x):
        amp=paramsin['amp'].value
        center=paramsin['center'].value
        sig=paramsin['sig'].value
        c1=-np.sqrt(3)
        c2=-np.sqrt(6)
        c3=2/np.sqrt(3)
        c4=np.sqrt(6)/3
        c5=np.sqrt(6)/4
        skew=paramsin['skew'].value
        kurt=paramsin['kurt'].value
        g=(x-center)/sig
        gaustot_gh = amp*np.exp(-.5*g**2)*(1+skew*(c1*g+c3*g**3) +
                                           kurt*(c5+c2*g**2+c4*(g**4)))
        return gaustot_gh

    def gh_fit_with_free_v_sigma_params(self, vhist, b):
        p_gh = lmfit.Parameters()
        p_gh.add('amp', value=np.max(vhist), vary=True)
        p_gh.add('center',
                 value=np.percentile(b, 50),
                 min=np.min(b),
                 max=np.max(b))
        #p_gh.add('sig',value=(b[w]-np.min(b))/2,min=b[w]-np.min(b),max=np.max(b)-b[w]);
        p_gh.add('sig',
                 value=np.abs(np.percentile(b,68)-np.percentile(b, 50)),
                 min=0,
                 max=700)
        p_gh.add('skew', value=0, vary=True, min=None, max=None)
        p_gh.add('kurt', value=0, vary=True, min=None, max=None)
        gausserr_gh = lambda p,x,y: self.gaussfunc_gh(p,x)-y
        fitout_gh = lmfit.minimize(gausserr_gh, p_gh, args=(b,vhist))

#unused        fitted_p_gh = fitout_gh.params
        pars_gh = [fitout_gh.params['amp'].value,
                   fitout_gh.params['center'].value,
                   fitout_gh.params['sig'].value,
                   fitout_gh.params['skew'].value,
                   fitout_gh.params['kurt'].value]
#unused        fit_gh=self.gaussfunc_gh(fitted_p_gh,b)
#unused        resid_gh=fit_gh-vhist
        return pars_gh[1], pars_gh[2]

    def conv_moments(self, vhis, b, dv_obs):
        L_obs_nn = np.where(vhis > 0, vhis, 0.0)

        mu0=np.sum(L_obs_nn, dtype=np.double)*dv_obs
        mu1=np.sum(b*L_obs_nn, dtype=np.double)*dv_obs
        mu2=np.sum((b**2)*L_obs_nn, dtype=np.double)*dv_obs
#unused        mu3=np.sum((b**3)*L_obs_nn, dtype=np.double)*dv_obs
#unused        mu4=np.sum((b**4)*L_obs_nn, dtype=np.double)*dv_obs

        v = mu1 / mu0
        s = np.sqrt(mu0 * mu2 - mu1 ** 2) / mu0
        return v, s

    def gh_fit_with_free_v_sigma_params_fortran(self, vhis,histbinsize,wbin):
        gamm, vmi, sgi, chi2 = pyfort_GaussHerm.gaussfit(veltemp=vhis,
                nvhist=wbin, nvmax=wbin, dvhist=histbinsize, nmc=10)
        return vmi, sgi

    def gh_expand_around_losvd_mean_and_std_deviation(self, orblib, nkin,xedg):
        mod_losvd=dyn.kinematics.Histogram(y=orblib[np.newaxis,:,:], xedg=xedg)
        mod_v = mod_losvd.get_mean()
        mod_sig = mod_losvd.get_sigma()
        return mod_v[0, nkin], mod_sig[0, nkin]

    def comps_aphist(self, conversion):
        if conversion not in self.conversions:
            text = f'Unknown conversion {conversion}, ' \
                   f'must be one of {self.conversions}.'
            self.logger.error(text)
            raise ValueError(text)
        self.logger.info(f'{conversion=}')
        n_orbs, bins, k_bins = self.losvd_histograms.y.shape
        #w = int(head1[0]) #bins
        #n = int(head1[1]) #apertures
        histbinsize = self.losvd_histograms.dx
        self.logger.info(f'shape histbin: {histbinsize.shape}.')

        b = (np.arange(bins , dtype=float) - ((bins+1)/2)) * histbinsize

        comp_kinem_moments = astropy.table.Table({'ap_id':range(k_bins)},
                                                 dtype=[int],
                                                 meta={'conversion':conversion})
        for comp in self.comps:
            # create losvd histograms for component
            orb_sel = np.array([comp in s for s in self.decomp['component']],
                               dtype=bool)
            losvd = np.dot(self.losvd_histograms.y[orb_sel,:,:].T,
                           self.weights[orb_sel]).T

            v_calculate = np.zeros(k_bins, dtype=float)
            s_calculate = np.zeros(k_bins, dtype=float)
            lsb = np.zeros(k_bins, dtype=float)

            #not calculating h3 and h4 at the moment
            for i in range(0, k_bins):
                vhis = losvd[:, i] #losvd for each kinematic bin
                # #diag start
                # if i==0 or True:
                #     with open(wdir+'vhis.out', 'a') as f:
                #         f.write(f'{vhis}')
                #     with open(wdir+'b.out', 'a') as f:
                #         f.write(f'{b}')
                # #diag end
                dv_obs= (np.max(b)-np.min(b))/np.double(len(b))

                if conversion=='gh_fit_with_free_v_sigma_params':
                    #other options are 'losvd_vsig', 'fortran', 'moments'
                    v_i, sigma_i = self.gh_fit_with_free_v_sigma_params(vhis,
                                                                        b)
                    # #diag start
                    # if i==0 or True:
                    #     with open(wdir+'v_i.out', 'a') as f:
                    #         f.write(f'{v_i}')
                    #     with open(wdir+'sigma_i.out', 'a') as f:
                    #         f.write(f'{sigma_i}')
                    # #diag end

                if conversion=='gh_expand_around_losvd_mean_and_std_deviation':
                    v_i, sigma_i = \
                        self.gh_expand_around_losvd_mean_and_std_deviation(
                            losvd,
                            i,
                            xedg=self.losvd_histograms.xedg)

                if conversion=='gh_fit_with_free_v_sigma_params_fortran':
                    v_i, sigma_i = \
                        self.gh_fit_with_free_v_sigma_params_fortran(
                            vhis,
                            histbinsize,
                            (bins-1)/2)

                if conversion=='moments':
                    v_i, sigma_i = self.conv_moments(vhis, b, dv_obs)

                L_obs_nn = np.where(vhis > 0, vhis, 0.0)
                mu0=np.sum(L_obs_nn, dtype=np.double)*dv_obs

                lsb[i] = mu0
                v_calculate[i] = v_i
                s_calculate[i] = sigma_i

            comp_kinem_moments.add_columns([lsb, v_calculate, s_calculate],
                                           names=[f'{comp}_lsb',
                                                  f'{comp}_v',
                                                  f'{comp}_sig'])
        return comp_kinem_moments

    # def create_orbital_component_files_giu(self,
    #                                        ocut=None,
    #                                        Rmax_arcs=None,
    #                                        xrange=None):
    def decompose_orbits(self, ocut=None, Rmax_arcs=None, xrange=None):

        if not ocut:
            ocut = [0.8, 0.25, -0.25]  #selection in lambda_z following Santucci+22
        self.logger.info(f'cut lines are: {ocut}.')

        if not xrange:
            xrange = [0.0, Rmax_arcs]

        # file4 = self.model.directory + 'nn_orb.out'
        file2 = self.model.directory_noml + 'datfil/orblib.dat_orbclass.out'  #orbitlibraries
        file3 = self.model.directory_noml + 'datfil/orblibbox.dat_orbclass.out'
        file3_test = os.path.isfile(file3)
        if not file3_test:
            file3= '%s' % file2

        nre = self.config.settings.orblib_settings['nE']
        nrth = self.config.settings.orblib_settings['nI2']
        nrrad = self.config.settings.orblib_settings['nI3']
        ndither = self.config.settings.orblib_settings['dithering']
        distance = self.config.all_models.system.distMPc
        conversion_factor = distance*1.0e6*1.49598e8

        norb = int(nre * nrth * nrrad)

#unused        nrow = norb
        ncol = int(ndither ** 3)
        #print('norb', norb)
        orbclass1 = np.genfromtxt(file2).T
        orbclass1 = orbclass1.reshape((5,ncol,norb), order='F')
        orbclass2 = np.genfromtxt(file3).T
        orbclass2 = orbclass1.reshape((5,ncol,norb), order='F')

        #print('norb, ndither', np.max(norb), ndither)
        # norbout, ener, i2, i3, regul, orbtype, orbw, lcut = \
        #     np.genfromtxt(file4,
        #                   skip_header=1,
        #                   usecols=(0,1,2,3,4,5,6,7),
        #                   unpack=True)
        orbw = self.weights
        n_orbs = len(orbw)

        #print('ener, i2, i3', np.max(ener), np.max(i2), np.max(i3))
        #print('Maxmin and Minimum(ener)', np.max(ener), np.min(ener))

        orbclass = np.dstack((orbclass1, orbclass1, orbclass2))
        self.logger.info(f'{len(orbclass) = }.')
        orbclass1a = np.copy(orbclass1)
        orbclass1a[0:3, :, :] *= -1  # the reverse rotating orbits of orbclass

        for i in range(int(0), norb):
            orbclass[:, :, i * 2] = orbclass1[:, :, i]
            orbclass[:, :, i * 2 + 1] = orbclass1a[:, :, i]

        ## define circularity of each orbit [nditcher^3, norb]
        lz = (orbclass[2, :, :] / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))  # lambda_z = lz/(r * Vrms)
#unused        lx = (orbclass[0, :, :] / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))  # lambda_x = lx/(r * Vrms)
#unused        l = (np.sqrt(np.sum(orbclass[0:3, :, :] ** 2, axis=0)) / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))
#unused        r = (orbclass[3, :, :] / conversion_factor)  # from km to kpc

        # average values for the orbits in the same bundle (ndither^3).
        # Only include the orbits within Rmax_arcs

        rm = np.sum(orbclass[3, :, :] / conversion_factor, axis=0) / ndither ** 3
#unused        lzm = np.sum(np.abs(lz), axis=0) / ndither ** 3
#unused        lxm=np.sum(lx, axis=0) / ndither ** 3
        #print("check 1", lzm, lxm)
        #s = np.ravel(np.where((rm > xrange[0]) & (rm < xrange[1])))

        # flip the sign of lz to confirm total(lz) > 0
        t = np.ravel(np.argsort(rm))

        yy = np.max(np.ravel(np.where(np.cumsum(orbw[t]) <= 0.5)))
        k = t[0:yy]
        if np.sum(np.sum(lz[:, k], axis=0) / (ndither ** 3) * orbw[k]) < 0:
            lz *= -1.0

        lzm_sign= np.sum(lz, axis=0) / ndither ** 3
#unused        lxm_sign= np.sum(lx, axis=0) / ndither ** 3
        #print("check 2 - sign", lzm_sign, lxm_sign)

        self.logger.info(f'Decomposing {n_orbs} orbits...')
        decomp = astropy.table.Table({'id':range(n_orbs),
                                      'component':['']*n_orbs},
                                     dtype=[int, 'U256'])

        # map components
        comp_map = np.zeros(n_orbs, dtype=int)
        # cold component
        comp_map[np.ravel(np.where(lzm_sign >= ocut[0]))] += \
            2**self.comps.index('thin_d')
        # warm component
        comp_map[np.ravel(np.where((lzm_sign > ocut[1])
                                 & (lzm_sign < ocut[0])))] += \
            2**self.comps.index('warm_d')
        # hot component
        comp_map[np.ravel(np.where((lzm_sign > ocut[2])
                                 & (lzm_sign < ocut[1])))] += \
            2**self.comps.index('bulge')
        # disk component
        comp_map[np.ravel(np.where(lzm_sign > ocut[1]))] += \
            2**self.comps.index('disk')
        # whole component
        comp_map += 2**self.comps.index('all')
        for i in np.ravel(np.where(comp_map > 0)):
            for k, comp in enumerate(self.comps):
                if comp_map[i] & (1 << k):
                    decomp['component'][i] += f'|{comp}|'

        return decomp

    def plot_comps_giu(self,
                       # savedata=True,
                       xlim=None,
                       ylim=None,
                       # Re=None,
                       conversion=None,
                       comp_kinem_moments=None,
                       figtype='.png'):

        conversion = comp_kinem_moments.meta['conversion'] \
                     if 'conversion' in comp_kinem_moments.meta.keys() else ''
        self.logger.info(f'Plotting decomposition for {conversion=}.')

        # read kinematic data and weights
        ## COLD COMPONENT
        flux_thin = comp_kinem_moments['thin_d_lsb']
        vel_thin = comp_kinem_moments['thin_d_v']
        sig_thin = comp_kinem_moments['thin_d_sig']
        wthin = self.weights[['thin_d' in s for s in self.decomp['component']]]

        ## WARM COMPONENT
        flux_thick = comp_kinem_moments['warm_d_lsb']
        vel_thick = comp_kinem_moments['warm_d_v']
        sig_thick = comp_kinem_moments['warm_d_sig']
        wthick=self.weights[['warm_d' in s for s in self.decomp['component']]]

        ## CC COMPONENT
        flux_disk = comp_kinem_moments['disk_lsb']
        vel_disk = comp_kinem_moments['disk_v']
        sig_disk = comp_kinem_moments['disk_sig']
        wdisk = self.weights[['disk' in s for s in self.decomp['component']]]

        ## HOT_cr COMPONENT
        flux_bulge = comp_kinem_moments['bulge_lsb']
        vel_bulge = comp_kinem_moments['bulge_v']
        sig_bulge = comp_kinem_moments['bulge_sig']
        wbulge = self.weights[['bulge' in s for s in self.decomp['component']]]

        ###WHOLE component
        flux_all = comp_kinem_moments['all_lsb']
        vel_all = comp_kinem_moments['all_v']
        sig_all = comp_kinem_moments['all_sig']
        wall = self.weights[['all' in s for s in self.decomp['component']]]

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
        file_name = f'comps_kin_test_s22_{conversion}_{kin_name}'
        table_file_name = self.model.directory + file_name + '.ecsv'
        plot_file_name = self.plotdir + file_name + figtype
        comps_kin.write(f'{table_file_name}',
                        format='ascii.ecsv',
                        overwrite=True)
        self.logger.info('Component grid kinematics written to '
                         f'{table_file_name}.')

        self.logger.debug(f'{conversion}: {vmax=}, {smax=}, {smin=}.')
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

    def run_dec(self):

        Rmax_arcs=15#gal_infos['Rmax[arcsec]'][i]
        xrange = None

        orblib = self.model.get_orblib()
        orblib.read_losvd_histograms()
        losvd_histograms, _, n_orbs, projected_masses = \
            orblib.losvd_histograms[self.kin_set], orblib.intrinsic_masses, \
            orblib.n_orbs, orblib.projected_masses[self.kin_set]
        self.logger.debug(f'{type(n_orbs)=}, {losvd_histograms.y.shape=}, '
                          f'{projected_masses.shape=}.')
        _ = self.model.get_weights(orblib)
        self.weights = self.model.weights
        #create the files with the orbits selected for each components
        decomp = self.decompose_orbits(Rmax_arcs=Rmax_arcs, xrange=xrange)
        return losvd_histograms, projected_masses, decomp
