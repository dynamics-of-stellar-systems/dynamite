import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from plotbin.display_pixels import display_pixels
import lmfit
import dynamite as dyn
from dynamite import pyfort_GaussHerm

class Decomposition:
    def __init__(self, config=None, model=None, kin_set=0):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None:
            text = f'{__class__.__name__} needs configuration object, ' \
                   'None provided.'
            self.logger.error(text)
            raise ValueError(text)
        self.config = config
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
        if model is None:
            # Select the best model for decomposition
            best_model_idx = config.all_models.get_best_n_models_idx(n=1)[0]
            self.model = config.all_models.get_model_from_row(best_model_idx)
        self.results_directory = \
            config.settings.io_settings['output_directory'] + 'figure_nn/'
        if not os.path.exists(self.results_directory):
            os.makedirs(self.results_directory)
            self.logger.debug(f'Created directory {self.results_directory}')
        else:
            self.logger.debug('Using existing directory '
                              f'{self.results_directory}')
        self.losvd_histograms, self.proj_mass = self.run_dec()
        self.logger.info('Orbits read and velocity histogram created.')
        self.comps = ['disk', 'thin_d', 'warm_d', 'bulge', 'all']
        self.conversions = ['gh_expand_around_losvd_mean_and_std_deviation',
                            # 'gh_fit_with_free_v_sigma_params',
                            'gh_fit_with_free_v_sigma_params_fortran',
                            'moments']
        # #diag start
        # with open(self.results_directory + 'losvd.out', 'w') as f:
        #     f.write(f'{self.losvd_histograms=}\n'
        #             f'{self.losvd_histograms.y.shape=}\n'
        #             f'{self.losvd_histograms.dx=}\n'
        #             f'{self.losvd_histograms.xedg=}\n'
        #             f'{list(self.losvd_histograms.y.flatten())=}\n'
        #             f'{list(self.proj_mass.flatten())=}')
        # #diag end

    def plot_decomp(self,xlim, ylim, conversion):
        self.comps_aphist(conversion)
        self.logger.info('Components done')
        self.plot_comps_giu(xlim=xlim, ylim=ylim, conversion=conversion)
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
        wdir = self.results_directory
        n_orbs, bins, k_bins = self.losvd_histograms.y.shape
        #w = int(head1[0]) #bins
        #n = int(head1[1]) #apertures
        histbinsize = self.losvd_histograms.dx
        self.logger.info(f'shape histbin: {histbinsize.shape}.')

        b = (np.arange(bins , dtype=float) - ((bins+1)/2)) * histbinsize

        xedg = self.losvd_histograms.xedg
        losvd_orbs = self.losvd_histograms.y
        self.logger.info(f'{losvd_orbs.shape = }.')
        for comp in self.comps:
            file=wdir+comp+'_orb_s22.out'
            norbout, orbw = np.genfromtxt(file,
                                          skip_header=1,
                                          usecols=(0,6),
                                          unpack=True)
            # norm_w=np.sum(orbw)
            #print('tot orb_weight', norm_w)
            orb_sel=norbout
            #print('selected orbs',len(orb_sel))
            losvd=np.copy(losvd_orbs)
            #print('before zeroing some orbits', np.sum(np.sum(losvd,0)))
            check=0
            #if k==3: print('orbs', orb_sel)
            for orb in range(n_orbs):

                if orb in orb_sel:
                    #print('orb #', orb)
                    lk=np.where(norbout==orb)[0][0]
                    for i in range(losvd.shape[1]):
                        losvd[orb,i,:]=losvd[orb,i,:]*orbw[lk]
                    #if check==1:print('check orbit taken',losvd[orb,:,0],losvd[orb,:,1])
                else:
                    for i in range(losvd.shape[1]):
                        losvd[orb,i,:]=0.0
                check+=1

            losvd=np.sum(losvd, 0) #sum the losvd in all the orbits
            v_calculate = np.zeros(k_bins, dtype=float)
            s_calculate = np.zeros(k_bins, dtype=float)
            lsb = np.zeros(k_bins, dtype=float)

            #not calculating h3 and h4 at the moment
            self.logger.info(f'{conversion=}')
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
                            xedg)

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

            file_o=wdir+comp+'_vs_kinem_'+conversion+'.out'

            with open(file_o, 'w') as outfile:
                outfile.write(("%13.3f" % k_bins) + ("%8i" % 2 ) + '\n')
                for i in range(0, k_bins):
                    if i==0:
                        self.logger.info(f'{comp = }, '
                                         f'{lsb[i] = }, '
                                         f'{v_calculate[i] = }, '
                                         f'{s_calculate[i] = }.')
                    outfile.write(("%12i" % (i+1))  + ("%15.7f" % lsb[i]) +  ("%15.7f" % lsb[i]) +
                                  ("%15.7f" % v_calculate[i]) + ("%15.7f" % v_calculate[i])  +
                                  ("%15.7f" % 5.0)  + ("%15.7f" % s_calculate[i]) + ("%15.7f" % s_calculate[i]) +
                                  ("%15.7f" % 5.0)  + ("%15.7f" % 0.0) + ("%15.7f" % 0.0)+ ("%15.7f" % 0.1) +
                                  ("%15.7f" % 0.0)  + ("%15.7f" % 0.0) + ("%15.7f" % 0.1) + '\n')

    def create_orbital_component_files_giu(self,
                                           ocut=None,
                                           Rmax_arcs=None,
                                           xrange=None):

        if not ocut:
            ocut = [0.8, 0.25, -0.25]  #selection in lambda_z following Santucci+22
        self.logger.info(f'cut lines are: {ocut}.')

        if not xrange:
            xrange = [0.0, Rmax_arcs]

        file4 = self.model.directory + 'nn_orb.out'
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
        norbout, ener, i2, i3, regul, orbtype, orbw, lcut = \
            np.genfromtxt(file4,
                          skip_header=1,
                          usecols=(0,1,2,3,4,5,6,7),
                          unpack=True)

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

        self.logger.info('**Create nn_orb.out for different bulk of orbits**')

        self.logger.info(f'#orbs: {len(norbout)}.')
        ### COLD COMPONENT
        hlz=np.ravel(np.where((lzm_sign) >= ocut[0]))
        self.logger.info(f'check 3 - cold comp: {len(hlz)}.')
        #print("and len(ocut)", len(ocut))

        if len(ocut) == 1:
            hlz=np.ravel(np.where((lzm_sign) >= ocut[0]))
        orbw_hlz = np.copy(orbw)

        nohlz= np.ravel(np.where((lzm_sign) < ocut[0]))
        if len(ocut) == 1:
            nohlz=np.ravel(np.where((lzm_sign) < ocut[0]))
        orbw_hlz[nohlz] = -999

        with open(self.results_directory+'thin_d_orb_s22.out', 'w') as outfile:
            norbw = len(orbw)
            outfile.write(("%12i" % len(hlz)) + '\n')
            for i in range(0, norbw):
                if orbw_hlz[i]!=-999:
                    outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                              ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                              ("%8i" % orbtype[i])  + ("%13.5e" % orbw_hlz[i]) +
                                  ("%13.5e" % lcut[i]) +  '\n')

        ### WARM COMPONENT
        wlz= np.ravel(np.where((((lzm_sign) > ocut[1])) & ((lzm_sign)< ocut[0])))
        orbw_wlz = np.copy(orbw)
        nowlz= np.ravel(np.where(((lzm_sign) > ocut[0]) | ((lzm_sign) <= ocut[1])))
        orbw_wlz[nowlz] = -999
        self.logger.info(f'check 3 - warm comp: {len(wlz)}.')
        with open(self.results_directory+'warm_d_orb_s22.out', 'w') as outfile:
            norbw = len(orbw)
            outfile.write(("%12i" % len(wlz)) + '\n')
            for i in range(0, norbw):
                if orbw_wlz[i]!=-999:
                    outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                              ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                              ("%8i" % orbtype[i])  + ("%13.5e" % orbw_wlz[i]) +
                                  ("%13.5e" % lcut[i]) +'\n')

        ### HOT_noCR COMPONENT
        zlz= np.ravel(np.where((((lzm_sign) > ocut[2])) & ((lzm_sign)< ocut[1])))

        orbw_zlz = np.copy(orbw)

        nozlz= np.ravel(np.where(((lzm_sign) > ocut[1]) | ((lzm_sign) <= ocut[2])))
        orbw_zlz[nozlz] = -999
        self.logger.info(f'check 3 - hot comp: {len(zlz)}.')
        with open(self.results_directory+'bulge_orb_s22.out', 'w') as outfile:
            norbw = len(orbw)
            outfile.write(("%12i" % len(zlz)) + '\n')
            for i in range(0, norbw):
                if orbw_zlz[i]!=-999:
                    outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                              ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                              ("%8i" % orbtype[i])  + ("%13.5e" % orbw_zlz[i]) +
                                  ("%13.5e" % lcut[i]) +'\n')

     ### DISK COMPONENT
        dlz= np.ravel(np.where(((lzm_sign)> ocut[1])))
        orbw_dlz = np.copy(orbw)

        nodlz= np.ravel(np.where(((lzm_sign) <= ocut[1]) ))
        orbw_dlz[nodlz] = -999
        self.logger.info(f'check 3 - disk comp: {len(dlz)}.')

        with open(self.results_directory + 'disk_orb_s22.out', 'w') as outfile:
            norbw = len(orbw)
            outfile.write(("%12i" % len(dlz)) + '\n')
            for i in range(0, norbw):
                if orbw_dlz[i]!=-999:
                    outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                              ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                              ("%8i" % orbtype[i])  + ("%13.5e" % orbw_dlz[i]) +
                                  ("%13.5e" % lcut[i]) +'\n')

      ### Whole COMPONENT
        orbw_allz = np.copy(orbw)

        self.logger.info(f'check 3 - whole comp: {len(orbw_allz)}.')
        with open(self.results_directory + 'all_orb_s22.out', 'w') as outfile:
            norbw = len(orbw)
            outfile.write(("%12i" % len(orbw)) + '\n')
            for i in range(0, norbw):
                if orbw_allz[i]!=-999:
                    outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                              ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                              ("%8i" % orbtype[i])  + ("%13.5e" % orbw_allz[i]) +
                                  ("%13.5e" % lcut[i]) +'\n')

#unused        tot_orb=np.sum(orbw_hlz[orbw_hlz!=-999])+ np.sum(orbw_wlz[orbw_wlz!=-999])+np.sum(orbw_zlz[orbw_zlz!=-999])
#unused        tot_orb2=np.sum(orbw_zlz[orbw_zlz!=-999])+ np.sum(orbw_dlz[orbw_dlz!=-999])
#unused        tot_orb_tot=np.sum(orbw_allz)
        #print('*************************')
        #print('total flux 1 and 2 (disk+bulge), and total', tot_orb, tot_orb2, tot_orb_tot)
        #print('*************************')

        #print('*************************')
        #print('total weights of hlz (thin), wlz (thick), zlz (bulge), dlz (disk), whole components are:',
        #      np.sum(orbw_hlz[orbw_hlz!=-999])/tot_orb,
        #      np.sum(orbw_wlz[orbw_wlz!=-999])/tot_orb,
        #      np.sum(orbw_zlz[orbw_zlz!=-999])/tot_orb,
        #      np.sum(orbw_dlz[orbw_dlz!=-999])/tot_orb,
        #      np.sum(orbw_allz)/tot_orb)
        #print('************************')
        #print('*************************')
        #print('total weights of hlz (cold), wlz (warm), zlz (hot), clz (CC) components are:',
        #     np.sum(orbw_hlz), np.sum(orbw_wlz),np.sum(orbw_zlz), np.sum(orbw_clz))
        #print('************************')

    def plot_comps_giu(self,
                       # savedata=True,
                       xlim=None,
                       ylim=None,
                       # Re=None,
                       conversion=None):
        if conversion not in self.conversions:
            text = f'Unknown conversion {conversion}, ' \
                   f'must be one of {self.conversions}.'
            self.logger.error(text)
            raise ValueError(text)

        wdir = self.results_directory

        #READING KINEMATIC DATA GENERATED BY triax_makevs_from_aphist
        ##############################################################
        ## COLD COMPONENT
        kinem_file = wdir + 'thin_d_vs_kinem_'+conversion+'.out'
        weight_file = wdir + 'thin_d_orb_s22.out'
        kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
        weight_matrix = np.genfromtxt(weight_file, skip_header=1)
        ident, flux, flux_thin, vel_thin, veld, dvel, sig_thin, sigd, dsig, \
            h3_thin, h3d, dh3, h4_thin, h4d, dh4 = kinem_matrix.T

        n,n1,n2,n3,n0,nc,wthin,lcutw = weight_matrix.T


        ## WARM COMPONENT
        kinem_file = wdir + 'warm_d_vs_kinem_'+conversion+'.out'
        weight_file = wdir + 'warm_d_orb_s22.out'
        kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
        weight_matrix = np.genfromtxt(weight_file, skip_header=1)
        ident, flux, flux_thick, vel_thick, veld, dvel, sig_thick, sigd, dsig, \
            h3_thick, h3d, dh3, h4_thick, h4d, dh4 = kinem_matrix.T
        n,n1,n2,n3,n0,nc,wthick,lcutthi = weight_matrix.T


       ### CC COMPONENT
        kinem_file = wdir + 'disk_vs_kinem_'+conversion+'.out'
        weight_file = wdir + 'disk_orb_s22.out'
        kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
        weight_matrix = np.genfromtxt(weight_file, skip_header=1)
        ident, flux, flux_disk, vel_disk, veld, dvel, sig_disk, sigd, dsig, \
            h3_disk, h3d, dh3, h4_disk, h4d, dh4 = kinem_matrix.T
        n,n1,n2,n3,n0,nc,wdisk,lcutdisk = weight_matrix.T
        #print('disk flux', np.sum(flux_disk))

        ## HOT_cr COMPONENT
        kinem_file = wdir + 'bulge_vs_kinem_'+conversion+'.out'
        weight_file = wdir + 'bulge_orb_s22.out'
        kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
        weight_matrix = np.genfromtxt(weight_file, skip_header=1)
        ident, flux, flux_bulge, vel_bulge, veld, dvel, sig_bulge, sigd, dsig,\
            h3_bulge, h3d, dh3, h4_bulge, h4d, dh4 = kinem_matrix.T
        n,n1,n2,n3,n0,nc,wbulge,lcutbul = weight_matrix.T
        #print('bulge flux', np.sum(flux_bulge))

        ###WHOLE component
        kinem_file = wdir + 'all_vs_kinem_'+conversion+'.out'
        weight_file = wdir + 'all_orb_s22.out'
        kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
        weight_matrix = np.genfromtxt(weight_file, skip_header=1)
        ident, flux, flux_all, vel_all, veld, dvel, sig_all, sigd, dsig, \
            h3_all, h3d, dh3, h4_all, h4d, dh4 = kinem_matrix.T
        n,n1,n2,n3,n0,nc,wall,lcutal = weight_matrix.T

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

#         input_directory = self.config.settings.io_settings['input_directory']

#         ##############################################################
#         # Read aperture1.dat
#         # The angle that is saved in this file is measured counter clock-wise
#         # from the galaxy major axis to the X-axis of the input data.
#         lines = [line.rstrip('\n').split()
#                  for line in open(input_directory + 'aperture.dat')]
# #unused        strhead = lines[0]
#         minx =float(lines[1][0])
#         miny = float(lines[1][1])
#         sx = float(lines[2][0])
#         sy = float(lines[2][1])
#         sy = sy + miny
#         angle_deg = float(lines[3][0])  # - 90.0 + 180
#         nx = int(lines[4][0])
#         ny = int(lines[4][1])

#         dx = sx / nx

#         self.logger.info('Pixel grid dimension is dx,nx,ny,angle: '
#                          f'{dx}, {nx}, {ny}, {angle_deg}.')
#         grid = np.zeros((nx, ny), dtype=int)

#         xr = (np.arange(nx, dtype=float) * dx + minx + 0.5 * dx)
#         yc = (np.arange(ny, dtype=float) * dx + miny + 0.5 * dx)

#         xi = np.outer(xr, (yc * 0 + 1))
#         xt = xi.T.flatten()
#         yi = np.outer((xr * 0 + 1), yc)
#         yt = yi.T.flatten()

#         xi=xt
#         yi=yt
#         ##############################################################
#         # Read bins1.dat
#         lines_bins = [line.rstrip('\n').split()
#                       for line in open(input_directory + 'bins.dat')]
#         i = 0
#         str_head = []
#         i_var = []
#         grid = []
#         while i < len(lines_bins):
#             for x in lines_bins[i]:
#                 if i == 0:
#                     str_head.append(str(x))
#                 if i == 1:
#                     i_var.append(int(x))
#                 if i > 1:
#                     grid.append(int(x))
#             i += 1
#         str_head = str(str_head[0])
#         i_var = int(i_var[0])
#         grid = np.ravel(np.array(grid))
#         # bins start counting at 1 in fortran and at 0 in idl:
#         grid = grid - 1
#         ##############################################################

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

        with open(wdir+'kin_comps_test_s22_'+conversion+'.dat', 'w') as outfile:
            outfile.write('x/arcs  y/arcs  SB_thin_disk  SB_thick_disk  '
                          'SB_disk  SB_bulge SB_whole vel_thin_disk  '
                          'vel_thick_disk  vel_disk  vel_bugle vel_whole '
                          'sig_thin_disk  sig_thick_disk sig_disk  sig_bulge '
                          'sig_whole' + '\n')
            for j in range(0, len(s)):
                outfile.write(("%10.4f" % xi_t[j])  + ("%10.4f" % yi_t[j]) +
                    ("%12.3e" % tthin[s[j]]) + ("%12.3e" % tthick[s[j]]) +
                    ("%12.3e" % tdisk[s[j]])+ ("%12.3e" % tbulge[s[j]]) +
                    ("%12.3e" % tall[s[j]]) +
                    ("%12.3e" % vel_thin[grid[s[j]]]) +
                    ("%12.3e" % vel_thick[grid[s[j]]])+
                    ("%12.3e" % vel_disk[grid[s[j]]]) +
                    ("%12.3e" % vel_bulge[grid[s[j]]])+
                    ("%12.3e" % vel_all[grid[s[j]]]) +
                    ("%12.3e" % sig_thin[grid[s[j]]]) +
                    ("%12.3e" % sig_thick[grid[s[j]]])+
                    ("%12.3e" % sig_disk[grid[s[j]]]) +
                    ("%12.3e" % sig_bulge[grid[s[j]]])+
                    ("%12.3e" % sig_all[grid[s[j]]]) +'\n')

        self.logger.info(f'{conversion}: {vmax=}, {smax=}, {smin=}.')
        # print(np.max(th),np.min(th),np.max(tw),np.min(tw),np.max(tz),
        #       np.min(tz),np.max(tc),np.min(tc))

        ### PLOT THE RESULTS
        # Plot settings
        figfile=f'{wdir}comps_kin_test_s22_{conversion}_kinset{self.kin_set}.pdf'
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
        plt.savefig(figfile)
        plt.close()

    def run_dec(self):

        ocut = [0.8, 0.25,-0.25] #cuts in lambda_z for the components (as in Santucci+22)

        Rmax_arcs=15#gal_infos['Rmax[arcsec]'][i]
        xrange = None

        orblib = self.model.get_orblib()
        orblib.read_losvd_histograms()
        losvd_histograms, _, n_orbs, projected_masses = \
            orblib.losvd_histograms[self.kin_set], orblib.intrinsic_masses, \
            orblib.n_orbs, orblib.projected_masses[self.kin_set]
        self.logger.debug(f'{type(losvd_histograms)=}, {type(n_orbs)=}, '
                          f'{type(projected_masses)=}, '
                          f'{losvd_histograms.y.shape=}, '
                          f'{projected_masses.shape=}.')
        #create the files with the orbits selected for each components
        self.create_orbital_component_files_giu(ocut=ocut,
                                                Rmax_arcs=Rmax_arcs,
                                                xrange=xrange)
        return losvd_histograms, projected_masses
