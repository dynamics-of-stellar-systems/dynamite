import logging
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt
from plotbin import sauron_colormap as pb_sauron_colormap
import dynamite as dyn
from plotbin import display_pixels
# from loess.loess_2d import loess_2d
import physical_system as physys


class Plotter():

    def __init__(self,
                 system=None,
                 settings=None,
                 parspace=None,
                 all_models=None):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.system = system
        self.settings = settings
        self.parspace = parspace
        self.all_models = all_models
        self.input_directory = settings.io_settings['input_directory']
        pb_sauron_colormap.register_sauron_colormap()

        # parspace is a list
        # par0 = parspace[0]
        # par0.fixed is True/False ...

        # self.settings.io_settings has the main output directory
        # for individual models, you can create a model.SchwarzschildModel object
        # and use its get_model_directory method

        # code below copied from schw_schwmodplot.py
        '''
        if not rootname: rootname = ['']
        if not kin_plot: kin_plot=0
        if not orbit_plot: orbit_plot=0
        if not mass_plot: mass_plot=0
        '''
        # <-- dont need to implement "which plot to make" switches yet,
        #     just individual plotting methods

        '''
        figdir = w_dir + object + '/figure_nn'+rootname[0]+'/'
        if not os.path.isdir(figdir): file_mkdir(figdir)

        pars, chilim, freeroot = schwparinfo(w_dir=w_dir, rootname=rootname, FIRSTITER=0, object=object)

        a = read_file_element(w_dir + object + '/infil/kin_data.dat', [1], [1])
        Nobs = long(a[0])

        b = read_file_element(w_dir + object + '/infil/nn'+ str(rootname[0]) +'.in', [6], [1])
        nGH = long(b[0])

        fchiname1 = w_dir + object + '/griddata/' + str(rootname[0]) + '_chi2.cat'
        '''


        '''
        Nf, npar, mpar, mtime, fls = read_chi2_file(fchiname1)

        chi2 = mpar[npar, :]         # the direct chi^2 from modelling fitting h1, h2, h3, h4..
        chikin = mpar[npar + 1, :]   # comparision between final modle and data, v, sigma, h3, h4

        chi2 = chikin                # in the following, we only use chikin

        # 1 sigma confidence level
        chlim = math.sqrt(2 * Nobs * nGH)
        '''
        # chlim can be found in e.g.
        # self.settings.parameter_space_settings.generator_settings['threshold_del_chi2']

        '''
        #### Normalize the chi2 anf chikin, the minimum chi2 per freedom should be ~ 1

        chi2pmin = (min(chi2)/ (nGH*Nobs))
        chi2 = chi2 * 1.0/ chi2pmin

        par = np.zeros((npar, Nf), dtype=float)
        par[0, :] = mpar[npar - 1, :]
        par[1:npar, :] = mpar[0:npar - 1, :]

        fix = np.ravel(np.where(pars.fixed == 1))
        nofix = np.ravel(np.where(pars.fixed == 0))
        nnofix = len(nofix)

        ml = par[0, :]
        qmin = par[1, :]

        s0 = np.ravel(np.where((chi2 > 0) & (chikin > 0)))
        ns0 = len(s0)
        nf = np.copy(ns0)

        par = par[:, s0]
        chi2 = chi2[s0]
        '''

    def make_chi2_plot(self,parset):

        # parset should be provided that it is clear what values the
        # remaining parameters are marginalized over

        val1=self.all_models.table['m-bh']
        val2=self.all_models.table['ml']
        val3=self.all_models.table['chi2'] - np.min(self.all_models.table['chi2'])

        mass = self.parspace.get_parameter_from_name('m-bh')
        ml = self.parspace.get_parameter_from_name('ml')       

        #print(ml.LaTeX)

        mbh=np.unique(val1)
        ml=np.unique(val2)

        mbh_array0,ml_array=np.meshgrid(mbh,ml)
        chi2_array=mbh_array0*0
        mbh_array=mbh_array0*0

        #not needed at the moment.. maybe later
        #model_dir=self.settings.io_settings['output_directory']+'models/'

        for j in range(len(mbh)):
        #for j in range(2):
            
            for k in range(len(ml)):
            #for k in range(2):

                #change the parameter setting to the parameters of the grid
                parset['ml'] = ml_array[k,j]
                parset['m-bh'] = mbh_array0[k,j]

                # TBD: I need to set up a model in order to re-create the
                # model directory. Can that be changed? To be discussed in 
                # dynamite meeting. Idea: add model directory to all_models
                # table?
                # model = dyn.model.LegacySchwarzschildModel(
                #             system=self.system,
                #             settings=self.settings,
                #             parspace=self.parspace,
                #             parset=parset)
                model = self.all_models.get_model_from_parset(parset)
                # param_fname=str(model.get_model_directory())+'nn.in'
                param_fname = model.get_model_directory() + 'nn.in'

                # extract the ml scale factor from nnls input file. Needs
                # to be changed for non-legacy version. 
                scale_factor=np.genfromtxt(param_fname, skip_header=10,max_rows=1,usecols=0)                

                #each black hole input mass needs to be multiplied with the scale factor.
                mbh_array[k,j]=mbh_array0[k,j]*scale_factor**2

                chi2_array[k,j]=val3[((val1==mbh_array0[k,j]) & (val2==ml_array[k,j]))]

        fig = plt.figure(figsize=(6, 4))
        ax=plt.subplot(1,1,1)

        levels2=np.array([2.3, 6.17, 11.8,11.8*2,11.8*4,11.8*8,11.8*16,11.8*32,11.8*64,11.8**128])
        #levels2=np.array([10**1,10**1.5,10**2,10**2.5,10**3,10**3.5,10**4])
        
        plt.plot(mbh_array,ml_array,marker='o',color='black',linestyle='', markersize=6)
        CS=plt.contour(mbh_array,ml_array,chi2_array,levels=levels2,colors='k')
        plt.contour(mbh_array,ml_array,chi2_array,levels=np.array([11.8*4]),colors='red')
        ax.set_xscale("log")
        plt.gca().set_xlabel(mass.LaTeX)
        plt.gca().set_ylabel('M/L')
        plt.tight_layout()
        plt.show()
        
        return fig


    def plot_kinematic_maps(self, model, kin_set=None):
        """
        Show kinematic map of the best-fitting model, with v, sigma, h3, h4
        Taken from schw_kin.py.
        Note: kin_set should be the index of the data set,
        e.g. kin_set=0 , kin_set=1
        """
        kinem_fname = model.get_model_directory() + 'nn_kinem.out'
        body_kinem = np.genfromtxt(kinem_fname, skip_header=1)

        stars = \
          self.system.get_component_from_class(physys.TriaxialVisibleComponent)
        n_kin = len(stars.kinematic_data)
        if kin_set >= n_kin:
            text = f'kin_set must be < {n_kin}, but it is {kin_set}'
            self.logger.error(text)
            raise ValueError(text)

        first_bin = sum(k.n_apertures for k in stars.kinematic_data[:kin_set])
        n_bins = stars.kinematic_data[kin_set].n_apertures
        body_kinem = body_kinem[first_bin:first_bin+n_bins]
        self.logger.debug(f'kin_set={kin_set}, plotting bins '
                          f'{first_bin} through {first_bin+n_bins-1}')
        # if kin_set==0:
        #     n_bins=stars.kinematic_data[0].n_apertures
        #     body_kinem=body_kinem[0:n_bins,:]
        #     self.logger.info(f'first_bin=0, last_bin={n_bins}')

        # elif kin_set==1:
        #     n_bins1=stars.kinematic_data[0].n_apertures
        #     n_bins2=stars.kinematic_data[1].n_apertures
        #     body_kinem=body_kinem[n_bins1:n_bins1+n_bins2,:]
        #     self.logger.info(f'first_bin={n_bins1}, last_bin={n_bins1+n_bins2}')

        # else:
        #     text = f'kin_set must be 0 or 1, not {kin_set}'
        #     self.logger.error(text)
        #     raise ValueError(text)

        if self.settings.weight_solver_settings['number_GH'] == 2:
            id_num, fluxm, flux, velm, vel, dvel, sigm, sig, dsig = body_kinem.T

            #to not need to change the plotting routine below, higher moments are set to 0
            h3m, h3, dh3, h4m, h4, dh4 = vel*0, vel*0, vel*0+0.4, vel*0, vel*0, vel*0+0.4

        if self.settings.weight_solver_settings['number_GH'] == 4:
            id_num, fluxm, flux, velm, vel, dvel, sigm, sig, dsig, h3m, h3, dh3, h4m, h4, dh4 = body_kinem.T

        if self.settings.weight_solver_settings['number_GH'] == 6:
            id_num, fluxm, flux, velm, vel, dvel, sigm, sig, dsig, h3m, h3, dh3, h4m, h4, dh4, h5m, h5, dh5, h6m, h6, dh6 = body_kinem.T

            #still ToDO: Add the kinematic map plots for h5 and h6


        vmax = max(np.hstack((velm, vel)))
        smax = max(np.hstack((sigm, sig)))
        smin = min(sig)
        h3max = max(np.hstack((h3m, h3)))
        h3min = min(np.hstack((h3m, h3)))
        h4max = max(np.hstack((h4m, h4)))
        h4min = min(np.hstack((h4m, h4)))

        # Read aperture.dat
        # The angle that is saved in this file is measured counter clock-wise
        # from the galaxy major axis to the X-axis of the input data.

        aperture_fname = stars.kinematic_data[kin_set].aperturefile
        aperture_fname = self.input_directory + aperture_fname

        lines = [line.rstrip('\n').split() for line in open(aperture_fname)]
        strhead = lines[0]
        minx = np.float(lines[1][0])
        miny = np.float(lines[1][1])
        sx = np.float(lines[2][0])
        sy = np.float(lines[2][1])
        maxx = sx + minx
        sy = sy + miny
        angle_deg = np.float(lines[3][0])  # - 90.0 + 180
        nx = np.int(lines[4][0])
        ny = np.int(lines[4][1])
        dx = sx / nx

        self.logger.debug(f"Pixel grid dimension is dx={dx},nx={nx},ny={ny}")
        grid = np.zeros((nx, ny), dtype=int)

        xr = np.arange(nx, dtype=float) * dx + minx + 0.5 * dx
        yc = np.arange(ny, dtype=float) * dx + miny + 0.5 * dx

        xi = np.outer(xr, (yc * 0 + 1))
        xt = xi.T.flatten()
        yi = np.outer((xr * 0 + 1), yc)
        yt = yi.T.flatten()

        radeg = 57.2958
        self.logger.debug(f'PA: {angle_deg}')
        xi = xt
        yi = yt

        # read bins.dat

        bin_fname = stars.kinematic_data[kin_set].binfile
        bin_fname = self.input_directory + bin_fname
        lines_bins = [line.rstrip('\n').split() for line in open(bin_fname)]
        i = 0
        str_head = []
        i_var = []
        grid = []
        while i < len(lines_bins):
            for x in lines_bins[i]:
                if i == 0:
                    str_head.append(str(x))
                if i == 1:
                    i_var.append(np.int(x))
                if i > 1:
                    grid.append(np.int(x))
            i += 1
        str_head = str(str_head[0])
        i_var = int(i_var[0])
        grid = np.ravel(np.array(grid))

        # bins start counting at 1 in fortran and at 0 in idl:
        grid = grid - 1

        # Only select the pixels that have a bin associated with them.
        s = np.ravel(np.where((grid >= 0)))
        fhist, fbinedge = np.histogram(grid[s], bins=len(flux))
        flux = flux / fhist
        fluxm = fluxm / fhist

        ### plot settings

        minf = min(np.array(list(map(np.log10, flux[grid[s]] / max(flux)))))
        maxf = max(np.array(list(map(np.log10, flux[grid[s]] / max(flux)))))
        minfm = min(np.array(list(map(np.log10, fluxm[grid[s]] / max(fluxm)))))
        maxfm = max(np.array(list(map(np.log10, fluxm[grid[s]] / max(fluxm)))))

        # The galaxy has NOT already rotated with PA to make major axis aligned with x

        # filename3 = figdir + object + '_kin4.pdf'

        # fig = plt.figure(figsize=(27, 12))
        fig = plt.figure(figsize=(27, 15))
        fig.suptitle(f'Python {self.version_p()}, '
            f'gfortran {self.version_f()}, Thomas\' Mac, '
            f'Random seed: {self.settings.orblib_settings["random_seed"]}',
            size=24)
        plt.subplots_adjust(hspace=0.7,
                            wspace=0.01,
                            left=0.01,
                            bottom=0.05,
                            top=0.99,
                            # top=0.80,
                            right=0.99)
        sauron_colormap = plt.get_cmap('sauron')
        sauron_r_colormap = plt.get_cmap('sauron_r')

        kw_display_pixels = dict(pixelsize=dx,
                                 angle=angle_deg,
                                 colorbar=True,
                                 nticks=7,
                                 cmap='sauron')
        x, y = xi[s], yi[s]

        ### PLOT THE REAL DATA
        plt.subplot(3, 5, 1)
        c = np.array(list(map(np.log10, flux[grid[s]] / max(flux))))
        display_pixels.display_pixels(x, y, c,
                                          vmin=minf, vmax=maxf,
                                          label='log10(SB)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 2)
        display_pixels.display_pixels(x, y, vel[grid[s]],
                                          vmin=-1.0 * vmax, vmax=vmax,
                                          label='Velocity',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 3)
        display_pixels.display_pixels(x, y, sig[grid[s]],
                                          vmin=smin, vmax=smax,
                                          label='Sigma',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 4)
        display_pixels.display_pixels(x, y, h3[grid[s]],
                                          vmin=h3min, vmax=h3max,
                                          label=r'$h_{3}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 5)
        display_pixels.display_pixels(x, y, h4[grid[s]],
                                          vmin=h4min, vmax=h4max,
                                          label=r'$h_{4}$',
                                          **kw_display_pixels)

        ### PLOT THE MODEL DATA
        plt.subplot(3, 5, 6)
        c = np.array(list(map(np.log10, fluxm[grid[s]] / max(fluxm))))
        display_pixels.display_pixels(x, y, c,
                                          vmin=minfm, vmax=maxfm,
                                          label='log10(SB) (model)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 7)
        display_pixels.display_pixels(x, y, velm[grid[s]],
                                          vmin=-1.0 * vmax, vmax=vmax,
                                          label='Velocity (model)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 8)
        display_pixels.display_pixels(x, y, sigm[grid[s]],
                                          vmin=smin, vmax=smax,
                                          label='Sigma (model)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 9)
        display_pixels.display_pixels(x, y, h3m[grid[s]],
                                          vmin=h3min, vmax=h3max,
                                          label=r'$h_{3}$'+' (model)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 10)
        display_pixels.display_pixels(x, y, h4m[grid[s]],
                                          vmin=h4min, vmax=h4max,
                                          label=r'$h_{4}$'+' (model)',
                                          **kw_display_pixels)


        ### PLOT THE ERROR NORMALISED RESIDUALS
        plt.subplot(3, 5, 11)
        c = (fluxm[grid[s]] - flux[grid[s]]) / flux[grid[s]]
        display_pixels.display_pixels(x, y, c,
                                          vmin=-0.05, vmax=0.05,
                                          label=r'$\delta_{flux}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 12)
        c = (velm[grid[s]] - vel[grid[s]]) / dvel[grid[s]]
        display_pixels.display_pixels(x, y, c,
                                          vmin=-10, vmax=10,
                                          label=r'$\delta_{vel}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 13)
        c = (sigm[grid[s]] - sig[grid[s]]) / dsig[grid[s]]
        display_pixels.display_pixels(x, y, c,
                                          vmin=-10, vmax=10,
                                          label=r'$\delta_{sigma}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 14)
        c = (h3m[grid[s]] - h3[grid[s]]) / dh3[grid[s]]
        display_pixels.display_pixels(x, y, c,
                                          vmin=-1, vmax=1,
                                          label=r'$\delta_{h_3}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 15)
        c = (h4m[grid[s]] - h4[grid[s]]) / dh4[grid[s]]
        display_pixels.display_pixels(x, y, c,
                                          vmin=-1, vmax=1,
                                          label=r'$\delta_{h_4}$',
                                          **kw_display_pixels)
        fig.subplots_adjust(left=0.07, wspace=0.3, hspace=0.01)
        kwtext = dict(size=20, ha='center', va='center', rotation=90.)
        fig.text(0.05, 0.9, 'Data', **kwtext)
        fig.text(0.05, 0.53, 'Model', **kwtext)
        fig.text(0.05, 0.2, 'Residual', **kwtext)
        return fig

    def version_p(self):
        return sys.version.split()[0]
    
    def version_f(self):
        v = subprocess.run("gfortran --version", capture_output=True, shell=True, \
            check=True).stdout.decode('utf-8').split(sep='\n')[0].split()[-1]
        return v















    # etc...
