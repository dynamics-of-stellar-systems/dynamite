import logging
import subprocess
import sys
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
from plotbin import sauron_colormap as pb_sauron_colormap
import dynamite as dyn
import kinematics
import weight_solvers
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

    def make_chi2_vs_model_id_plot(self, which_chi2=None):
        """
        Generates a (kin)chi2 vs. model id plot.

        Parameters
        ----------
        which_chi2 : STR, optional
            Determines whether chi2 or kinchi2 is used. If None, the setting
            in the configuration file's parameter settings is used.
            Must be None, 'chi2', or 'kinchi2'. The default is None.

        Raises
        ------
        ValueError
            If which_chi2 is not one of None, 'chi2', or 'kinchi2'.

        Returns
        -------
        fig : matplotlib.pyplot.figure
            Figure instance.

        """
        if which_chi2==None:
            which_chi2 = self.settings.parameter_space_settings['which_chi2']
        if which_chi2 != 'chi2' and which_chi2 != 'kinchi2':
            text = 'which_chi2 needs to be chi2 or kinchi2, ' \
                   f'but it is {which_chi2}'
            self.logger.error(text)
            raise ValueError(text)
        n_models = len(self.all_models.table)
        fig = plt.figure()
        plt.plot([i for i in range(n_models)],
                  self.all_models.table[which_chi2],
                  'rx')
        plt.gca().set_title(f'{which_chi2} vs. model id')
        plt.xlabel('model id')
        plt.ylabel(which_chi2)
        fig.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        self.logger.info(f'{which_chi2} vs. model id plot created '
                         f'({n_models} models).')
        return fig

    def make_chi2_plot(self, which_chi2=None):
        """
        This implementation is still EXPERIMENTAL. Don't use unless you
        know what you are doing...
        """

        if which_chi2==None:
            which_chi2 = self.settings.parameter_space_settings['which_chi2']
        if which_chi2 != 'chi2' and which_chi2 != 'kinchi2':
            text = 'which_chi2 needs to be chi2 or kinchi2, ' \
                   f'but it is {which_chi2}'
            self.logger.error(text)
            raise ValueError(text)
        self.logger.info(f'Making chi2 plot scaled according to {which_chi2}')

        #Note: it could be a nice feature to exclude the first 50, 100 (specified by the user) or so models in case the values were really off there or
        #alternatively based on too big Delta chi2
        #TODO: after each iteration create chi2 plot
        pars=self.parspace
        val=self.all_models.table

        #only use models that are finished
        val=val[val['all_done']==True]

        #because of the large parameter range dh properties and black hole are plotted in log
        val['c-dh']=np.log10(val['c-dh'])
        val['f-dh']=np.log10(val['f-dh'])
        val['m-bh']=np.log10(val['m-bh'])
        #TBD: add black hole scaling
        #TBD: check the definition of chi2

        #get number and names of parameters that are not fixed
        nofix_sel=[]
        nofix_name=[]
        nofix_latex=[]

        for i in np.arange(len(pars)):
            if pars[i].fixed==False:

                pars[i].name
                nofix_sel.append(i)
                if pars[i].name == 'ml':
                    # nofix_name=np.insert(nofix_name,0,'ml')         #Note: in the old version ml was in the first column, remove hard-coding
                    # nofix_latex=np.insert(nofix_latex,0,'$Y_{r}$')
                    nofix_name.insert(0, 'ml')
                    nofix_latex.insert(0, '$Y_{r}$')
                else:
                    nofix_name.append(pars[i].name)
                    nofix_latex.append(pars[i].LaTeX)

        nnofix=len(nofix_sel)

        nf=len(val)

        nGH=4 #self.weight_solver_settings['number_GH']
        #stars = self.system.get_component_from_class( \
        #                                    physys.TriaxialVisibleComponent)
        #kinematics = stars.kinematic_data.data
        #print(len(kinematics))
        Nobs=353

        self.logger.info(f'nGH={nGH}, Nobs={Nobs}')
        #this is from previous code
        ## 1 sigma confidence level
        #chlim = np.sqrt(2 * Nobs * nGH)
        ##read the chi2 and for some reason normalization, see schw_schwmodplot
        chi2pmin=np.min(val[which_chi2])
        #chi2 = self.all_models.table['chi2']* 1.0/ chi2pmin

        #sabine's code
        chlim = np.sqrt(2 * Nobs * nGH)

        chi2=val[which_chi2]
        chi2-=chi2pmin

        #start of the plotting
        colormap = plt.get_cmap('Spectral')

        fig = plt.figure(figsize=(14, 14))

        #loop over each parameter pair
        for i in range(0, nnofix - 1):
            for j in range(nnofix-1, i, -1):

                xtit = ''
                ytit = ''

                if i==0 : ytit = nofix_latex[j]
                xtit = nofix_latex[i]

                pltnum = (nnofix-1-j) * (nnofix-1) + i+1
                ax = fig.add_subplot(nnofix-1, nnofix-1, pltnum)

                ax.plot(val[nofix_name[i]],val[nofix_name[j]], 'D', color='black', markersize=1.0)
                ax.set_xlabel(xtit, fontsize=12)
                ax.set_ylabel(ytit, fontsize=12)

                #color Delta chi2
                for k in range(nf - 1, -1, -1):
                    #NOTE: I have re-written this part as the old routine did not work properly for me. Would be good to double-check
                    if chi2[k]/chlim<=3: #only significant chi2 values

                        color = colormap(chi2[k]/chlim * 240) #colours the significant chi2
                        markersize = 5-(chi2[k]/chlim) #smaller chi2 become bigger :)
                        ax.plot((val[nofix_name[i]])[k], (val[nofix_name[j]])[k], 'o', markersize=markersize, color=color)

                    if chi2[k]==0:
                        ax.plot((val[nofix_name[i]])[k], (val[nofix_name[j]])[k], 'x', markersize=8, color='black')

        #Note: not so nice is that the red best-fit points are sometimes overplotted with other models. Can probably be improved

        ## add the colorbar with the inverse sauron colormap
        axcb = fig.add_axes([0.55, 0.35, 0.2, 0.02])
        cb = mpl.colorbar.ColorbarBase(axcb, cmap=plt.get_cmap('Spectral'), norm=mpl.colors.Normalize(vmin=0., vmax=3),orientation='horizontal')
        return fig

    def make_contour_plot(self):
        # first version written by sabine, will add in the weekend
        #
        pass

    def plot_kinematic_maps(self, model=None, kin_set=0, cbar_lims='data'):
        """
        Show kinematic map of a model with v, sigma, h3, h4.
        If model=None, select best fitting model so far.
        Taken from schw_kin.py.
        Note: kin_set should be the index of the data set,
        e.g. kin_set=0 , kin_set=1. Default is kin_set=0.
        If kin_set=='all', kinematic maps for all kinematics are plotted
        and a list of (fig,kin_set_name) is returned where fig are figure
        objects and kin_set_name are the names of the kinematics sets.
        """
        stars = \
          self.system.get_component_from_class(physys.TriaxialVisibleComponent)
        n_kin = len(stars.kinematic_data)
        #########################################
        if kin_set == 'all':
            self.logger.info(f'Plotting kinematic maps for {n_kin} kin_sets.')
            figures = []
            for i in range(n_kin):
                fig = self.plot_kinematic_maps(model=model,
                                               kin_set=i,
                                               cbar_lims=cbar_lims)
                figures.append((fig, stars.kinematic_data[i].name))
            return figures # returns a list of (fig,kin_name) tuples
        #########################################
        if kin_set >= n_kin:
            text = f'kin_set must be < {n_kin}, but it is {kin_set}'
            self.logger.error(text)
            raise ValueError(text)
        self.logger.info(f'Plotting kinematic maps for kin_set no {kin_set}: '
                         f'{stars.kinematic_data[kin_set].name}')

        if model is None:
            which_chi2 = self.settings.parameter_space_settings['which_chi2']
            models_done = np.where(self.all_models.table['all_done'])
            min_chi2 = min(m[which_chi2]
                           for m in self.all_models.table[models_done])
            t = self.all_models.table.copy(copy_data=True) # deep copy!
            t.add_index(which_chi2)
            model_id = t.loc_indices[min_chi2]
            model = self.all_models.get_model_from_row(model_id)

        # currently this only works for GaussHermite's and LegacyWeightSolver
        kin_type = type(stars.kinematic_data[kin_set])
        if kin_type is not kinematics.GaussHermite:
            self.logger.info(f'kinematic maps cannot be plot for {kin_type} - '
                             'only GaussHermite')
            fig = plt.figure(figsize=(27, 12))
            return fig
        weight_solver = model.get_weights()
        ws_type = type(weight_solver)
        if ws_type is not weight_solvers.LegacyWeightSolver:
            self.logger.info('kinematic maps cannot be plot for weight solver '
                             f'{ws_type} - only LegacyWeightSolver')
            fig = plt.figure(figsize=(27, 12))
            return fig

        kinem_fname = model.get_model_directory() + 'nn_kinem.out'
        body_kinem = np.genfromtxt(kinem_fname, skip_header=1)

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

        text = '`cbar_lims` must be one of `model`, `data` or `combined`'
        # assert cbar_lims in ['model', 'data', 'combined'], text
        if not cbar_lims in ['model', 'data', 'combined']:
            self.logger.error(text)
            raise AssertionError(text)
        if cbar_lims=='model':
            vmax = np.max(np.abs(velm))
            smax, smin = np.max(sigm), np.min(sigm)
            h3max, h3min = np.max(h3m), np.min(h3m)
            h4max, h4min = np.max(h4m), np.min(h4m)
        elif cbar_lims=='data':
            vmax = np.max(np.abs(vel))
            smax, smin = np.max(sig), np.min(sig)
            h3max, h3min = np.max(h3), np.min(h3)
            h4max, h4min = np.max(h4), np.min(h4)
        elif cbar_lims=='combined':
            tmp = np.hstack((velm, vel))
            vmax = np.max(np.abs(tmp))
            tmp = np.hstack((sigm, sig))
            smax, smin = np.max(tmp), np.min(tmp)
            tmp = np.hstack((h3m, h3))
            h3max, h3min = np.max(tmp), np.min(tmp)
            tmp = np.hstack((h4m, h4))
            h4max, h4min = np.max(tmp), np.min(tmp)
        else:
            self.logger.error('unknown choice of `cbar_lims`')

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

        fig = plt.figure(figsize=(27, 12))
        plt.subplots_adjust(hspace=0.7,
                            wspace=0.01,
                            left=0.01,
                            bottom=0.05,
                            top=0.99,
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
