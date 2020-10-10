import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from plotbin import sauron_colormap as pb_sauron_colormap
from plotbin import cap_display_pixels

class Plotter(object):

    def __init__(self,
                 system=None,
                 settings=None,
                 parspace=None,
                 all_models=None):
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
        # <-- dont need to implement "which plot to make" switches yet, just indeividual plotting methods

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
        pass

    def make_chi2_plot(self):
        # ...
        #
        pass

    def plot_kinematic_maps(self, model):
        """
        Show kinematic map of the best-fitting model, with v, sigma, h3, h4
        Taken from schw_kin.py.
        """
        kinem_fname = model.get_model_directory() + 'nn_kinem.out'
        body_kinem = np.genfromtxt(kinem_fname, skip_header=1)
        id, fluxm, flux, velm, vel, dvel, sigm, sig, dsig, h3m, h3, dh3, h4m, h4, dh4 = body_kinem.T

        vmax = max(np.hstack((velm, vel)))
        smax = max(np.hstack((sigm, sig)))
        smin = min(sig)
        h3max = max(np.hstack((h3m, h3)));
        h3min = min(np.hstack((h3m, h3)))
        h4max = max(np.hstack((h4m, h4)));
        h4min = min(np.hstack((h4m, h4)))

        # Read aperture.dat
        # The angle that is saved in this file is measured counter clock-wise
        # from the galaxy major axis to the X-axis of the input data.
        stars = self.system.get_component_from_name('stars')
        aperture_fname = stars.kinematic_data[0].aperturefile
        aperture_fname = self.input_directory + aperture_fname

        lines = [line.rstrip('\n').split() for line in open(aperture_fname)]
        strhead = lines[0]
        minx = np.float(lines[1][0]);
        miny = np.float(lines[1][1])
        sx = np.float(lines[2][0]);
        sy = np.float(lines[2][1])
        maxx = sx + minx;
        sy = sy + miny
        angle_deg = np.float(lines[3][0])  # - 90.0 + 180
        nx = np.int(lines[4][0]);
        ny = np.int(lines[4][1])
        dx = sx / nx

        # print("Pixel grid dimension is dx,nx,ny,", dx, nx, ny)
        grid = np.zeros((nx, ny), dtype=int)

        xr = np.arange(nx, dtype=float) * dx + minx + 0.5 * dx
        yc = np.arange(ny, dtype=float) * dx + miny + 0.5 * dx

        xi = np.outer(xr, (yc * 0 + 1));
        xt = xi.T.flatten()
        yi = np.outer((xr * 0 + 1), yc);
        yt = yi.T.flatten()

        radeg = 57.2958
        #print('PA: ', angle_deg)
        xi = np.cos(angle_deg / radeg) * xt - np.sin(angle_deg / radeg) * yt
        yi = np.sin(angle_deg / radeg) * xt + np.cos(angle_deg / radeg) * yt

        # read bins.dat
        stars = self.system.get_component_from_name('stars')
        bin_fname = stars.kinematic_data[0].binfile
        bin_fname = self.input_directory + bin_fname
        lines_bins = [line.rstrip('\n').split() for line in open(bin_fname)]
        i = 0;
        str_head = [];
        i_var = [];
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

        # The galaxy has already rotated with PA to make major axis aligned with x
        angle_deg = 0
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
        cap_display_pixels.display_pixels(x, y, c,
                                          vmin=minf, vmax=maxf,
                                          label='log10(SB)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 2)
        cap_display_pixels.display_pixels(x, y, vel[grid[s]],
                                          vmin=-1.0 * vmax, vmax=vmax,
                                          label='Velocity',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 3)
        cap_display_pixels.display_pixels(x, y, sig[grid[s]],
                                          vmin=smin, vmax=smax,
                                          label='Sigma',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 4)
        cap_display_pixels.display_pixels(x, y, h3[grid[s]],
                                          vmin=h3min, vmax=h3max,
                                          label=r'$h_{3}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 5)
        cap_display_pixels.display_pixels(x, y, h4[grid[s]],
                                          vmin=h4min, vmax=h4max,
                                          label=r'$h_{4}$',
                                          **kw_display_pixels)

        ### PLOT THE MODEL DATA
        plt.subplot(3, 5, 6)
        c = np.array(list(map(np.log10, fluxm[grid[s]] / max(fluxm))))
        cap_display_pixels.display_pixels(x, y, c,
                                          vmin=minfm, vmax=maxfm,
                                          label='log10(SB) (model)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 7)
        cap_display_pixels.display_pixels(x, y, velm[grid[s]],
                                          vmin=-1.0 * vmax, vmax=vmax,
                                          label='Velocity (model)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 8)
        cap_display_pixels.display_pixels(x, y, sigm[grid[s]],
                                          vmin=smin, vmax=smax,
                                          label='Sigma (model)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 9)
        cap_display_pixels.display_pixels(x, y, h3m[grid[s]],
                                          vmin=h3min, vmax=h3max,
                                          label=r'$h_{3}$'+' (model)',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 10)
        cap_display_pixels.display_pixels(x, y, h4m[grid[s]],
                                          vmin=h4min, vmax=h4max,
                                          label=r'$h_{4}$'+' (model)',
                                          **kw_display_pixels)


        ### PLOT THE ERROR NORMALISED RESIDUALS
        plt.subplot(3, 5, 11)
        c = (fluxm[grid[s]] - flux[grid[s]]) / flux[grid[s]]
        cap_display_pixels.display_pixels(x, y, c,
                                          vmin=-0.05, vmax=0.05,
                                          label=r'$\delta_{flux}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 12)
        c = (velm[grid[s]] - vel[grid[s]]) / dvel[grid[s]]
        cap_display_pixels.display_pixels(x, y, c,
                                          vmin=-10, vmax=10,
                                          label=r'$\delta_{vel}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 13)
        c = (sigm[grid[s]] - sig[grid[s]]) / dsig[grid[s]]
        cap_display_pixels.display_pixels(x, y, c,
                                          vmin=-10, vmax=10,
                                          label=r'$\delta_{sigma}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 14)
        c = (h3m[grid[s]] - h3[grid[s]]) / dh3[grid[s]]
        cap_display_pixels.display_pixels(x, y, c,
                                          vmin=-1, vmax=1,
                                          label=r'$\delta_{h_3}$',
                                          **kw_display_pixels)
        plt.subplot(3, 5, 15)
        c = (h4m[grid[s]] - h4[grid[s]]) / dh4[grid[s]]
        cap_display_pixels.display_pixels(x, y, c,
                                          vmin=-1, vmax=1,
                                          label=r'$\delta_{h_4}$',
                                          **kw_display_pixels)
        fig.subplots_adjust(left=0.07, wspace=0.3, hspace=0.01)
        kwtext = dict(size=20, ha='center', va='center', rotation=90.)
        fig.text(0.05, 0.9, 'Data', **kwtext)
        fig.text(0.05, 0.53, 'Model', **kwtext)
        fig.text(0.05, 0.2, 'Residual', **kwtext)
        return fig















    # etc...
