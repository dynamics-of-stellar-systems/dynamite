import os
import logging
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from plotbin.display_pixels import display_pixels
import cmasher as cmr
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

    The methodology in this class has been contributed by Ling Zhu and
    Giulia Santucci. Please cite Zhu+18, MNRAS 473, 3000 and
    Santucci+22, ApJ 930, 153 if used.

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
        set (e.g. kin_set=0, kin_set=1). The default is 0.
    ocut : list of floats, optional
        The cuts in lambda_z. The default is None, which translates to
        ocut=[0.8, 0.25, -0.25, -0.8], the selection in lambda_z
        following Santucci+22.
    decomp_table : bool, optional
        If True, write a table mapping each orbit to its respective
        component(s). The default is False.
    comps_weights : bool, optional
        If True, write a table of aggregated weights in each component.
        The default is False.

    Raises
    ------
    ValueError
        if no config object is given or the kin_set does not exist.

    """
    def __init__(self,
                 config=None,
                 model=None,
                 kin_set=0,
                 ocut=None,
                 decomp_table=False,
                 comps_weights=False):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None:
            text = f'{__class__.__name__} needs configuration object, ' \
                   'None provided.'
            self.logger.error(text)
            raise ValueError(text)
        self.config = config
        if model is None:
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
        self.orblib = self.model.get_orblib()
        self.orblib.read_losvd_histograms()
        self.losvd_histograms = self.orblib.losvd_histograms[self.kin_set]
        self.proj_mass = self.orblib.projected_masses[self.kin_set]
        self.logger.debug(f'{self.losvd_histograms.y.shape=}, '
                          f'{self.proj_mass.shape=}.')
        # Get orbit weights and store them in self.model.weights
        _ = self.model.get_weights(self.orblib)
        # Do the decomposition
        self.comps=['thin_d', 'thick_d', 'disk',
                    'cr_thin_d', 'cr_thick_d', 'cr_disk', 'bulge', 'all']
        # Important: the 'all' component needs to be the last one in the list!
        if ocut is not None:
            self.ocut = ocut
        else:
            self.ocut = [  0.8,     0.25,   -0.25,        -0.8        ]
        #             thin_d  thick_d   bulge    cr_thick_d   cr_thin_d
        self.decomp = self.decompose_orbits()
        self.logger.info('Orbits read and velocity histogram created.')
        if decomp_table:
            file_name = self.model.directory + 'decomp_table.ecsv'
            self.decomp.write(file_name, format='ascii.ecsv', overwrite=True)
            self.logger.info('Orbit decomposition information written to '
                             f'{file_name}.')
        if comps_weights:
            self.comps_weights()

    def plot_decomp(self,
                    xlim,
                    ylim,
                    v_sigma_option='fit',
                    comps_plot='all',
                    individual_colorbars=False,
                    figtype='.png',
                    dpi=100):
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
        comps_plot : dict or string 'all', optional
            If 'all', all components will be in the decomposition plot.
            Specific components can be selected by passing a dictionary, e.g.,
            comps_plot = {'thin_d': True, 'thick_d': True, 'disk': True,
                          'cr_thin_d': False, 'cr_thick_d': False,
                          'cr_disk: False', 'bulge': False, 'all': False} will
            only create the plots for 'thin_d', 'thick_d', and 'disk'. `False`
            entries can be omitted in the dictionary. The default is 'all'.
        individual_colorbars : bool or dict, optional
            If True, then the sb (surface brightness), vel (velocity), and
            sig (velocity dispersion) colorbars adapt to their respective
            value ranges. This can be useful for identifying structures
            invisible otherwise.
            If False, the sb, vel, and sig colorbars will be the same for all
            components.
            The individual colorbars are accessed by passing a dict. For
            example: {'sb': True, 'vel': False, 'sig': False} will only
            adapt the sb colorbar to the respective component's value range.
            'False' entries can be omitted.
            The default is individual_colorbars=False.
        figtype : str, optional
            Determines the file format and extension to use when saving the
            figure. The default is '.png'.
        dpi : float, optional
            The resolution of saved figures (if not overridden later). The
            default is 100 dpi.

        Returns
        -------
        None.

        """
        mpl.rcParams['savefig.dpi'] = dpi
        comp_kinem_moments = self.comps_aphist(v_sigma_option)
        self.logger.info('Component data done. '
                         f'Plotting decomposition for {v_sigma_option=}.')
        weights = self.model.weights
        comps = self.decomp.meta["comps"]

        if comps_plot == 'all':
            comps_plot = {comp: True for comp in comps}
        for comp in comps:
            if comp not in comps_plot:
                comps_plot[comp] = False
        self.logger.info(f'Plotting data for components {comps_plot}.')

        if type(individual_colorbars) is bool:
            switch = individual_colorbars
            individual_colorbars = {k: switch for k in ['sb', 'vel', 'sig']}
        else:
            for k in ['sb', 'vel', 'sig']:
                if k not in individual_colorbars:
                    individual_colorbars[k] = False

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

        vel = []
        sig = []
        t = []
        min_flux = {}
        max_flux = {}
        max_vel = {}
        min_sig = {}
        max_sig = {}
        last_comps_idx = len(comps) - 1
        for c_idx, comp in enumerate(comps):
            labels = [col for col in comp_kinem_moments.colnames
                          if col.startswith(comp)]
            flux = comp_kinem_moments[labels[0]]
            w = weights[[f'|{comp}|' in s for s in self.decomp['component']]]
            fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux))
            flux = flux / fhist
            tt = flux[grid]*1.
            tt = tt * np.sum(w)/np.sum(tt)
            t.append(tt.copy())
            vel.append(comp_kinem_moments[labels[1]])
            sig.append(comp_kinem_moments[labels[2]])
            min_flux[comp] = np.nanmin(np.log10(tt[tt != 0]))
            max_flux[comp] = np.nanmax(np.log10(tt[tt != 0]))
            max_vel[comp] = max(np.nanmax(vel[c_idx]), -np.nanmin(vel[c_idx]))
            min_sig[comp] = np.nanmin(np.array(sig[c_idx])[np.array(sig[c_idx]) > 0])
            max_sig[comp] = np.nanmax(np.array(sig[c_idx])[np.array(sig[c_idx]) > 0])
            if c_idx == last_comps_idx:  # the last item MUST be 'all'
                # minf = np.nanmin(np.log10(tt))
                # maxf = np.nanmax(np.log10(tt[tt !=0]))
                # maxv = max_vel[comp]
                # mins = min_sig[comp]
                # maxs = max_sig[comp]
                totalf = np.sum(tt)  # tt here refers to the 'all' comp

        t = t/totalf
        for comp in comps:
            if not individual_colorbars['sb']:
                min_flux[comp] = min(min_flux[c] for c in comps)
                max_flux[comp] = max(max_flux[c] for c in comps)
            if not individual_colorbars['vel']:
                max_vel[comp] = max(max_vel[c] for c in comps)
            if not individual_colorbars['sig']:
                min_sig[comp] = min(min_sig[c] for c in comps)
                max_sig[comp] = max(max_sig[c] for c in comps)

        # if not individual_sb_colorbars:
        #     for comp in comps:
        #         min_flux[comp] = minf
        #         max_flux[comp] = maxf

        # vmax = np.nanmax(vel)
        # sig_t = np.array(sig)

        # smax = np.nanmax(sig_t[sig_t > 0])
        # smin = np.nanmin(sig_t[sig_t > 0])

        xi_t=(xi[s])
        yi_t=(yi[s])

        table = {'x/arcs':xi_t,'y/arcs':yi_t}
        for c_idx, comp in enumerate(comps):
            labels = [col for col in comp_kinem_moments.colnames
                          if col.startswith(comp)]
            table.update({labels[0]:t[c_idx][s],
                         labels[1]:vel[c_idx][grid[s]],
                         labels[2]:sig[c_idx][grid[s]]})
        comps_kin = astropy.table.Table(table)

        kin_name = stars.kinematic_data[self.kin_set].name
        file_name = f'comps_kin_{v_sigma_option}_{kin_name}'
        table_file_name = self.model.directory + file_name + '.ecsv'
        plot_file_name = self.config.settings.io_settings['plot_directory'] \
                         + file_name \
                         + figtype
        comps_kin.write(table_file_name, format='ascii.ecsv', overwrite=True)
        self.logger.info('Component grid kinematics written to '
                         f'{table_file_name}.')

        self.logger.debug(f'{v_sigma_option}: {min_flux=}, {max_flux=}, '
                          f'{max_vel=}, {min_sig=}, {max_sig=}.')

        c_skipped = len([comp for comp in comps_plot if not comps_plot[comp]])
        LL = len(comps) - c_skipped
        map1 = cmr.get_sub_cmap('twilight_shifted', 0.05, 0.6)
        map2 = cmr.get_sub_cmap('twilight_shifted', 0.05, 0.95)
        # titles = ['THIN DISK','THICK DISK','DISK','BULGE','ALL']
        # compon = np.array(['thin_d','thick_d','disk','bulge','all'])
        titles = [c.replace('_disk', ' disk').replace('_d', ' disk').replace('_',' ')
                  for c in comps]
        compon = np.array(comps)
        kwtext = dict(size=20, ha='center', va='center', rotation=90.)
        kw_display1 = dict(pixelsize=dx, colorbar=True,
                                  nticks=7, cmap=map1)
        kw_display2 = dict(pixelsize=dx, colorbar=True,
                                  nticks=7, cmap=map2)

        plt.figure(figsize=(16, int((LL+2)*3)*ylim/xlim))
        plt.subplots_adjust(hspace=0.7, wspace=0.01, left=0.01,
                            bottom=0.05, top=0.99, right=0.99)

        i_plot = 0
        for c_idx, comp in enumerate(comps):
            if not comps_plot[comp]:
                continue
            ax = plt.subplot(LL, 3, 3*i_plot+1)
            if i_plot == 0:
                ax.set_title('surface brightness (log)',fontsize=20,pad=20)
            display_pixels(xi_t,
                           yi_t,
                           np.log10(t[c_idx][s])-max_flux[comp],
                           vmin=min_flux[comp]-max_flux[comp],
                           vmax=0,
                           **kw_display1)
            ax.text(-0.32, 0.5, titles[np.where(compon==comp)[0][0]],
                    **kwtext, transform=ax.transAxes)
            ax.set_ylabel('arcsec')
            if i_plot == LL - 1:
                ax.set_xlabel('arcsec')

            ax = plt.subplot(LL, 3, 3*i_plot+2)
            if i_plot == 0:
                plt.title('velocity',fontsize=20,pad=20)
            display_pixels(xi_t, yi_t, vel[c_idx][grid[s]],
                           vmin=-1.0*max_vel[comp], vmax=max_vel[comp],
                           **kw_display2)
            if i_plot == LL - 1:
                ax.set_xlabel('arcsec')

            ax = plt.subplot(LL, 3, 3*i_plot+3)
            if i_plot == 0:
                plt.title('velocity dispersion',fontsize=20,pad=20)
            display_pixels(xi_t, yi_t, sig[c_idx][grid[s]],
                           vmin=min_sig[comp], vmax=max_sig[comp],
                           **kw_display1)
            if i_plot == LL - 1:
                ax.set_xlabel('arcsec')
            i_plot += 1

        plt.tight_layout()
        plt.savefig(plot_file_name)
        self.logger.info(f'Component plots written to {plot_file_name}.')
        plt.close()

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
            self.logger.info(f'Component {comp}...')
            # calculate flux and losvd histograms for component
            orb_sel = np.array([f'|{comp}|' in s for s in self.decomp['component']],
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
            # Important: the sequence of lsb - v - sig matters!
        return comp_flux_v_sigma

    def decompose_orbits(self, ocut=None):
        """Decompose orbits based on lambda_z.

        Parameters
        ----------
        ocut : DEPRECATED, will be removed in the next major release.
            Use ocut= when instatiating the Decomposition object.'

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
        comps = self.comps
        if ocut is None:
            ocut = self.ocut
        else:
            self.logger.warning('Argument ocut is DEPRECATED and will be '
                                'removed in the next major release. Use '
                                'ocut= when instatiating the '
                                f'{__class__.__name__} object.')
        self.logger.debug(f'Cut lines are: {ocut}.')
        file2 = self.model.directory_noml + 'datfil/orblib.dat_orbclass.out'
        file3 = self.model.directory_noml + 'datfil/orblibbox.dat_orbclass.out'
        file3_test = os.path.isfile(file3)
        if not file3_test:
            file3 = '%s' % file2

        n_orb = self.config.settings.orblib_settings['nE'] * \
                self.config.settings.orblib_settings['nI2'] * \
                self.config.settings.orblib_settings['nI3']
        n_dither = self.config.settings.orblib_settings['dithering']
        conversion_factor=self.config.all_models.system.distMPc*1.0e6*1.49598e8

        ncol = n_dither ** 3
        orbclass1=self.orblib.read_orbit_property_file_base(file2, ncol, n_orb)
        orbclass2=self.orblib.read_orbit_property_file_base(file3, ncol, n_orb)

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
        lz = (orbclass[2, :, :] / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))

        # Average values for the orbits in the same bundle (n_dither^3).
        # Only include the orbits within Rmax_arcs
        rm = np.sum(orbclass[3, :, :]/conversion_factor, axis=0) / n_dither**3

        # flip the sign of lz to confirm total(lz) > 0
        t = np.ravel(np.argsort(rm))

        yy = np.max(np.ravel(np.where(np.cumsum(orbw[t]) <= 0.5)))
        k = t[0:yy]
        if np.sum(np.sum(lz[:, k], axis=0) / (n_dither ** 3) * orbw[k]) < 0:
            lz *= -1.0

        lzm_sign= np.sum(lz, axis=0) / n_dither ** 3

        self.logger.info(f'Decomposing {n_orbs} orbits into {comps=}...')
        decomp = astropy.table.Table({'id':range(n_orbs),
                                      'component':['']*n_orbs},
                                     dtype=[int, 'U256'],
                                     meta={'comps':comps})
        # map components
        comp_map = np.zeros(n_orbs, dtype=int)
        # cold component (thin disk)
        comp_map[np.ravel(np.where(lzm_sign >= ocut[0]))] += \
            2**comps.index('thin_d')
        # warm component (thick disk)
        comp_map[np.ravel(np.where((lzm_sign > ocut[1])
                                 & (lzm_sign < ocut[0])))] += \
            2**comps.index('thick_d')
        # hot component (bulge)
        comp_map[np.ravel(np.where((lzm_sign > ocut[2])
                                 & (lzm_sign <= ocut[1])))] += \
            2**comps.index('bulge') # was lzm_sign<ocut[1]
        # disk component (disk)
        comp_map[np.ravel(np.where(lzm_sign > ocut[1]))] += \
            2**comps.index('disk')
        # counter-rotating cold component (cr thin disk)
        comp_map[np.ravel(np.where(lzm_sign <= ocut[3]))] += \
            2**comps.index('cr_thin_d')
        # counter-rotating warm component (cr thick disk)
        comp_map[np.ravel(np.where((lzm_sign > ocut[3])
                                 & (lzm_sign <= ocut[2])))] += \
            2**comps.index('cr_thick_d')
        # counter-rotating disk (cr disk)
        comp_map[np.ravel(np.where((lzm_sign <= ocut[2])))] += \
            2**comps.index('cr_disk')
        # whole component (all)
        comp_map += 2**comps.index('all')
        for i in np.ravel(np.where(comp_map > 0)):
            for k, comp in enumerate(comps):
                if comp_map[i] & (1 << k):
                    decomp['component'][i] += f'|{comp}|'
        return decomp

    def comps_weights(self):
        """ Write a table of aggregated weights in each component.
        """
        weights = self.model.weights
        comps = self.decomp.meta["comps"]
        comps_weights = []
        for comp in comps:
            w = weights[[f'|{comp}|' in s for s in self.decomp['component']]]
            comps_weights.append(sum(w))
        weights_table = astropy.table.Table([comps, comps_weights],
                                            names=('component', 'weight'),
                                            dtype=['U256', float])
        file_name = self.model.directory + 'comps_weights.ecsv'
        weights_table.write(file_name, format='ascii.ecsv', overwrite=True)
        self.logger.info(f'Component aggregate weights written to {file_name}.')



class Analysis:
    """Class to hold results' analysis methods.

    This class contains methods that help analyzing DYANMITE results and can
    be called, e.g. by plotting routines.

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    model : a ``dyn.model.Model`` object, optional, default: best model so far
    kin_set : int, optional, default: 0

    """

    def __init__(self, config, model=None, kin_set=0):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None:
            text = f'{__class__.__name__} needs configuration object, ' \
                   'None provided.'
            self.logger.error(text)
            raise ValueError(text)
        if model is None:
            model_index = config.all_models.get_best_n_models_idx(1)[0]
            model = config.all_models.get_model_from_row(model_index)
        self.config = config
        self.model = model
        self.kin_set = kin_set

    def get_gh_model_kinematic_maps(self,
                                    model=None,
                                    kin_set=None,
                                    v_sigma_option='fit',
                                    kinematics_as='table',
                                    weights=None):
        """
        Generates an astropy table in the model directory that holds the
        model's data for creating Gauss-Hermite kinematic maps:
        flux, v, sigma, h3 ... h<number_GH>.
        v and sigma are either directly calculated from the model's losvd
        histograms or from fitting a Gaussian in each aperture.

        Parameters
        ----------
        model : a ``dyn.model.Model`` object, optional
            The default is the Analysis object's model.
        kin_set : int, optional
            Which kinematics set to use. The default is the
            Analysis object's kin_set.
        v_sigma_option : str, optional
            If 'fit', v_mean and v_sigma are calculated based on fitting
            Gaussians, if 'moments', v_mean and v_sigma are calculated
            directly from the model's losvd histograms. The default is 'fit'.
        kinematics_as : str, optional
            If 'table', return ``gh_table``, the model's kinematics as an
            astropy table, if 'file', write the table to disk in ascii.ecsv
            format and return its full path ``f_name``, if 'both', write the
            table to disk and return a tuple ``(gh_table, f_name)``.
            The default is 'table'.
        weights : ``numpy.array`` like, optional
            Orbital weights to use. The default is ``None`` and will
            determine the weights via ``model.get_weights(orblib)``.

        Raises
        ------
        ValueError
            if v_sigma_option or kinematics_as are invalid.

        Returns
        -------
        gh_table : astropy table (if kinematics_as='table')
            The astropy table holding the model's gh kinematics.
        f_name : str (if kinematics_as='file')
            The file name (full path) of the astropy table holding the model's
            gh kinematics.
        (gh_table, f_name) : tuple  (if kinematics_as='both')

        """
        if model is None:
            model = self.model
        if kin_set is None:
            kin_set = self.kin_set
        if v_sigma_option not in ['moments', 'fit']:
            txt = f"{v_sigma_option=} but must be either 'fit' or 'moments'."
            self.logger.error(txt)
            raise ValueError(txt)
        if kinematics_as not in ['table', 'file', 'both']:
            txt = f"{kinematics_as=} but must be either 'table', 'file', or " \
                   "'both'."
            self.logger.error(txt)
            raise ValueError(txt)
        stars = self.config.system.get_component_from_class(
                                dyn.physical_system.TriaxialVisibleComponent)
        kin_name = stars.kinematic_data[kin_set].name
        self.logger.info('Getting projected masses and losvds for '
                         f'model {model.directory}.')
        orblib = model.get_orblib()
        if weights is None:
            _ = model.get_weights(orblib)
            weights = model.weights
        # get losvd_histograms and projected masses:
        orblib.read_losvd_histograms()
        # get all orbits' losvds; orbits_losvd.shape = n_orb,n_vbin,n_aperture
        orbits_losvd = orblib.losvd_histograms[kin_set].y[:,:,]
        # weighted sum of orbits_losvd; model_losvd.shape = 1,n_vbin,n_aperture
        model_losvd = np.dot(orbits_losvd.T, weights).T[np.newaxis]
        #model_losvd /= np.sum(model_losvd, 0) # normalisation not necessary
        model_proj_masses = np.dot(orblib.projected_masses[kin_set].T,
                                   weights) # .shape = n_aperture
        # calculate v_mean and v_sigma values from the losvd histograms
        model_losvd_hist = \
            dyn.kinematics.Histogram(xedg=orblib.losvd_histograms[kin_set].xedg,
                                     y=model_losvd,
                                     normalise=False)
        if v_sigma_option == 'moments':
            v_mean = np.squeeze(model_losvd_hist.get_mean()) # from distr.
            v_sigma = np.squeeze(model_losvd_hist.get_sigma()) # from distr.
            v_sig_text = 'losvd moments'
        elif v_sigma_option == 'fit':
            v_mean, v_sigma = model_losvd_hist.get_mean_sigma_gaussfit()
            v_mean = np.squeeze(v_mean)
            v_sigma = np.squeeze(v_sigma)
            v_sig_text = 'fitted Gaussians'
        else:
            pass
        self.logger.debug(f'Calculated v_mean and v_sigma from {v_sig_text} '
                          f'for {len(v_mean)} apertures.')
        gh_table = astropy.table.Table([model_proj_masses, v_mean, v_sigma],
                                       names = ['flux', 'v', 'sigma'],
                                       meta={'v_sigma_option': v_sigma_option,
                                             'kin_set': kin_name})
        weight_solver_settings = self.config.settings.weight_solver_settings
        n_gh = weight_solver_settings['number_GH']
        if n_gh > 2:
            # calculate the model's gh expansion coefficients
            gh = dyn.kinematics.GaussHermite()
            gh.data = gh_table
            model_gh_coefficients = gh.transform_orblib_to_observables(
                            losvd_histograms=model_losvd_hist,
                            weight_solver_settings=weight_solver_settings)
            # unscale by projected masses (see weight solver)
            model_gh_coefficients = \
                (np.squeeze(model_gh_coefficients).T / model_proj_masses).T
            # add the gh coefficients to the astropy table
            col_names = [f'h{i}' for i in range(3,n_gh+1)]
            tab_data = list(model_gh_coefficients[:,2:].T)
            gh_table.add_columns(tab_data, names=col_names)
        if kinematics_as == 'table':
            return gh_table
        f_name = f'{model.directory}model_gh_kins_' + \
                 f'{kin_name}_{v_sigma_option}.ecsv'
        gh_table.write(f_name, format='ascii.ecsv', overwrite=True)
        self.logger.info(f'Model gh kinematics {kin_name}, {n_gh=} '
                         f'written to {f_name}.')
        if kinematics_as == 'file':
            return f_name
        return (gh_table, f_name)
