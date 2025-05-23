import os
import logging
import numpy as np
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
        self.plot_comps(xlim=xlim, ylim=ylim,
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
            self.logger.info(f'Component {comp}...')
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

        comps=['thin_d', 'warm_d', 'disk', 'bulge', 'all']
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

    def plot_comps(self,
                   xlim,
                   ylim,
                   comp_kinem_moments,
                   figtype='.png'):
        """ Generate decomposition plots.

        Parameters
        ----------
        xlim : float
            restricts plot x-coordinates to abs(x) <= xlim.
        ylim : float
            restricts plot y-coordinates to abs(y) <= ylim.
        comp_kinem_moments : astropy table
            The table columns are: aperture index (starting with 0), followed
            by three columns per component holding the flux, mean velocity,
            and velocity dispersion.
            The chosen v_sigma_option is in the table meta data.
        figtype : str, optional
            Determines the file format and extension to use when saving the
            figure. The default is '.png'.

        Returns
        -------
        None.

        """

        v_sigma_option = comp_kinem_moments.meta['v_sigma_option'] \
                         if 'v_sigma_option' in comp_kinem_moments.meta.keys()\
                         else ''
        self.logger.info(f'Plotting decomposition for {v_sigma_option=}.')

        weights = self.model.weights
        comps = self.decomp.meta["comps"]

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

        quant = ['_lsb', '_v', '_sig']
        vel = []
        sig = []
        t = []
        totalf = 0
        for i in range(len(comps)):
            labels = [comps[i] + qq for qq in quant]
            flux = comp_kinem_moments[labels[0]]
            w = weights[[comps[i] in s for s in self.decomp['component']]]
            fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux))
            flux = flux / fhist
            tt = flux[grid]*1.
            tt = tt * np.sum(w)/np.sum(tt)
            t.append(tt.copy())
            if comps[i] in ['thin_d', 'warm_d', 'bulge']:
                totalf += np.sum(tt)
                if comps[i] == 'thin_d':
                    fluxtot = tt
                else:
                    fluxtot += tt
            vel.append(comp_kinem_moments[labels[1]])
            sig.append(comp_kinem_moments[labels[2]])

        t = t/totalf

        vmax = np.nanmax(vel)
        sig_t = np.array(sig)

        smax = np.nanmax(sig_t[sig_t > 0])
        smin = np.nanmin(sig_t[sig_t > 0])

        minf=np.nanmin(np.log10(fluxtot))
        maxf=np.nanmax(np.log10(fluxtot[fluxtot !=0]))
        xi_t=(xi[s])
        yi_t=(yi[s])

        table = {'x/arcs':xi_t,'y/arcs':yi_t}
        for i in range(len(comps)):
            labels = [comps[i] + qq for qq in quant]
            table.update({labels[0]:t[i][s],
                         labels[1]:vel[i][grid[s]],
                         labels[2]:sig[i][grid[s]]})
        comps_kin = astropy.table.Table(table)

        kin_name = stars.kinematic_data[self.kin_set].name
        file_name = f'comps_kin_{v_sigma_option}_{kin_name}'
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

        LL = len(comps)
        map1 = cmr.get_sub_cmap('twilight_shifted', 0.05, 0.6)
        map2 = cmr.get_sub_cmap('twilight_shifted', 0.05, 0.95)
        titles = ['THIN DISK','THICK DISK','DISK','BULGE','ALL']
        compon = np.array(['thin_d','warm_d','disk','bulge','all'])
        kwtext = dict(size=20, ha='center', va='center', rotation=90.)
        kw_display1 = dict(pixelsize=dx, colorbar=True,
                                  nticks=7, cmap=map1)
        kw_display2 = dict(pixelsize=dx, colorbar=True,
                                  nticks=7, cmap=map2)

        plt.figure(figsize=(16, int((LL+2)*3)*ylim/xlim))
        plt.subplots_adjust(hspace=0.7, wspace=0.01, left=0.01,
                            bottom=0.05, top=0.99, right=0.99)

        for ii in range(len(comps)):
            ax = plt.subplot(LL, 3, 3*ii+1)
            if ii == 0:
                ax.set_title('surface brightness (log)',fontsize=20,pad=20)
            display_pixels(xi_t, yi_t, np.log10(t[ii][s])-maxf,
                           vmin=minf-maxf, vmax=0, **kw_display1)
            ax.text(-0.2, 0.5, titles[np.where(compon==comps[ii])[0][0]],
                    **kwtext, transform=ax.transAxes)

            plt.subplot(LL, 3, 3*ii+2)
            if ii == 0:
                plt.title('velocity',fontsize=20,pad=20)
            display_pixels(xi_t, yi_t, vel[ii][grid[s]],
                           vmin=-1.0*vmax, vmax=vmax, **kw_display2)

            plt.subplot(LL, 3, 3*ii+3)
            if ii == 0:
                plt.title('velocity dispersion',fontsize=20,pad=20)
            display_pixels(xi_t, yi_t, sig[ii][grid[s]],
                           vmin=smin, vmax=smax, **kw_display1)

        plt.tight_layout()
        plt.savefig(plot_file_name)
        self.logger.info(f'Component plots written to {plot_file_name}.')
        plt.close()


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

    def get_flux_for_selected_orbit_bundles(self,
                                            model=None,
                                            kin_set=None,
                                            weights=None,
                                            bundle_mapping=None,
                                            flux_as='table'):
        """
        Generates an astropy table that holds the weighted contribution of the
        orbit bundles defined in bundle_mapping to the model's projected mass
        in each aperture.

        Parameters
        ----------
        model : a ``dyn.model.Model`` object, optional
            The default is the Analysis object's model.
        kin_set : int, optional
            Which kinematics set to use for the spatial binning. The default
            is the Analysis object's kin_set.
        weights : ``numpy.array`` like, shape=(n_orbits,), optional
            Orbital weights to use. The default is ``None`` and will
            determine the weights via ``model.get_weights(orblib)``.
        bundle_mapping : ``numpy.array`` like, shape=(n_bundles, n_orbits),
                         mandatory
            The mapping of orbits to orbit bundles. The values in the array
            indicate what weight a specific orbit contributes to a specific
            bundle. Note that an orbit may contribute to more than one bundle
            because - depending on the ``dithering`` configuration setting -
            each orbit may in fact represent a family of similar orbits (not
            to be mistaken for the bundles in the sense of this method).
        flux_as : str, optional
            If 'table', return ``map_table``, if 'file', write the table to
            the model directory in ascii.ecsv format and return its full path
            ``f_name``, if 'both', write the table to disk and return a tuple
            ``(gh_table, f_name)``. The default is 'table'.

        Raises
        ------
        ValueError
            if flux_as is invalid.

        Returns
        -------
        map_table : astropy table (if flux_as='table')
            The astropy table holding the model's weighted contribution of the
            orbit bundles defined in bundle_mapping to the model's projected
            mass in each aperture. The table columns are flux_000, ...,
            flux_n_b, flux_all, where flux_xxx is the weighted contribution of
            orbit bundle xxx, flux_all the weighted contribution of all orbit
            bundles, and n_b is the number of orbit bundles.
        f_name : str (if flux_as='file')
            The file name (full path) of the astropy table holding the model's
            gh kinematics.
        (gh_table, f_name) : tuple  (if flux_as='both')

        """
        if flux_as not in ['table', 'file', 'both']:
            txt = f"{flux_as=} but must be either 'table', 'file', or 'both'."
            self.logger.error(txt)
            raise ValueError(txt)
        if model is None:
            model = self.model
        if kin_set is None:
            kin_set = self.kin_set
        stars = self.config.system.get_unique_triaxial_visible_component()
        kin_name = stars.kinematic_data[kin_set].name
        self.logger.info('Getting model projected masses.')
        orblib = model.get_orblib()
        if weights is None:
            _ = model.get_weights(orblib)
            weights = model.weights
        # Get projected masses if necessary
        if not hasattr(orblib, 'projected_masses'):
            orblib.read_losvd_histograms()
        # flux.shape = (n_bundles, n_aperture)
        flux = np.matmul(bundle_mapping, orblib.projected_masses[kin_set])
        flux_all = np.sum(flux, axis=0)
        n_bins = flux.shape[0]
        map_table = astropy.table.Table(
            np.hstack((flux.T, flux_all[:,np.newaxis])),
            names = [f'flux_{i:03d}' for i in range(n_bins)] + ['flux_all'],
            meta={'kin_set': kin_name})
        if flux_as == 'table':
            return map_table
        f_name = f'{model.directory}orbit_bundle_flux_{kin_name}.ecsv'
        map_table.write(f_name, format='ascii.ecsv', overwrite=True)
        self.logger.info('Flux for orbit bundles binned for kinematics '
                         f'{kin_name} written to {f_name}.')
        if flux_as == 'file':
            return f_name
        return (map_table, f_name)


    def get_gh_model_kinematic_maps(self,
                                    model=None,
                                    kin_set=None,
                                    v_sigma_option='fit',
                                    kinematics_as='table',
                                    weights=None):
        """
        Generates an astropy table that holds the
        model's data for creating Gauss-Hermite kinematic maps:
        v, sigma, h3 ... h<number_GH>.
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
            astropy table, if 'file', write the table to the model directory in
            ascii.ecsv format and return its full path ``f_name``, if 'both',
            write the table to disk and return a tuple ``(gh_table, f_name)``.
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
        self.logger.info('Getting model projected masses and losvds.')
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
