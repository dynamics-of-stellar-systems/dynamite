
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy

import yaml
import pymc as pm

import cmasher

from vorbin.voronoi_2d_binning import voronoi_2d_binning

from dynamite import constants


class Coloring:
    """Class to hold coloring-related routines

    This class provides methods for Voronoi binning of orbits in the
    radius-circularity phase space, fitting Bayesian models to the
    observed data, and calculating orbital decomposition of the coloring
    data. It also includes methods for plotting the results of the Voronoi
    binning and the orbital decomposition.

    The methods and code in this DYNAMITE class are inspired and partly
    adapted from Zhu et al., 2020, MNRAS, 496, 1579 and Zhu et al., 2022,
    A&A, 664, A115. Many thanks to Ling Zhu for sharing her code, which was
    instrumental for developing this class.
    The Voronoi orbit binning makes use of the vorbin package by Michele
    Cappellari (M. Cappellari, Y. Copin, 2003, MNRAS 342, 345).

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
        the configuration object containing the model and system settings.
    minr : float or str, optional
        the minimum radius [kpc] considered in the binning. If 'auto', minr
        is set to zero. The default is 'auto'.
    maxr : float or str, optional
        the maximum radius [kpc] considered in the binning. If 'auto', maxr
        is set to the config file's ``orblib_settings: logrmax`` value,
        converted to kpc. The default is 'auto'.
    nr : int, optional
        number of radial bins. The default is 50.
    nl : int, optional
        number of circularity bins. The default is 61.

    """

    def __init__(self, config, minr='auto', maxr='auto', nr=50, nl=61):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.config = config
        arc_kpc = constants.ARC_KPC(config.system.distMPc)
        if minr == 'auto':
            minr = 0.0
        if maxr == 'auto':
            maxr = 10 ** config.settings.orblib_settings['logrmax'] * arc_kpc
        self.minr = minr
        self.maxr = maxr
        self.minr_arcsec = self.minr / arc_kpc
        self.maxr_arcsec = self.maxr / arc_kpc
        self.nr = nr
        self.nl = nl
        self.logger.info(f'Coloring initialized with minr={self.minr} kpc, '
                         f'maxr={self.maxr} kpc, nr={self.nr}, nl={self.nl}.')

    def bin_phase_space(self,
                        model=None,
                        vor_weight=0.01,
                        vor_ignore_zeros=False,
                        make_diagnostic_plots=False,
                        extra_diagnostic_output=False,
                        cvt=False,
                        wvt=False,
                        use_cache=True):
        """
        Perform Voronoi binning of orbits in the radius-circularity phase
        space. The goal is to group the "original" n_orbits orbit bundles into
        fewer n_bundle "Voronoi" bundles with each of these Voronoi bundles
        accounting for a weight of at least ``vor_weight``. The resulting
        orbit bundle mapping is returned and written to the model directory.
        This method uses the ``vorbin`` package by Michele Cappellari
        (M. Cappellari, Y. Copin, 2003, MNRAS 342, 345).

        Parameters
        ----------
        model : a ``dynamite.model.Model`` object, optional
            the model used for the binning. If ``None``, choose the minimum
            :math:`\\chi^2` model (the configuration setting ``which_chi2``
            determines the :math:`\\chi^2` type). The default is ``None``.
        vor_weight : float, optional
            the target total orbital weight in each Voronoi bin.
            The default is 0.01.
        vor_ignore_zeros : bool, optional
            If ``True``, then radius-circularity bins that have zero total
            weight will be ignored in the binning process. NOTE: this
            feature is still EXPERIMENTAL due to insufficient testing, may
            still be buggy, and is therefore NOT RECOMMENDED.
            The default is ``False``.
        make_diagnostic_plots : bool, optional
            If ``True``, both vorbin and DYNAMITE will produce diagnostic
            plots visualizing the binning result. The default is ``False``.
        extra_diagnostic_output : bool, optional
            If ``True``, vorbin will print details on the binning to stdout.
            The default is ``False``.
        cvt : bool, optional
            If ``True``, the Voronoi binning will be performed using the
            CVT (Centroidal Voronoi Tessellation) algorithm. For details,
            refer to the ``vorbin`` documentation. The default is ``False``.
        wvt : bool, optional
            If ``True``, the Voronoi binning will be performed using the
            WVT (Weighted Voronoi Tessellation) algorithm (Diehl & Statler
            2006, MNRAS, 368, 497). For details, refer to the ``vorbin``
            documentation. The default is ``False``.
        use_cache : bool, optional
            If ``True``, the method will use cached results if available.
            If the Voronoi binning has already been performed, the results
            will be loaded from the cache in the model directory instead of
            recalculating them. The default is ``True``.

        Returns
        -------
        (vor_bundle_mapping, phase_space_binning) : tuple, where
        vor_bundle_mapping : np.array of shape (n_bundle, n_orbits)
            Mapping between the "original" orbit bundles and the Voronoi
            orbit bundles: vor_bundle_mapping(i_bundle, i_orbit) is the
            fraction of i_orbit assigned to i_bundle, multiplied by i_orbit's
            weight.
        phase_space_binning : dict
            'in': np.array of shape (3, self.nr * self.nl), the binning input:
            bin r, bin lambda_z, bin total weight
            'out': np.array of shape (3, n_bundle), the Voronoi binning output:
            weighted Voronoi bin centroid coordinates r_bar, lambda_bar
            and Voronoi bin total weights
            'map': np.array of shape (self.nr * self.nl,) the phase space
            mapping: Voronoi bin numbers for each input bin

        """
        # save parameters relevant for caching
        params = dict(locals())
        params.pop('self', None)  # remove self from the parameters
        params.pop('model', None)  # remove model from the parameters
        params['code_version'] = '0.1'  # version of this method
        params['nr'] = self.nr
        params['nl'] = self.nl
        params['minr'] = self.minr
        params['maxr'] = self.maxr
        compare_keys = [k for k in params]
        params['f_name'] = f_name = None
        if model is None:
            best_model_idx = self.config.all_models.get_best_n_models_idx(1)[0]
            model = self.config.all_models.get_model_from_row(best_model_idx)
        f_metadata = model.directory + 'voronoi_orbit_bundles.yaml'
        bundle_metadata = []

        if use_cache:
            try:
                with open(f_metadata, 'r') as f:
                    bundle_metadata = yaml.safe_load(f)
                self.logger.info('Voronoi binning metadata read from '
                                f'{f_metadata}.')
            except FileNotFoundError:
                self.logger.info('No Voronoi binning metadata found.')
            for dat in bundle_metadata:
                if set(dat.keys()) == set(params.keys()) and \
                   all(dat[k] == params[k] for k in compare_keys):
                    f_name = dat['f_name']
                    self.logger.info('Reading Voronoi orbit bundle data from '
                        f'existing file {model.directory}{f_name}.')
                    self.logger.debug(f'Existing metadata: {dat}.')
                    break
            if f_name is not None:  # existing Voronoi binning data found
                data = np.load(model.directory + f_name)
                vor_bundle_mapping = data['vor_bundle_mapping']
                phase_space_binning = \
                    dict((iom, data[iom]) for iom in ('in', 'out', 'map'))
                return vor_bundle_mapping, phase_space_binning  ############
            else:
                self.logger.info('Performing phase space Voronoi binning.')

        orblib = model.get_orblib()
        _ = model.get_weights(orblib)
        weights = model.weights

        orblib.get_projection_tensor(minr=self.minr,
                                     maxr=self.maxr,
                                     nr=self.nr,
                                     nl=self.nl,
                                     force_lambda_z=True,
                                     r_scale='linear')
        # pick entry [2] = fraction of orbits in each r,l bin (ALL orbit types)
        # indices r, l, b = radius, lambda_z, original_orbit_bundle
        projection_tensor = orblib.projection_tensor[2]

        # get the orbit distribution in the r, lambda space
        # = total weight in the r, l bins
        orbit_distribution = projection_tensor.dot(weights)

        # build the input for vorbin
        dr = (self.maxr - self.minr) / self.nr
        input_bins_r = np.linspace(self.minr + dr / 2, self.maxr - dr / 2,
                                   num=self.nr)
        dl = 2 / self.nl
        input_bins_l = np.linspace(-1 + dl / 2, 1 - dl / 2, num=self.nl)
        input_grid = np.meshgrid(input_bins_r,input_bins_l)
        vor_in = np.vstack((np.array([m.ravel() for m in input_grid]),
                            orbit_distribution.T.ravel()))
        if vor_ignore_zeros:
            nonzero_weight_bins = vor_in[2,:] > 0
            vor_in = vor_in[:,nonzero_weight_bins]

        def _sum_weights(index, signal, noise):
            """
            Function to calculate the "S/N" (here: sum of orbital weights) of
            a Voronoi bin comprising radius, lambda_z bins given in index

            Parameters
            ----------
            index : int vector
                integer vector of length N containing indices of the r, l bins
                for which the combined orbital weight has to be returned.
                The indices refer to elements of the vectors signal and noise..
            signal : float vector
                vector of length M>N with the orbital weights in all r,l bins
            noise : float vector
                vector of length M>N with the noise of all r,l bins, ignored
                (provided for compatibility with vorbin).

            Returns
            -------
            float
                sum of all orbital weights in r,l bins given in index

            """
            return np.sum(signal[index])

        binning, _, _, r_bar, l_bar, bin_weights, _, _ = \
            voronoi_2d_binning(*vor_in,
                               noise=np.ones_like(vor_in[0]),
                               target_sn=vor_weight,
                               plot=make_diagnostic_plots,
                               sn_func=_sum_weights,
                               quiet=not extra_diagnostic_output,
                               cvt=cvt,
                               wvt=wvt)
        phase_space_binning = {'in':vor_in,
                               'out':np.vstack((r_bar, l_bar, bin_weights)),
                               'map':binning}  # Vor bin numbers of r, l bins
        n_orbits = projection_tensor.shape[-1]  # "original" orbit bundles
        n_bundle = phase_space_binning['out'].shape[-1]  # Voronoi bundles
        self.logger.info(f'{n_orbits} original orbit bundles aggregated into '
                         f'{n_bundle} Voronoi bundles.')
        if make_diagnostic_plots:
            if vor_ignore_zeros:
                vor_bin_map = iter(phase_space_binning['map'])
                vor_map_dense = [next(vor_bin_map) if nz else np.nan
                                 for nz in nonzero_weight_bins]
            else:
                vor_map_dense = phase_space_binning['map']
            vor_map_dense = np.array(vor_map_dense).reshape(self.nl, self.nr)
            plt.figure()
            plt.imshow(orbit_distribution.T,
                       interpolation='spline16',
                       origin='lower',
                       extent=(self.minr, self.maxr, -1, 1),
                       aspect='auto',
                       cmap='Greys',
                       alpha=0.9)
            plt.pcolormesh(*input_grid,
                           vor_map_dense,
                           shading='nearest',
                           cmap='tab20',
                           alpha=0.3)
            plt.xlabel('$r$ [kpc]')
            plt.ylabel('Circularity $\\lambda_z$')
            plt.colorbar(label='Voronoi bin id')
            plt.savefig(f'{self.config.settings.io_settings["plot_directory"]}'
                        'voronoi_orbit_bundles.png')

        # map "original" orbit bundles to Voronoi bundles
        # calculate weighted orbit fractions in the input bins: "flatten"
        # r and l indices of the projection tensor making r the fastest moving
        # index, like in phase_space_binning['map']==binning and multiply
        # element-wise with the orbit weights; orb_frac.shape=(nr*nl, n_orbits)
        orb_frac = projection_tensor.todense().reshape((self.nr * self.nl, -1),
                                                       order='F') * weights
        if vor_ignore_zeros:  # eliminate bins with zero total weight
            orb_frac = orb_frac[nonzero_weight_bins]
        # each Voronoi bin number represents a Voronoi (orbit) bundle
        # for each "original" (orbit) bundle, calculate its fractions that go
        # into each Voronoi bundle
        vor_bundle_mapping = np.zeros((n_bundle, n_orbits))
        for i in range(orb_frac.shape[0]):  # collect weighted orbit fractions
            vor_bundle_mapping[binning[i]] += orb_frac[i]
        # vor_bundle_mapping[binning[:]] += orb_frac[:]

        # save the Voronoi binning results to the model directory
        # Voronoi binning metadata already exists? -> use that filename
        if bundle_metadata == []:  # otherwise, it has been read already
            try:
                with open(f_metadata, 'r') as f:
                    bundle_metadata = yaml.safe_load(f)
                self.logger.info('Voronoi binning metadata read from '
                                f'{f_metadata}.')
            except FileNotFoundError:
                self.logger.info('No Voronoi binning metadata exists.')
        for dat in bundle_metadata:
            if set(dat.keys()) == set(params.keys()) and \
               all(dat[k] == params[k] for k in compare_keys):
                self.logger.info('Voronoi binning metadata exists, will '
                                 'overwrite bundle data.')
                self.logger.debug(f'Existing metadata: {dat}.')
                f_name = dat['f_name']
                break
        else:
            # no existing metadata found, create anew
            idx = len(bundle_metadata)
            while True:
                f_name = f'voronoi_orbit_bundles_{idx:03d}.npz'
                if not os.path.isfile(model.directory + f_name):
                    break
                idx += 1
            params['f_name'] = f_name
            self.logger.debug(f'Creating metadata: {params}.')
            bundle_metadata.append(params)
            with open(f_metadata, 'w') as f:
                yaml.dump(bundle_metadata, f)
            self.logger.info('Voronoi binning metadata saved in '
                             f'{f_metadata}.')
        np.savez(model.directory + f_name,
                 vor_bundle_mapping=vor_bundle_mapping,
                 **phase_space_binning)
        self.logger.info('Voronoi orbit bundle mapping saved in '
                         f'{model.directory + f_name}.')

        return vor_bundle_mapping, phase_space_binning

    def fit_bayesian(self,
                     prior_dist,
                     prior_par,
                     flux_data_norm,
                     obs_data,
                     sample):
        """Fit orbit bundle color maps to the observed data using Bayesian
        inference.

        Employs Probabilistic Programming with PyMC to fit the observed
        quantity (e.g., age or metallicity) based on the normalized flux data
        for the spatial bins. The model uses a truncated normal or lognormal
        distribution for the prior of the observed quantity and a Student's t
        distribution with fixed sigma (Half-Cauchy distributed with beta=5)
        and nu (Exponential distributed with parameter 1/30) parameters for
        the likelihood of the observed data.
        This method uses the Markov chain Monte Carlo (MCMC) sampling
        algorithm NUTS (No-U-Turn Sampler), initialized with the ADVI
        (Automatic Differentiation Variational Inference) method.
        The code has been strongly inspired by a similar implementation
        provided by Ling Zhu as used in e.g., Zhu et al. 2022, A&A, 664, A115.

        Parameters
        ----------
        prior_dist : str
            Distribution for the prior of the observed quantity (e.g., age or
            metallicity). Currently implemented:
            'normal': truncated normal distribution
            'lognormal': truncated lognormal distribution
            Other distributions will raise a NotImplementedError.
        prior_par : dict
            Parameters for the prior distribution.
            For 'normal', it should contain 'mu', 'sigma', 'lower', and
            'upper' keys, where 'mu' and 'sigma' are the mean and standard
            deviation of the normal distribution, and 'lower' and 'upper' are
            the truncation limits.
            For 'lognormal', it should contain 'mu', 'sigma', 'lower', and
            'upper' keys, where 'mu' and 'sigma' are the parameters of the
            lognormal distribution, and 'lower' and 'upper' are the truncation
            limits.
            The 'mu' values should be a 1d numpy array of length equal to the
            number of orbit bundles fitted to the observed data. The rest of
            the parameters should be scalars.
        flux_data_norm : np.array of shape (n_spatial_bins, n_bundle)
            Normalized flux data for the spatial bins, where each column
            corresponds to an orbit bundle and each row corresponds to a
            spatial bin. This data is used to compute the expected value of
            the observed quantity based on the fitted model.
            The normalization is done such that the sum of fluxes in each
            spatial bin is equal to 1.
            This is typically the result of orbit bundle maps calculated from a
            Voronoi binning process of the phase space.
        obs_data : np.array of shape (n_spatial_bins,)
            Observed data for the spatial bins, which is the quantity to be
            fitted. This data is used as the observed values in the Bayesian
            inference process. Typical examples are age or metallicity.
        sample : dict
            Parameters for the sampling process, including:
            - 'n_draws': int, number of draws for the MCMC sampling.
            - 'n_tune': int, number of tuning steps for the MCMC sampling.
            - 'advi_init': int, number of initialization steps for ADVI.

        Returns
        -------
        (model, trace) : tuple, where
            model : pymc.Model
                The PyMC model object containing the Bayesian model
                specification.
            trace : arviz.data.inference_data.InferenceData
                The trace of the MCMC sampling, containing the posterior
                distributions of the model parameters.

        Raises
        ------
        NotImplementedError
            If the distribution specified in `prior_dist` is not implemented.
        """
        self.logger.info('Fitting Bayesian model to the observed data.')
        with pm.Model() as model:
            if prior_dist == 'normal':  # truncated normal distribution
                qty = pm.TruncatedNormal(name='qty',
                                         mu=prior_par['mu'],
                                         sigma=prior_par['sigma'],
                                         lower=prior_par['lower'],
                                         upper=prior_par['upper'],
                                         shape=prior_par['mu'].size)
            elif prior_dist == 'lognormal':  # truncated lognormal distribution
                LogNormalDist = pm.LogNormal.dist(mu=prior_par['mu'],
                                                  sigma=prior_par['sigma'],
                                                  shape=prior_par['mu'].size)
                qty = pm.Truncated(name='qty',
                                   dist=LogNormalDist,
                                   lower=prior_par['lower'],
                                   upper=prior_par['upper'],
                                   shape=prior_par['mu'].size)
            else:
                txt = f'Prior distribution {prior_dist} not implemented.'
                self.logger.error(txt)
                raise NotImplementedError(txt)
            student_t_sigma = pm.HalfCauchy('student_t_sigma', beta=5)
            student_t_nu = pm.Exponential('student_t_nu', 1/30)
            student_t_mu = pm.math.dot(flux_data_norm, qty)
            obs = pm.StudentT('obs',
                              mu=student_t_mu,
                              sigma=student_t_sigma,
                              nu=student_t_nu,
                              observed=obs_data)
            trace = pm.sample(draws=sample['n_draws'],
                              tune=sample['n_tune'],
                              step=None,
                              init='advi',
                              n_init=sample['advi_init'])
        self.logger.info('Bayesian model fitting completed. Shape of '
                         'posterior: (n_chain, n_draw, n_dim_0) = '
                         f'{trace.posterior["qty"].shape}.')
        return model, trace

    def orbit_bundle_plot(self,
                          orbit_distribution=None,
                          model=None,
                          phase_space_mapping=None,
                          figtype='.png',
                          dpi=100):
        """Create a plot of the orbit bundle distribution in the
        radius-circularity phase space.

        Parameters
        ----------
        orbit_distribution : np.array of shape (self.nr, self.nl), optional
            Probability density of stellar orbits in the (r, lambda_z) bins.
            If None, the orbit distribution is calculated from the model
            parameter. The default is None.
        model : a ``dynamite.model.Model`` object, optional
            The model used for the binning if orbit_distribution is not
            specified or None. This parameter is ignored if orbit_distribution
            is given. If model is None, choose the minimum
            :math:`\\chi^2` model (the configuration setting ``which_chi2``
            determines the :math:`\\chi^2` type). The default is None.
        phase_space_mapping : np.array of shape (self.nr * self.nl,)
            The phase space mapping: Voronoi bin numbers for each input bin.
            Can directly use the bin_phase_space() output
            phase_space_binning['map'].
        figtype : str, optional
            Determines the file format of the saved figure, by default '.png'.
        dpi : int, optional
            The resolution of saved figure, by default 100.

        Returns
        -------
        matplotlib.figure.Figure
            The created figure object.

        Raises
        ------
        ValueError
            If phase_space_binning is None.
        """
        self.logger.info('Creating orbit bundle plot.')
        if orbit_distribution is None:
            self.logger.info('No orbit distribution provided, '
                             'calculating from model parameter.')
            if model is None:
                best_model_idx = \
                    self.config.all_models.get_best_n_models_idx(1)[0]
                model = \
                    self.config.all_models.get_model_from_row(best_model_idx)
            orblib = model.get_orblib()
            orblib.get_projection_tensor(minr=self.minr,
                                         maxr=self.maxr,
                                         r_scale='linear',
                                         nr=self.nr,
                                         nl=self.nl,
                                         force_lambda_z=True)
            _ = model.get_weights(orblib)
            orbit_distribution = orblib.projection_tensor[2].dot(model.weights)
        if phase_space_mapping is None:
            txt = 'No phase space mapping provided, cannot create ' \
                  'orbit bundle plot.'
            self.logger.error(txt)
            raise ValueError(txt)
        vor_map_dense = np.array(phase_space_mapping).reshape(self.nl, self.nr)
        fig,ax = plt.subplots(figsize=(7, 5))
        plt.imshow(orbit_distribution.T,
                interpolation='spline16',
                origin='lower',
                extent=(self.minr_arcsec, self.maxr_arcsec, -1, 1),
                aspect='auto',
                cmap='Greys',
                alpha=0.9)
        plt.colorbar(label='Orbit probability density')
        plt.xlabel(r'$r$ [arcsec]')
        plt.ylabel(r'Circularity $\lambda_z$')
        if self.maxr >= 1:
            r_unit = 'kpc'
            factor = constants.ARC_KPC(self.config.system.distMPc)
        else:
            r_unit = 'pc'
            factor = constants.ARC_KPC(self.config.system.distMPc) * 1e3
        secax = ax.secondary_xaxis('top', functions=(lambda x: x * factor,
                                                     lambda x: x / factor))
        secax.set_xlabel(f'$r$ [{r_unit}]')
        dr = (self.maxr_arcsec - self.minr_arcsec) / self.nr
        dl = 2 / self.nl
        for i_l, lam in enumerate(np.linspace(-1, 1 - dl, num=self.nl)):
            for i_r, r in enumerate(np.linspace(self.minr_arcsec,
                                                self.maxr_arcsec - dr,
                                                num=self.nr)):
                if i_r < self.nr - 1 \
                   and vor_map_dense[i_l, i_r] != vor_map_dense[i_l, i_r + 1]:
                    plt.plot([r + dr, r + dr],
                             [lam, lam + dl],
                             color='tab:red',
                             linewidth=0.5)
                if i_l < self.nl - 1 \
                   and vor_map_dense[i_l, i_r] != vor_map_dense[i_l + 1, i_r]:
                    plt.plot([r, r + dr],
                             [lam + dl, lam + dl],
                             color='tab:red',
                             linewidth=0.5)
        # Save the figure
        figname = self.config.settings.io_settings['plot_directory'] + \
            'orbit_bundle_plot' + figtype
        fig.savefig(figname, dpi=dpi)
        self.logger.info(f'Orbit bundle plot saved in {figname}.')
        return fig

    def get_color_orbital_decomp(self, models, vor_bundle_mappings, colors):
        """Calculate orbital decomposition of the coloring data

        This method calculates the orbital decomposition of the coloring data
        for a set of DYNAMITE models, given the Voronoi bundle mappings and
        colors. For each model, it computes the total orbit weight in each
        (r, lambda_z) phase space bin and the color distribution in those bins.
        The method assumes that all DYNAMITE models have the same number of
        colors and that the Voronoi orbital bundle mappings and colors are
        provided in the same order as the models. It returns an array of shape
        (n_colors + 1, nr, nl, n_models) with the total orbit weight and
        n_colors color distributions in each of the nr * nl (r, lambda_z) bins
        for each of the n_models DYNAMITE models.

        Parameters
        ----------
        models : list of ``dynamite.model.Model`` objects
            List of DYNAMITE models for which the orbital decomposition is
            calculated.
        vor_bundle_mappings : list of np.arrays
            List of Voronoi orbit bundle mappings for each model. Each mapping
            is a 2D numpy array where the first dimension corresponds to the
            Voronoi orbit bundles and the second dimension corresponds to the
            original orbit bundles: its shape is (n_bundle, n_orbits).
        colors : list of lists of np.arrays
            Coloring data in a list. There is one entry per model. Those
            entries are again lists; their length is n_colors. They consist
            of 1D numpy arrays of shape (n_bundle,), each of which holds data
            of one color in the Voronoi orbit bundles for that model. The
            order of the colors must be the same for each model. Example for
            2 models and 2 colors 'age' and 'metallicity':
            ``colors = [[age_model1, met_model1], [age_model2, met_model2]]``

        Returns
        -------
        orbit_weight_and_color_distribution :
            np.array of shape (n_colors + 1, nr, nl, n_models)
            Array containing both the total orbit weight and the color
            distribution in each (r, lambda_z) phase space bin for each model.
            The first slice along the first dimension contains the
            orbit weight, and the subsequent slices contain the color
            distributions. In (r, lambda_z) bins without any orbits, the
            color_distribution will contain np.nan values.

        Raises
        ------
        ValueError
            If the number of colors is not the same for all models, or if the
            number of models in the arguments do not match. The error message
            will indicate the issue.
        """
        if len(set(map(len,colors))) != 1:
            txt = 'All models need to have an equal number of colors.'
            self.logger.error(txt)
            raise ValueError(txt)
        if len(set(map(len,[models, vor_bundle_mappings, colors]))) != 1:
            txt = 'models, vor_bundle_mappings, and colors must be of equal ' \
                  'length.'
            self.logger.error(txt)
            raise ValueError(txt)
        n_colors = len(colors[0])
        n_models = len(models)
        self.logger.info(f'Calculating color orbital decomposition for '
                         f'{n_models} models with {n_colors} colors each.')
        # orbit_weight collects the total orbit weight in each (r, l) bin
        orbit_weight = np.zeros((self.nr, self.nl, n_models))
        # color_distribution will collect the colors of each phase space bin
        # (bins without any orbits get np.nan):
        color_distribution = np.full((n_colors, self.nr, self.nl, n_models),
                                     np.nan)
        for i_model, model in enumerate(models):  # for each model...
            # PART 1: get the total orbit weights in the (r, l) phase space
            orblib = model.get_orblib()
            orblib.get_projection_tensor(minr=self.minr,
                                         maxr=self.maxr,
                                         r_scale='linear',
                                         nr=self.nr,
                                         nl=self.nl,
                                         force_lambda_z=True)
            _ = model.get_weights(orblib)
            # Now, use projection_tensor[2] = fraction of orbits in each r, l
            # bin (based on ALL orbit types); its shape = (nr, nl, n_orbits).
            # Calculate total orbit weight in each (nr, nl) bin:
            orbit_weight[..., i_model] = \
                orblib.projection_tensor[2].dot(model.weights)
            # PART 2: get all the colors in the phase space
            # Identify (r, l) bins with nonzero total orbit weight:
            valid_rl = [(r, lam) for r in range(self.nr)
                                 for lam in range(self.nl)
                                 if orbit_weight[r, lam, i_model]]
            # get the colors of each original orbit bundle
            vor_bundle_mapping = vor_bundle_mappings[i_model]
            _, n_orbits = vor_bundle_mapping.shape
            orbit_color = np.zeros((n_colors, n_orbits))
            for i_orbit in range(n_orbits):
                # Need to catch sum(weights)=0 for weightewd average to work:
                if np.any(vor_bundle_mapping[:,i_orbit]):
                    for i_color, color in enumerate(colors[i_model]):
                        # weighted average in case dithering!=0
                        orbit_color[i_color, i_orbit] = \
                            np.average(color,
                                       weights=vor_bundle_mapping[:,i_orbit])
            # Add weights of original orbit bundles that contribute to the
            # phase space bins:
            color_distribution_weights = \
                orblib.projection_tensor[2] * model.weights  # element-wise
            for i_color in range(n_colors):
                for (i_r, i_l) in valid_rl:  # 'invalid' (nr, nl) stay np.nan
                    color_distribution[i_color, i_r, i_l, i_model] = \
                        np.average(orbit_color[i_color], weights=\
                            color_distribution_weights[i_r, i_l].todense())
            self.logger.info(f'Color distribution for model {i_model} done.')
        return np.concatenate((orbit_weight[np.newaxis], color_distribution),
                              axis=0)

    def coloring_decomp_plot(self,
                             orbit_data,
                             plot_labels=None,
                             colorbar_scale='linear',
                             rcut_kpc=3.5,
                             lz_disk=0.8,
                             lz_warm=0.5,
                             figtype='.png',
                             dpi=100):
        """Coloring decomposition plot, averaging the data of multiple models

        Create and save an orbital decomposition plot in the (r, lambda_z)
        phase space, consisting of multiple panels: the first panel shows
        the orbit probability density distribution, the subsequent panels
        show the color distributions. Dashed lines indicate the orbit-based
        division into four components: disk, warm, bulge, and hot inner
        stellar halo (can be switched off). The plot averages the data of
        multiple DYNAMITE models.

        Parameters
        ----------
        orbit_data : np.array of shape (n_colors + 1, nr, nl, n_models)
            Array containing both the total orbit weight and the color
            distribution in each (r, lambda_z) phase space bin for each of the
            n_models DYNAMITE models. The first slice along the first dimension
            contains the orbit weight, and the subsequent slices contain the
            color distributions. Can directly use the output of
            ``get_color_orbital_decomp()``.
        plot_labels : list of str or ``None``, optional
            Labels for the individual plots. Must have length of n_colors + 1.
            If ``None``, the default labels will be used:
            'Orbit PDF', 'Color 1', 'Color 2', etc. The default is ``None``.
        colorbar_scale : str or list of str, optional
            Scale of the colorbar, either 'linear' or 'log'. If a string, it
            applies to all plots. If a list of strings of length n_colors + 1,
            it sets the scale of each plot individually. The default is
            'linear'.
        Cuts for stellar components:
            - disk (lambda_z > lz_disk)
            - bulge (lambda_z < lz_disk, r < rcut_kpc)
            - warm (lz_warm < lambda_z < lz_disk, rcut_kpc < r)
            - hot inner stellar halo (lambda_z < lz_warm, rcut_kpc < r)
            If one of the cuts is ``None``, no lines will be drawn.
        rcut_kpc : float or ``None``, optional
            by default 3.5
        lz_disk : float or ``None``, optional
            by default 0.8
        lz_warm : float or ``None``, optional
            by default 0.5
        figtype : str, optional
            Determines the file format of the saved figure, by default '.png'.
        dpi : int, optional
            The resolution of saved figure, by default 100.

        Returns
        -------
        matplotlib.figure.Figure
            The created figure object.

        Raises
        ------
        ValueError
            If colorbar_scale is invalid. The error message will indicate the
            issue.
        """
        n_models = orbit_data.shape[-1]
        self.logger.info('Creating coloring decomposition plot '
                         f'averaging data of {n_models} models.')
        arctkpc = constants.ARC_KPC(self.config.system.distMPc)
        if rcut_kpc is None or lz_disk is None or lz_warm is None:
            self.logger.info('One of the cuts for stellar components is None, '
                             'no lines will be drawn in the plot.')
            rcut_arc = None
        else:
            self.logger.info(f'Using rcut={rcut_kpc} kpc, '
                             f'lz_disk={lz_disk}, lz_warm={lz_warm} for '
                             'the component cuts in the plot.')
            rcut_arc = rcut_kpc / arctkpc
        minlz, maxlz = -1, 1
        minr_arc = self.minr / arctkpc
        maxr_arc = self.maxr / arctkpc
        all_rl = [(r, lam) for r in range(self.nr) for lam in range(self.nl)]
        n_plots = len(orbit_data)
        fig = plt.figure(figsize=(16 / 3 * n_plots, 6))
        if plot_labels is None or not isinstance(plot_labels, list) or \
           len(plot_labels) != n_plots:
            plot_labels = ['Orbit PDF'] + \
                           [f'Color {i + 1}' for i in range(n_plots - 1)]
            self.logger.warning('plot_labels should be a list of labels '
                                f'for each plot, using default {plot_labels}.')
        if colorbar_scale == 'linear' or colorbar_scale == 'log':
            self.logger.info(f'Using {colorbar_scale} colorbar scale.')
            colorbar_scale = [colorbar_scale] * n_plots
        elif isinstance(colorbar_scale, list) and \
            len(colorbar_scale) == n_plots and \
            all([s == 'linear' or s == 'log' for s in colorbar_scale]):
            self.logger.info('Using individual colorbar scales for each plot.')
        else:
            txt = 'colorbar_scale needs to be either "linear" or "log", ' \
                  f'or a list of those with length {n_plots}.'
            self.logger.error(txt)
            raise ValueError(txt)
        for i_plot, plot_data in enumerate(orbit_data):
            label = plot_labels[i_plot]
            # The first entry in orbit_data is the orbits' weight, convert
            # to probability density distribution of the orbits in the (r, l)
            # phase space:
            if i_plot == 0:
                weight = plot_data
                values = np.sum(weight, axis=2)
                values = values / np.sum(values)
                values[values == 0.] = np.nan  # avoid zero values for colorbar
            else:  # Average over all models
                # (r, l) bins without orbits get a color value of np.nan
                values = np.full((self.nr, self.nl), np.nan)
                for i_rl in [i_rl for i_rl in all_rl if np.any(weight[i_rl])]:
                    values[i_rl] = np.average(plot_data[i_rl],
                                              axis=0,
                                              weights=weight[i_rl])
            ax = fig.add_subplot(1, n_plots, i_plot + 1)
            if i_plot == 0:
                if colorbar_scale[i_plot] == 'linear':
                    vmin = 0  # orbit_pdf_min=0
                else:
                    vmin = np.nanmin(values)
                cmap = cmasher.get_sub_cmap('gist_heat_r', 0.05, 1.)
            else:
                vmin = np.nanmin(values)
                cmap = cmasher.get_sub_cmap('cubehelix', 0, 0.9)
            im = ax.imshow(values.T,
                           cmap=cmap,
                           norm=colorbar_scale[i_plot],
                           vmin=vmin,
                           vmax=np.nanmax(values),
                           aspect='auto',
                           interpolation='none',
                           origin='lower',
                           extent=(minr_arc, maxr_arc, minlz, maxlz))
            fig.colorbar(im,
                         orientation='horizontal',
                         location='top',
                         pad=0.15,
                         label=label)
            if rcut_arc is not None:
                plt.plot((rcut_arc, rcut_arc), (-1, lz_disk), 'k--')
                plt.plot((minr_arc, maxr_arc), (lz_disk, lz_disk), 'k--')
                plt.plot((rcut_arc, maxr_arc), (lz_warm, lz_warm), 'k--')
            # ax.set_yticks([-1, 0, 1])
            # ax.set_xticks([0, 50, 100])
            ax.set_xlabel(r'$r$ [arcsec]')
            if self.maxr >= 1:
                r_unit = 'kpc'
                factor = arctkpc
            else:
                r_unit = 'pc'
                factor = arctkpc * 1e3
            secax = ax.secondary_xaxis('top', functions=(lambda x: x * factor,
                                                         lambda x: x / factor))
            secax.set_xlabel(f'$r$ [{r_unit}]')
            if i_plot == 0:  # y-axis label and components only for first plot
                ax.set_ylabel(r'Circularity $\lambda_z$')
                if rcut_arc is not None:
                    ax.text(rcut_arc * 1.1,
                            (lz_disk + 1) / 2,
                            'Disk',
                            verticalalignment='center')
                    ax.text(rcut_arc * 1.1,
                            (lz_warm + lz_disk) / 2,
                            'Warm',
                            verticalalignment='center')
                    ax.text(rcut_arc * 1.1,
                            lz_warm / 2,
                            'Hot inner stellar halo',
                            verticalalignment='center')
                    ax.text(rcut_arc * 0.2,
                            -0.25,
                            'Bulge',
                            verticalalignment='center')
        # Save the figure
        figname = self.config.settings.io_settings['plot_directory'] + \
                  f'coloring_decomp_plot_{n_models:02d}' + figtype
        fig.savefig(figname, dpi=dpi)
        self.logger.info(f'Coloring decomposition plot saved in {figname}.')
        return fig

    def circularity_color_plot(self,
                               weights,
                               color_data,
                               *,
                               c_label='Color',
                               c_scale='linear',
                               p_scale='linear',
                               n_color_bins=14,
                               interpolation='none',
                               disk_fraction=True,
                               lz_disk=0.8,
                               figtype='.png',
                               dpi=100):
        """Plots orbit probability distribution in the (color, circularity)
        plane and the fraction of disk orbits as a function of the color

        Create and save a plot showing the orbit probability distribution in
        the (color, circularity) plane, averaged over multiple DYNAMITE models.
        If disk_fraction is True, the cold orbit fraction as a function of the
        color is plotted in a separate panel above the main plot. The cold
        orbit fraction is defined as the fraction of orbits with circularity
        $\lambda_z$ greater than a specified threshold lz_disk (default is 0.8)
        within each color bin. A dashed line indicates the disk fraction
        threshold and a vertical dashed line indicates the color value at which
        the cold orbit fraction crosses 50%. The plot is saved in the specified
        file format and resolution.

        Parameters
        ----------
        weights : np.array of shape (nr, nl, n_models)
            Array with the total orbit weight in each (r, lambda_z) phase space
            bin for each DYNAMITE model. Can directly use the first slice of
            the output of ``get_color_orbital_decomp()``.
        color_data : np.array of shape (nr, nl, n_models)
            Array with the color distribution in each (r, lambda_z) phase
            space bin for each DYNAMITE model. Can directly use a
            slice of the output of ``get_color_orbital_decomp()``.
        c_label : str, optional
            Label for the color axis, by default 'Color'.
        c_scale : str, optional
            Scale of the color axis, either 'linear' or 'log', by default
            'linear'.
        p_scale : str, optional
            Scale of the probability density colorbar, either 'linear' or
            'log', by default 'linear'.
        n_color_bins : int, optional
            Number of color bins; determines the resolution along the color
            axis of both the orbit distribution and cold orbit fraction plots.
            The bins will be evenly and linearly spaced between the minimum and
            maximum age values in the orbit_data. The default is 14.
        interpolation : str, optional
            Interpolation method for the imshow plot of the orbit distribution.
            See
            https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
            for available options. The default is 'none'.
        disk_fraction : bool, optional
            Whether to plot the cold orbit fraction as a function of color in
            a separate panel above the main plot. If False, the plot will have
            only one pane showing the orbit probability distribution.
            The default is True.
        lz_disk : float, optional
            Circularity threshold for the disk fraction, above which orbits are
            considered as part of the disk. This is used to calculate the cold
            orbit fraction. The default value is 0.8.
        figtype : str, optional
            Determines the file format of the saved figure, by default '.png'.
        dpi : int, optional
            The resolution of saved figure, by default 100.

        Returns
        -------
        (figure, disk_50) : tuple, where
            figure : matplotlib.figure.Figure
                The created figure object.
            disk_50 : float or np.nan
                The color value at which the cold orbit fraction crosses 50%
                (disk fraction threshold). Calculated via linear interpolation
                between the two color bins that cross the threshold. Set to
                np.nan if the cold orbit fraction does not cross the threshold
                in the given color range or if disk_fraction is False.

        Raises
        ------
        ValueError
            If the shape of `orbit_data` is not compatible with the expected
            format, i.e., it should have at least two slices along the first
            dimension (one for orbit weight and one for age).
        """
        if c_scale != 'linear' and c_scale != 'log':
            txt = 'c_scale needs to be either "linear" or "log", ' \
                  f'but is {c_scale}.'
            self.logger.error(txt)
            raise ValueError(txt)
        if p_scale != 'linear' and p_scale != 'log':
            txt = 'p_scale needs to be either "linear" or "log", ' \
                  f'but is {p_scale}.'
            self.logger.error(txt)
            raise ValueError(txt)
        if weights.shape != color_data.shape:
            txt = 'Weights and color_data must have the same shape, but ' \
                  f'they are {weights.shape} and {color_data.shape}, resp.'
            self.logger.error(txt)
            raise ValueError(txt)
        n_models = weights.shape[-1]
        weights = np.ravel(weights)  # shape = (nr * nl * n_models)
        color_data = np.ravel(color_data)  # shape = (nr * nl * n_models)
        self.logger.info('Creating circularity vs. color '
                        f'plot averaging data of {n_models} models.')
        # Calculate circularity vs. age values
        min_c, max_c = np.nanmin(color_data), np.nanmax(color_data)
        lambda_z = np.zeros((self.nr, self.nl, n_models))
        for i_l in range(self.nl):
            lambda_z[:, i_l, :] = -1 + 2 / (self.nl - 1) * i_l
        lambda_z = np.ravel(lambda_z)
        values, color_edges, _ = \
            np.histogram2d(color_data,
                           lambda_z,
                           weights=weights,
                           bins=np.array([n_color_bins, self.nl]),
                           range=[[min_c, max_c], [-1, 1]])
        values = values / np.sum(values)
        if disk_fraction:  # Calculate cold orbit (disk) fraction
            c_cold = np.linspace(np.mean(color_edges[0:2]),
                                 np.mean(color_edges[-2:]),
                                 n_color_bins)
            f_cold = np.zeros_like(c_cold)
            # Make sure to include the full last bin
            color_edges[-1] *= 1 + np.finfo(np.float64).eps * 10
            for i in range(n_color_bins):
                f_cold[i] = np.sum(weights[(color_data >= color_edges[i]) \
                                           & (color_data < color_edges[i + 1])\
                                           & (lambda_z > lz_disk)])
                if f_cold[i] > 0:
                    f_cold[i] /= \
                        np.sum(weights[(color_data >= color_edges[i]) \
                                       & (color_data < color_edges[i + 1])])
            # Determine where f_disk crosses 50% (linear interpolation)
            disk_level = 0.5
            cross_idx = np.where(np.diff(np.sign(f_cold - disk_level)))[0]
            if len(cross_idx) > 0:  # does it cross 50% at all?
                idx = cross_idx[0]  # take first crossing
                x1 = c_cold[idx]
                x2 = c_cold[idx + 1]
                y1 = f_cold[idx] - disk_level
                y2 = f_cold[idx + 1] - disk_level
                disk_50 = x1 - y1 * (x2 - x1) / (y2 - y1)
            else:
                disk_50 = np.nan  # no crossing, set to NaN
                self.logger.warning('Cold orbit fraction does not cross 50% '
                    'in the given color range, setting disk_50 to NaN.')
            self.logger.info(f'Cold orbit fraction crosses {disk_level:.2f} '
                            f'at {disk_50:.2f} Gyr.')
        # Plot circularity vs. age
        fig = plt.figure(figsize=(6, 6))
        if disk_fraction:
            ax = fig.add_axes([0.15, 0.1, 0.80, 0.6])
            ax2 = fig.add_axes([0.15, 0.706, 0.6, 0.225])
        else:
            ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
        vmin = np.nanmin(values)
        vmax = np.nanmax(values)
        if p_scale == 'log':
            if vmin < 0:  # == 0 is ok, see below...
                self.logger.warning('Colorbar logscale might not make sense '
                                    'due to negative numbers.')
            vmin = np.nanmin(values[values > 0])
        cax = ax.imshow(values.T,
                        # cmap='Oranges',
                        cmap='gist_heat_r',
                        norm=p_scale,
                        interpolation=interpolation,
                        extent=[min_c, max_c, -1, 1],
                        origin='lower',
                        vmin=vmin,
                        vmax=vmax,
                        aspect='auto')
        ax.set_xlim(min_c, max_c)
        ax.set_yticks([-lz_disk, -0.5, 0, 0.5, lz_disk])
        if c_scale == 'log':
            ax.set_xscale('log')
        ax.set_xlabel(c_label)
        ax.set_ylabel('Circularity $\lambda_z$')
        fig.colorbar(cax,
                     orientation='vertical',
                     pad=0.1,
                     label=r'$r$')
        if disk_fraction:  # Plot cold orbit fraction vs. age
            ax.plot([min_c, max_c], [lz_disk, lz_disk], 'k--')
            ax.plot([disk_50, disk_50], [-1, 1], 'k--')
            if c_scale == 'log':
                ax2.set_xscale('log')
            ax2.set_xticks([])
            ax2.set_yticks([0, 0.5, 1.0])
            ax2.plot([disk_50, disk_50], [0, 1], 'k--')
            ax2.plot([min_c, max_c], [disk_level, disk_level], 'k--')
            ax2.plot(c_cold, f_cold, 'r-')
            ax2.set_ylabel(r'$f_{\mathrm{disk}}$')
            ax2.set_xlim(min_c, max_c)
        # Save the figure
        figname = self.config.settings.io_settings['plot_directory'] + \
            f'circularity_color_plot_{n_models:02d}' + figtype
        fig.savefig(figname, dpi=dpi)
        self.logger.info(f'Circularity color plot saved in {figname}.')
        return fig, disk_50

    def color_color_plot(self,
                         x_posterior,
                         y_posterior,
                         weights=1,
                         *,
                         x_label='Color 1',
                         y_label='Color 2',
                         x_scale='linear',
                         y_scale='linear',
                         n_smooth=100,
                         figtype='.png',
                         dpi=100):

        """Create a color vs color plot with smoothing

        Create and save a plot showing the relation of two colors (e.g. age vs
        metallicity, AMR) for a set of stellar bundles, averaged over multiple
        MCMC chains and draws. The plot includes the probability density
        distribution in the color values and a scatter plot with diamonds
        indicating the average color values of each orbit bundle and the
        diamond sizes being proportional to the orbit bundles' weights. The
        plot is saved in the specified file format and resolution.

        Parameters
        ----------
        x_posterior : xarray or np.array of shape (n_chain, n_draw, n_bundle)
            Posterior distribution of the color plotted on the x axis for the
            orbit bundles, where `n_chain` is the number of MCMC chains,
            `n_draw` is the number of draws per chain, and `n_bundle` is the
            number of orbit bundles.
        y_posterior : xarray or np.array of shape (n_chain, n_draw, n_bundle)
            Posterior distribution of the color plotted on the y axis for the
            orbit bundles, with the same shape as `x_posterior`.
        weights : np.array of shape (n_bundle,) or 1, optional
            Weight of each stellar bundle, used to determine the size of the
            symbols in the scatter plot. If 1, all symbols will have the same
            size. The default is 1.
        x_label : str, optional
            Label for the x axis, by default 'Color 1'.
        y_label : str, optional
            Label for the y axis, by default 'Color 2'.
        x_scale : str, optional
            Scale of the x axis, either 'linear' or 'log', by default 'linear'.
        y_scale : str, optional
            Scale of the y axis, either 'linear' or 'log', by default 'linear
        n_smooth : int, optional
            Number of points for each orbit bundle to smooth the posterior
            distributions over. If `n_smooth` is greater than the number of
            available points, it will be set to the number of available
            points. The default is 100.
        figtype : str, optional
            Determines the file format of the saved figure, by default '.png'.
        dpi : int, optional
            The resolution of saved figure, by default 100.

        Returns
        -------
        matplotlib.figure.Figure
            The created figure object.

        Raises
        ------
        ValueError
            If the shapes of `x_posterior` and `y_posterior` do not match,
            if the length of `weight` does not match the number of bundles
            in the posterior distributions, or if x_scale or y_scale are
            invalid. The error message will indicate the issue.
        """
        self.logger.info(f'Creating "{x_label}" vs "{y_label}" plot with '
                         f'smoothing factor {n_smooth}.')
        if x_posterior.shape != y_posterior.shape:
            txt = 'x and y posterior shapes do not match: ' \
                  f'{x_posterior.shape} vs {y_posterior.shape}.'
            self.logger.error(txt)
            raise ValueError(txt)
        n_chain, n_draw, n_bundle = x_posterior.shape
        if np.ndim(weights) == 0:
            if weights == 1:
                weights = np.ones(n_bundle)
            else:
                txt = f'If weights is a scalar, it must be 1, not {weights}.'
                self.logger.error(txt)
                raise ValueError(txt)
        elif len(weights) != n_bundle:
            txt = f'Weight array length {len(weights)} does not match ' \
                  f'the number of orbit bundles {n_bundle}.'
            self.logger.error(txt)
            raise ValueError(txt)
        if x_scale not in ['linear','log'] or y_scale not in ['linear','log']:
            txt = 'x_scale and y_scale need to be either "linear" or "log", ' \
                  f'but are "{x_scale}" and "{y_scale}".'
            self.logger.error(txt)
            raise ValueError(txt)
        # Smoothing
        idx_list = [(i, j) for i in range(n_chain) for j in range(n_draw)]
        if n_smooth > len(idx_list):
            n_smooth = len(idx_list)
            smooth_list = np.array(idx_list)
            x_smooth = np.array(x_posterior)
            y_smooth = np.array(y_posterior)
            self.logger.info(f'Smoothing requested for {n_smooth} points, '
                             f'but only {len(idx_list)} available, '
                             'using all available points.')
        else:
            rng = np.random.default_rng()
            # Select random values without duplicates
            smooth_idx = rng.choice(len(idx_list),size=n_smooth,replace=False)
            smooth_list = np.array(idx_list)[smooth_idx]
            x_smooth = np.array([x_posterior[c, d, b]
                                for (c, d) in smooth_list
                                for b in range(n_bundle)])
            y_smooth = np.array([y_posterior[c, d, b]
                                for (c, d) in smooth_list
                                for b in range(n_bundle)])
        # Calculate bundle mean values
        x_mean = x_posterior.mean(axis=(0,1))
        y_mean = y_posterior.mean(axis=(0,1))

        fig = plt.figure(figsize=(6, 5))

        min_x = min(np.nanmin(x_smooth), np.nanmin(x_mean))
        max_x = max(np.nanmax(x_smooth), np.nanmax(x_mean))
        min_y = min(np.nanmin(y_smooth), np.nanmin(y_mean))
        max_y = max(np.nanmax(y_smooth), np.nanmax(y_mean))

        # Density plot
        k = scipy.stats.gaussian_kde([x_smooth, y_smooth])
        n_bins = 50  # number of bins for the density plot
        xi, yi = np.mgrid[min_x:max_x:n_bins * 1j, min_y:max_y:n_bins * 1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()])).reshape(xi.shape)
        plt.pcolormesh(xi, yi, zi, shading='gouraud', cmap='gist_heat_r')
        # Optional: scatter plot of all smoothed points
        # weight_smooth = np.array([weights[b]
        #                           for (c, d) in smooth_list
        #                           for b in range(n_bundle)])
        # plt.scatter(age_smooth,
        #             met_smooth,
        #             s=weight_smooth / np.max(weight) * 200,
        #             marker='o',
        #             facecolor='black',
        #             alpha=0.2)
        # Scatter plot of orbital bundle mean values
        weights = weights / np.max(weights) * 200
        scatter_data = [x_mean, y_mean]
        scatter_args = {'marker': 'D', 'facecolor': 'none', 'linewidth': 1}
        plt.scatter(*scatter_data,
                    **scatter_args,
                    s=weights,
                    edgecolor='whitesmoke')
        plt.scatter(*scatter_data,
                    **scatter_args,
                    s=(np.sqrt(weights)-1)**2,
                    edgecolor='red')
        # Axes
        if x_scale == 'log':
            plt.xscale('log')
        if y_scale == 'log':
            plt.yscale('log')
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xlim(min_x, max_x)
        plt.ylim(min_y, max_y)

        # Save the figure
        figname = self.config.settings.io_settings['plot_directory'] + \
            'color_color_plot' + figtype
        fig.savefig(figname, dpi=dpi)
        self.logger.info(f'Color-color plot saved in {figname}.')
        return fig

    def color_maps(self,
                   colors,
                   model_data,
                   flux_norm,
                   cbar_lims='data',
                   figtype='.png',
                   dpi=100):
        """Plot color maps for observed and model data

        Create and save color maps for the observed and model data, along with
        the residuals defined by residual = (model - data) / data_error.
        The method generates a grid of subplots with the individual colors in
        the columns and observed data, model data, and residuals in the three
        rows, respectively. If the provided model data includes errors, the
        observation and model errors will be plotted next to the data. The
        color maps are generated using the `cmasher` library, and the color bar
        limits can be derived from the data or model values, or automatically
        set. The method raises a ValueError if the number of colors does not
        match the number of model data columns, or if the provided `cbar_lims`
        is not one of the expected values.

        Parameters
        ----------
        colors : dict
            Dictionary that defines the colors to be plotted. The keys are the
            color descriptors and must match the column names in the observed
            data (in any order). The values are the descriptions of those
            colors and are used in the plot titles. Note, that unlike for the
            observed data columns, the ORDER MATTERS for the model data:
            model_data needs to provide the data in the same order as the
            colors dictionary. Example: {'age': 'Age [Gyr]',
            'met': 'Metallicity'}.
        model_data : list of np.arrays of shape (n_bundle,)
            List of model data arrays where each array corresponds to a color
            in the `colors` dictionary. The length of the list determines
            whether errors are to be plottetd: If the length matches the
            number of colors in the `colors` dictionary, errors are not
            plotted. If the length is twice the number of colors, errors are
            plotted and error data is expected in every other column.
            Each array in the list should be a 1D numpy array of shape
            (n_bundle,), where `n_bundle` is the number of Voronoi orbit
            bundles in the model. The order of the arrays must match the
            order of the colors in the `colors` dictionary.
            Example with errors plotted: [age, dage, met, dmet], without errors
            plottetd: [age, met]. Each array has shape (n_bundle,).
        flux_norm : np.array of shape (n_spatial_bins, n_bundle)
            Normalized flux data for the spatial bins. Each column corresponds
            to an orbit bundle and each row corresponds to a spatial bin. It is
            normalized such that the sum of fluxes in each spatial bin is equal
            to 1. This is typically the result of orbit bundle maps calculated
            from a Voronoi binning process of the phase space.
        cbar_lims : str, optional
            Determines the limits for the color bar. Can be 'data', 'model',
            or 'auto'. If 'data', the limits are based on the observed data.
            If 'model', the limits are based on the model data. If 'auto',
            the limits adapt for each color. The default is 'data'.
        figtype : str, optional
            Determines the file format of the saved figure, by default '.png'.
        dpi : int, optional
            The resolution of saved figure, by default 100.

        Returns
        -------
        matplotlib.figure.Figure
            The created figure object containing the color maps.

        Raises
        ------
        ValueError
            If the number of colors is not compatible with the number of model
            data columns, or if the provided `cbar_lims` is not one of the
            expected values. The error message will indicate the issue.
        """
        n_colors = len(colors)
        txt = f' for {n_colors} colors: {list(colors.keys())}.'
        if len(model_data) == n_colors:
            self.logger.info('Color maps: plotting data' + txt)
            plot_errors = False
        elif len(model_data) == 2 * n_colors:
            self.logger.info('Color maps: plotting data and errors' + txt)
            plot_errors = True
        else:
            txt = 'Model data must be of same or twice the length of colors '
            txt += '({n_colors} or {2 * n_colors}), not {len(model_data)}.'
            self.logger.error(txt)
            raise ValueError(txt)
        # Generate column names for the colors in the observed data,
        # e.g. 'age', 'dage', 'met', 'dmet'
        if plot_errors:
            plot_cols = [p + col for col in colors for p in ('', 'd')]
        else:
            plot_cols = [col for col in colors]
        # Get the observed data
        if self.config.system.is_bar_disk_system():
            stars = self.config.system.get_unique_bar_component()
        else:
            stars = self.config.system.get_unique_triaxial_visible_component()
        pops = stars.population_data[0]
        map_plotter = pops.get_map_plotter()
        pops_data = pops.get_data()
        # Calculate color bar limits
        if cbar_lims == 'data':
            vmin = [np.min(pops_data[color]) for color in colors]
            vmax = [np.max(pops_data[color]) for color in colors]
        elif cbar_lims == 'model':
            idx = range(0, 2 * n_colors, 2) if plot_errors else range(n_colors)
            tmp = [np.dot(flux_norm, model_data[i_column]) for i_column in idx]
            vmin = [np.min(tmp[i_column]) for i_column in idx]
            vmax = [np.max(tmp[i_column]) for i_column in idx]
        elif cbar_lims == 'auto':
            vmin = vmax = [None] * n_colors
        else:
            txt = 'cbar_lims must be "data", "model", or "auto", not ' \
                  f'{cbar_lims}.'
            self.logger.error(txt)
            raise ValueError(txt)
        if plot_errors:  # adaptive limits for errors
            vmin = [v for b in [(value, None) for value in vmin] for v in b]
            vmax = [v for b in [(value, None) for value in vmax] for v in b]
        # Plot layout
        n_columns = 2 * n_colors if plot_errors else n_colors
        n_rows = 3  # data, model, and residual
        fig = plt.figure(figsize=(5 * n_columns, 2 * n_rows))
        fig.subplots_adjust(wspace=0.3)
        # Color bar and colormap settings
        kw_maps = {'colorbar': True,
                   'cmap': cmasher.get_sub_cmap('twilight_shifted', 0.05, 0.6)}
        # Plot the observed data in the first row
        for i_plot, plot_col in enumerate(plot_cols):
            ax = plt.subplot(n_rows, n_columns, i_plot + 1)
            map_plotter(pops_data[plot_col],
                        vmin=vmin[i_plot],
                        vmax=vmax[i_plot],
                        **kw_maps)
            if plot_errors:
                if i_plot % 2 == 0:
                    title = colors[plot_col]
                else:
                    title = colors[plot_cols[i_plot - 1]] + ' error'
            else:
                title = colors[plot_col]
            ax.set_title(title)
            if i_plot == 0:  # first column
                ax.set_ylabel('Data\n\ny [arcsec]')
        # Plot the model data (and errors if applicable) in the second row
        for i_column in range(n_columns):
            i_plot += 1
            ax = plt.subplot(n_rows, n_columns, i_plot + 1)
            if plot_errors:
                if i_column % 2 == 0:  # data: even columns
                    col_data_aperture = np.dot(flux_norm, model_data[i_column])
                else:  # errors: odd columns
                    col_data_aperture = \
                        np.dot(flux_norm, np.square(model_data[i_column]))
                    col_data_aperture = np.sqrt(col_data_aperture)
                    ax.set_xlabel('x [arcsec]')
            else:
                col_data_aperture = np.dot(flux_norm, model_data[i_column])
            map_plotter(col_data_aperture,
                        vmin=vmin[i_column],
                        vmax=vmax[i_column],
                        **kw_maps)
            if i_column == 0:  # first column
                ax.set_ylabel('Model\n\ny [arcsec]')
        # Plot the residuals in the third row
        for i_column in range(n_columns):
            i_plot += 1
            if (not plot_errors) or i_column % 2 == 0:  # data columns only
                ax = plt.subplot(n_rows, n_columns, i_plot + 1)
                # residual = (model - data) / data_error
                col_data_aperture = (np.dot(flux_norm, model_data[i_column]) \
                                     - pops_data[plot_cols[i_column]] \
                                    ) / pops_data['d' + plot_cols[i_column]]
                map_plotter(col_data_aperture,
                            vmin=None,  # adaptive limits for residuals
                            vmax=None,
                            **kw_maps)
                ax.set_xlabel('x [arcsec]')
                if i_column == 0:  # first column
                    ax.set_ylabel('Residuals\n\ny [arcsec]')
        # Save the figure
        figname = self.config.settings.io_settings['plot_directory'] + \
            'color_maps' + figtype
        fig.savefig(figname, dpi=dpi)
        self.logger.info(f'Color maps saved in {figname}.')
        return fig
