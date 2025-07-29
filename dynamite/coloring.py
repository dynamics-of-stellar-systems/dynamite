
import logging
import os
from turtle import distance
import matplotlib.pyplot as plt
import numpy as np

import yaml
import pymc as pm

import cmasher

from vorbin.voronoi_2d_binning import voronoi_2d_binning

from dynamite import constants


class Coloring:
    """Class to hold coloring-related routines

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
        if minr == 'auto':
            self.minr = 0.0
        if maxr == 'auto':
            arc_kpc = constants.ARC_KPC(config.system.distMPc)
            self.maxr = 10**config.settings.orblib_settings['logrmax']*arc_kpc
        self.nr = nr
        self.nl = nl
        self.logger.info(f'Coloring initialized with minr={self.minr}, '
                         f'maxr={self.maxr}, nr={self.nr}, nl={self.nl}.')

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
            self.logger.info('Voronoi binning metadata saved to '
                             f'{f_metadata}.')
        np.savez(model.directory + f_name,
                 vor_bundle_mapping=vor_bundle_mapping,
                 **phase_space_binning)
        self.logger.info('Voronoi orbit bundle mapping saved to '
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
            This is typically the result of orbit color maps calculated from a
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
        return model, trace

    def color_maps(self, model_data=None, flux_data_rel=None):
        if model_data is None:
            model_data = {}
        n_models = len(model_data)

        # Plot the observed data in the first row
        stars = self.config.system.get_unique_triaxial_visible_component()
        pops = stars.population_data[0]
        age, dage, met, dmet = [pops.get_data()[i]
                                for i in ('age', 'dage', 'met', 'dmet')]
        map_plotter = pops.get_map_plotter()
        n_cols = 4
        n_rows = 1 + n_models
        fig = plt.figure(figsize=(20, 2 * n_rows))
        fig.subplots_adjust(wspace=0.5)
        col_min = [min(a) for a in (age, dage, met, dmet)]
        col_max = [max(a) for a in (age, dage, met, dmet)]
        for i_plot, col in enumerate((age, dage, met, dmet)):
            ax = plt.subplot(n_rows, n_cols, i_plot + 1)
            map_plotter(col,
                        vmin=col_min[i_plot % n_cols],
                        vmax=col_max[i_plot % n_cols],
                        # label=col.name,
                        colorbar=True,
                        cmap=cmasher.get_sub_cmap('twilight_shifted',
                                                  0.05,
                                                  0.6))
            ax.set_title(f'{col.name}')
            if i_plot + 1 > (n_rows - 1) * n_cols:  # last row
                ax.set_xlabel('x [arcsec]')
            if i_plot % n_cols == 0:  # first column
                ax.set_ylabel('Observed\ny [arcsec]')

        # Plot the model data in the remaining rows
        for m_name, m_data in model_data.items():
            for i_col, col in enumerate(m_data):  # even columns are data,
                                                  # odd columns are errors
                i_plot += 1
                if col is None:
                    continue
                ax = plt.subplot(n_rows, n_cols, i_plot + 1)
                if i_col % 2 == 0:  # data
                    col_data_aperture = np.dot(flux_data_rel, col)
                    vmin = col_min[i_plot % n_cols]
                    vmax = col_max[i_plot % n_cols]
                else:  # error
                    col_data_aperture = np.dot(flux_data_rel, np.square(col))
                    col_data_aperture = np.sqrt(col_data_aperture)
                    vmin = vmax = None  # no limits for errors
                map_plotter(col_data_aperture,
                            vmin=vmin,
                            vmax=vmax,
                            colorbar=True,
                            cmap=cmasher.get_sub_cmap('twilight_shifted',
                                                      0.05,
                                                      0.6))
                if i_plot + 1 > (n_rows - 1) * n_cols:  # last row
                    ax.set_xlabel('x [arcsec]')
                if i_plot % n_cols == 0:  # first column
                    ax.set_ylabel(f'{m_name}\ny [arcsec]')
        return fig