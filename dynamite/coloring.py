
import logging
import matplotlib.pyplot as plt
import numpy as np

from vorbin.voronoi_2d_binning import voronoi_2d_binning


class Coloring:
    """Class to hold coloring-related routines

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object

    """

    def __init__(self, config):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.config = config

    def bin_phase_space(self,
                        model=None,
                        minr=None,
                        maxr=None,
                        nr=50,
                        nl=61,
                        vor_weight=0.01,
                        vor_ignore_zeros=False,
                        make_diagnostic_plots=False,
                        extra_diagnostic_output=False):
        """
        Perform Voronoi binning of orbits in the radius-circularity phase
        space. The goal is to group the "original" n_orbits orbit bundles into
        fewer n_bundle "Voronoi" bundles with each of these Voronoi bundles
        accounting for a weight of at least ``vor_weight``.
        This method uses the ``vorbin`` package by Michele Cappellari
        (M. Cappellari, Y. Copin, 2003, MNRAS 342, 345).

        Parameters
        ----------
        model : a ``dynamite.model.Model`` object, optional
            the model used for the binning. If ``None``, choose the minimum
            :math:`\\chi^2` model (the configuration setting ``which_chi2``
            determines the :math:`\\chi^2` type). The default is ``None``.
        minr : float, optional
            the minimum radius [kpc] considered in the binning. If ``None``,
            this is set to the minimum radius of the orbit library.
            The default is ``None.``
        maxr : float, optional
            the maximum radius [kpc] considered in the binning. If ``None``,
            this is set to the maximum radius of the orbit library.
            The default is ``None``.
        nr : int, optional
            number of radial bins. The default is 50.
        nl : int, optional
            number of circularity bins. The default is 61.
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

        Returns
        -------
        (vor_bundle_mapping, phase_space_binning) : tuple, where
        vor_bundle_mapping : np.array of shape (n_bundle, n_orbits)
            Mapping between the "original" orbit bundles and the Voronoi
            orbit bundles: vor_bundle_mapping(i_bundle, i_orbit) is the
            fraction of i_orbit assigned to i_bundle, multiplied by i_orbit's
            weight.
        phase_space_binning : dict
            phase_space_binning['in']: np.array of shape (3, nr*nl) holding the
            binning input: bin r coordinates, bin lambda_z coordinates,
            bin total weights
            phase_space_binning['out']: np.array of shape (3, n_bundle) holding
            the Voronoi binning output: weighted Voronoi bin centroid
            coordinates r_bar, lambda_bar and Voronoi bin total weights
            phase_space_binning['map']: np.array of shape (nr*nl,) holding the
            phase space mapping: Voronoi bin numbers for each input bin

        """
        if model is None:
            best_model_idx = self.config.all_models.get_best_n_models_idx(1)[0]
            model = self.config.all_models.get_model_from_row(best_model_idx)
        orblib = model.get_orblib()
        _ = model.get_weights(orblib)
        weights = model.weights
        if minr is None or maxr is None:
            # use the orbit library's limits for the radius limits
            if not hasattr(orblib, 'orb_properties'):
                orblib.read_orbit_property_file()
            orb_properties = orblib.orb_properties
            if minr is None:
                minr = np.min(orb_properties['r']).value
            if maxr is None:
                maxr = np.max(orb_properties['r']).value
        log_minr, log_maxr = np.log10(minr), np.log10(maxr)

        orblib.get_projection_tensor(minr=minr,
                                     maxr=maxr,
                                     nr=nr,
                                     nl=nl,
                                     force_lambda_z=True)
        # pick entry [2] = fraction of orbits in each r, l bin (ALL orbit types)
        # indices r, l, b = radius, lambda_z, original_orbit_bundle
        projection_tensor = orblib.projection_tensor[2]

        # get the orbit distribution in the r, lambda space
        # = total weight in the r, l bins
        orbit_distribution = projection_tensor.dot(weights)

        # build the input for vorbin
        dr = (log_maxr - log_minr) / nr
        input_bins_r = np.linspace(log_minr + dr / 2, log_maxr - dr / 2, num=nr)
        dl = 2 / nl
        input_bins_l = np.linspace(-1 + dl / 2, 1 - dl / 2, num=nl)
        input_grid = np.meshgrid(input_bins_r,input_bins_l)
        vor_in = np.vstack((np.array([m.ravel() for m in input_grid]),
                            orbit_distribution.T.ravel()))
        if vor_ignore_zeros:
            nonzero_weight_bins = vor_in[2,:] > 0
            vor_in = vor_in[:,nonzero_weight_bins]

        def sum_weights(index, signal, noise):
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
                vector of length M>N with the noise of all r,l bins, ignored.

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
                               sn_func=sum_weights,
                               quiet=not extra_diagnostic_output)
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
            vor_map_dense = np.array(vor_map_dense).reshape(nl, nr)
            plt.figure()
            plt.imshow(orbit_distribution.T,
                       interpolation='spline16',
                       origin='lower',
                       extent=(log_minr, log_maxr, -1, 1),
                       aspect='auto',
                       cmap='Greys',
                       alpha=0.9)
            plt.pcolormesh(*input_grid,
                           vor_map_dense,
                           shading='nearest',
                           cmap='tab20',
                           alpha=0.3)
            plt.xlabel('$\\log_{10}$ ($r$/kpc)')
            plt.ylabel('Circularity $\\lambda_z$')
            plt.colorbar(label='Voronoi bin id')
            plt.savefig(f'{self.config.settings.io_settings["plot_directory"]}'
                        'coloring_voronoi_binning.png')

        # map "original" orbit bundles to Voronoi bundles
        # calculate weighted orbit fractions in the input bins: "flatten"
        # r and l indices of the projection tensor making r the fastest moving
        # index, like in phase_space_binning['map']==binning and multiply
        # element-wise with the orbit weights; orb_frac.shape=(nr*nl, n_orbits)
        orb_frac = projection_tensor.todense().reshape((nr * nl, -1),
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

        return vor_bundle_mapping, phase_space_binning
