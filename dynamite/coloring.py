
import logging
import matplotlib.pyplot as plt
import numpy as np

from vorbin.voronoi_2d_binning import voronoi_2d_binning

import dynamite as dyn


class Coloring:

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
                        vor_ignore_zeros=True,
                        make_diagnostic_plots=False,
                        extra_diagnostic_output=False):
        """


        Parameters
        ----------
        model : TYPE, optional
            DESCRIPTION. The default is None.
        minr : TYPE, optional
            DESCRIPTION. The default is None.
        maxr : TYPE, optional
            DESCRIPTION. The default is None.
        nr : TYPE, optional
            DESCRIPTION. The default is 50.
        nl : TYPE, optional
            DESCRIPTION. The default is 61.
        vor_weight : TYPE, optional
            DESCRIPTION. The default is 0.01.
        vor_ignore_zeros : TYPE, optional
            DESCRIPTION. The default is True.
        make_diagnostic_plots : TYPE, optional
            DESCRIPTION. The default is False.
        extra_diagnostic_output : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        vor_bundle_mapping : TYPE
            DESCRIPTION.
        phase_space_binning : TYPE
            DESCRIPTION.

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
        # indices r, l, b = radius, lambda_z, orbit_bundle
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
            vor_in = vor_in[:,vor_in[2,:] > 0]  # eliminate zero-weight bins
            # voronoi_2d_binning(x=vor_in[0],
            #                    y=vor_in[1],
            #                    signal=vor_in[2],
        binning, _, _, r_bar, l_bar, bin_weights, _, _ = \
            voronoi_2d_binning(*vor_in,
                               noise=np.ones_like(vor_in[0]),
                               target_sn=vor_weight,
                               plot=make_diagnostic_plots,
                               quiet=not extra_diagnostic_output)
        phase_space_binning = {'in':vor_in,
                               'out':np.vstack((r_bar, l_bar, bin_weights)),
                               'map':binning}  # Vor bin numbers of r, l bins
        #log: 'x original orbit bundles divided into y Voronoi bundles'
        if make_diagnostic_plots:
            if vor_ignore_zeros:
                weight_nonzero = orbit_distribution.T.ravel() > 0
                vor_bin_map = iter(phase_space_binning['map'])
                vor_map_dense = [next(vor_bin_map) if nz else np.nan
                                 for nz in weight_nonzero]
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

        # map "original" orbit bundles to voronoi bundles
        # flatten the r and l indices of the projection tensor and make
        # r the fastest moving index, like in binning
        orbit_fractions = projection_tensor.todense().reshape((nr * nl, -1),
                                                              order='F')
        # each Voronoi bin number represents a Voronoi (orbit) bundle
        # for each "original" (orbit) bundle, calculate its fractions that go
        # into each Voronoi bundle
        n_orbits = projection_tensor.shape[-1]  # "original" orbit bundles
        n_bundle = phase_space_binning['out'].shape[-1]  # Voronoi bundles
        vor_bundle_mapping = np.zeros((n_bundle, n_orbits))
        vor_bundle_mapping[binning[:]] = orbit_fractions[:]

        return vor_bundle_mapping, phase_space_binning
