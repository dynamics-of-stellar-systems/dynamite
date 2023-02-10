import logging
import numpy as np
from scipy.optimize import curve_fit
import astropy
import dynamite as dyn

class Analysis:

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
                                    v_sigma_option='fit'):
        """
        Generates an astropy table in the model directory that holds the
        model's data for creating Gauss-Hermite kinematic maps:
        v, sigma, h3 ... h<number_GH>
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

        Raises
        ------
        ValueError
            if v_sigma_option is neither 'moments' nor 'fit'.

        Returns
        -------
        Nothing

        """
        if model is None:
            model = self.model
        if kin_set is None:
            kin_set = self.kin_set
        if v_sigma_option not in ['moments', 'fit']:
            txt = f"{v_sigma_option=} but must be either 'fit' or 'moments'."
            self.logger.error(txt)
            raise ValueError(txt)
        orblib = model.get_orblib()
        _ = model.get_weights(orblib)
        weights = model.weights
        orblib.read_losvd_histograms()
        # get all orbits' losvds; orbits_losvd.shape = n_orb,n_vbin,n_aperture
        orbits_losvd = orblib.losvd_histograms[kin_set].y[:,:,]
        # weighted sum of orbits_losvd; model_losvd.shape = 1,n_vbin,n_aperture
        model_losvd = np.dot(orbits_losvd.T, weights).T[np.newaxis,:]
        #model_losvd /= np.sum(model_losvd, 0) # normalisation not necessary
        # calculate v_mean and v_sigma values from the losvd histograms
        # also needed as initial conditions if v_sigma_option=='fit'
        model_losvd_hist = \
            dyn.kinematics.Histogram(xedg=orblib.losvd_histograms[kin_set].xedg,
                                     y=model_losvd,
                                     normalise=False)
        v_mean = np.squeeze(model_losvd_hist.get_mean()) # v from distribution
        v_sigma = np.squeeze(model_losvd_hist.get_sigma()) # sigma from distr.
        if v_sigma_option == 'fit': # fit a gaussian in each aperture
            def gauss(x, a, mean, sigma):
                return a*np.exp(-(x-mean)**2/(2.*sigma**2))
            for aperture in range(model_losvd.shape[-1]):
                p_initial = [1., v_mean[aperture], v_sigma[aperture]]
                p_opt, _ = curve_fit(gauss,
                                     orblib.losvd_histograms[kin_set].x,
                                     model_losvd[0,:,aperture],
                                     p0=p_initial)
                v_mean[aperture] = p_opt[1] # overwrite v_mean from distr
                v_sigma[aperture] = p_opt[2] # overwrite v_sigma from distr
        # calculate the model's gh expansion coefficients
        n_gh = self.config.settings.weight_solver_settings['number_GH']
        gh = dyn.kinematics.GaussHermite()
        model_gh_coefficients = gh.get_gh_expansion_coefficients(
            v_mu=v_mean,
            v_sig=v_sigma,
            vel_hist=model_losvd_hist,
            max_order=n_gh)
        # put the coefficients into an astropy table
        col_names = ['v', 'sigma']
        tab_data = [v_mean, v_sigma]
        if n_gh > 2:
            col_names += [f'h{i}' for i in range(3,n_gh+1)]
            tab_data += list(np.squeeze(model_gh_coefficients)[:,3:].T)
        gh_table = astropy.table.Table(tab_data,
                                       names=col_names,
                                       meta={'v_sigma_option': v_sigma_option})
        f_name = f'{model.directory}gh_kinematics_{v_sigma_option}.ecsv'
        gh_table.write(f_name, format='ascii.ecsv', overwrite=True)

    def get_projection_tensor_for_orbit_distributions(self):
        pass

    def get_orbit_distributions(self):
        pass
