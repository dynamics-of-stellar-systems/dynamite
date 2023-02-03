

class Analysis:

    def __init__(self, config, model):
        pass

    def get_gh_model_kinematic_maps(self, ..., v_sigma_option='fit'):
        if v_sigma_option = 'fit':
            v, sigma = ...
            # for each spatial bin, do some scipy.optimize to fit a Gaussian 
            # to the LOSVD in that bin
         if v_sigma_option = 'moments':
            # extract the orblib from the model
            v = orblib.get_mean_v()
            sigma = orblib.get_sigma_v()
        # extract gh kinematic object
        gh_kins = kinematics.get_gh_coeffs( ... )
        # methods are in https://dynamics.univie.ac.at/dynamite_docs/tutorial_notebooks/6_orbits_and_weights.html
        table = ... # [v, sigma, h3, h4, ...]
        if v_sigma_option = 'fit':
            filename = ...
        if v_sigma_option = 'moments':
            filename = ...
        astropy.save(table, filename)

    def get_projection_tensor_for_orbit_distributions():
        pass

    def get_orbit_distributions():
        pass