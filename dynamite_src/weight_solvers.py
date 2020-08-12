from . import dynamics

class WeightSolver(object):

    def __init__(self,
                 weight_solver_args={}
                 ):
        self.weight_solver_args = weight_solver_args

    def set_kinematics(self, system):
        self.kinematics = []
        for component in system.cmp_list:
            self.kinematics += component.kinematic_data

    def set_orb_lib(self, orb_lib):
        self.orb_lib = orb_lib

    def get_observables_from_orbits(self):
        self.orb_observed = []
        for kin_data0 in self.kinematics:
            orb_obs0 = kin_data0.transform_orbits_to_observables(self.orb_lib)
            self.orb_observed += [orb_obs0]

    def solve(self):
        # placeholder function to solve for weights given
        # self.kinematics.values and self.orb_observed
        # return wts, chi2
        return 0, 0


class NNLS(WeightSolver):

    def __init__(self,
                 **kwargs):
        super(NNLS, self).__init__(*kwargs)     # initialise parent class

    def solve(self):
        # actual code to do NNLS
        # return MAP_weight
        return 0, 0


class LatentPosterior(WeightSolver):

    def __init__(self,
                 **kwargs):
        super(LatentPosterior, self).__init__(*kwargs) # initialise parent class

    def solve(self):
        # actual code to sample posterior of weights
        return 0, 0


class LegacyWeightSolver(WeightSolver):

    def __init__(self,
                 mod_dir=None,
                 settings=None):
        self.mod_dir = mod_dir
        self.settings = settings

    def solve(self):
        orb_wts, chi2 = 0., 0.
        return orb_wts, chi2







# end
