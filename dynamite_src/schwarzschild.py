from . import weight_solvers as ws
from . import dynamics as dyn

class SchwarzschildModelLoop(object):

    def __init__(self,
                 system=None,
                 parset_list=None,
                 weight_solver=None,
                 orblib_pars = {'nE':10, 'nI2':5, 'nI3':5}
                 ):
        self.models = []
        for parset0 in parset_list:
            mod0 = SchwarzschildModel(
                system=system,
                parset=parset0,
                weight_solver=weight_solver,
                orblib_pars = orblib_pars)
            self.models += [mod0]


class SchwarzschildModel(object):

    def __init__(self,
                 system=None,
                 parset=None,
                 weight_solver=None,
                 orblib_pars = {'nE':10, 'nI2':5, 'nI3':5}
                 ):

        potential = dyn.Potential(system, parset)
        orb_lib = dyn.OrbitLibrary(
            potential=potential,
            **orblib_pars)
        weight_solver.set_orb_lib(orb_lib)
        orb_wts, chi2 = weight_solver.solve()

        # do colouring
        orb_labels = []
        for cmp in system.cmp_list:
            for pop_data0 in cmp.population_data:
                orb_labels0 = pop_data0.colouring_recipe(orb_lib,
                                                          orb_wts)
                orb_labels += [orb_labels0]

        # store resuult
        self.parset = parset
        self.orb_lib = orb_lib
        self.orb_wts = orb_wts
        self.chi2 = chi2
        self.orb_labels = orb_labels


# end
