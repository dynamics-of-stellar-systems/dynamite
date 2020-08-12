from . import weight_solvers as ws
from . import dynamics as dyn

class SchwarzschildModelLoop(object):

    def __init__(self,
                 system=None,
                 parset_list=None,
                 weight_solver=None,
                 config=None
                 ):
        self.models = []
        for parset0 in parset_list:
            if config.legacy_mode:
                mod0 = LegacySchwarzschildModel(
                    system=system,
                    parset=parset0,
                    weight_solver=weight_solver,
                    orblib_pars = orblib_pars)
            else:
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


class LegacySchwarzschildModel(object):

    def __init__(self,
                 system=None,
                 config=None,
                 parset=None):
        self.system = system
        self.config = config
        self.parset = parset

        mod_dir = self.get_model_directory()

        orb_lib = dyn.LegacyOrbitLibrary(
            system=system,
            settings=config.orblib_settings,
            parset=parset)

        # In principle, a weight solver needs (i) orb_lib, (ii) the data (which
        # is stored the system object) and (iii) weight_solver_settings.
        # The LegacyWeightSolver however accesses these quantities all through
        # files which are saved in mod_dir - so we can juct pass mod_dir
        weight_solver = ws.LegacyWeightSolver(
            mod_dir=mod_dir,
            settings=config.weight_solver_settings)
        orb_wts, chi2 = weight_solver.solve()

        # store result
        self.orb_wts = orb_wts
        self.chi2 = chi2

    def get_model_directory(self):
        out_dir = self.config.output_settings['directory']
        print(out_dir)
        return out_dir


# end
