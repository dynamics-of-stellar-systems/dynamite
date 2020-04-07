class Parameter(object):

    def __init__(self,
                 name=None,
                 lo=None,
                 hi=None,
                 step=None,
                 fixed=False,
                 value=None):
        self.name = name
        self.lo = lo
        self.hi = hi
        self.step = step
        self.fixed = fixed
        self.value = value


class ParameterSpace(object):

    def __init__(self, system):

        self.n_par = system.n_par
        self.n_par_fixed = 0
        for cmp in system.cmp_list:
            for par in cmp.parameters:
                self.n_par_fixed += par.fixed
        self.n_par_free = self.n_par - self.n_par_fixed


class ParameterGenerator(object):

    def __init__(self,
                 param_list=[]):
        self.param_list = param_list

    def generate(self, current_models):
        # placeholder function
        # return new_parameter_list
        return []


class GridSearch(ParameterGenerator):

    def generate(self, current_models=None):
        # actual code to do grid search
        # return new_parameter_list
        return []


class GaussianProcessEmulator(ParameterGenerator):

    def generate(self, current_models):
        # actual code to do gaussian process emulation
        # return new_parameter_list
        return []



# end
