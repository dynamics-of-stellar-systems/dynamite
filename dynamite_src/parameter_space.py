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

    def generate(self, current_models, n_new):
        # placeholder function to generate a list of "n_new" parameters
        # return new_parameter_list
        stop = self.check_stopping_critera()
        if stop:
            return []
        # else ...
        return []

    def check_stopping_critera(self, current_models):
        stop_generic = self.check_generic_stopping_critera()
        stop_specific = self.check_specific_stopping_critera()
        stop = stop_generic or stop_specific
        return stop

    def check_generic_stopping_critera(self, current_models):
        stop = True # or false
        return stop


class GridSearch(ParameterGenerator):

    def generate(self, current_models=None, n):
        # actual code to do grid search
        # return new_parameter_list of length n_new
        return []

    def check_specific_stopping_critera(self, current_models):
        stop = True # or false
        return stop


class GaussianProcessEmulator(ParameterGenerator):

    def generate(self, current_models, n):
        # actual code to do gaussian process emulation
        # return new_parameter_list of length n_new
        return []

    def check_specific_stopping_critera(self, current_models):
        stop = True # or false
        return stop


# end
