import yaml

class Parameter(object):

    def __init__(self,
                 desc=None,
                 lo=None,
                 hi=None,
                 step=None,
                 fixed=False,
                 value=None,
                 minstep=None,
                 LaTeX=None,
                 sformat="%g"):
#        self.name = name
        self.desc = desc
        self.lo = lo
        self.hi = hi
        self.step = step
        self.fixed = fixed
        self.value = value
        self.minstep = minstep
        self.LaTeX = LaTeX
        self.sformat = sformat

class StellarParameters(object):
    def __init__(self, *, name, **kwargs):
        self.q = None
        self.p = None
        self.u = None
        if (name == 'q'):
            self.q = Parameter(**kwargs)
        elif (name == 'p'):
            self.p = Parameter(**kwargs)
        elif (name == 'u'):
            self.u = Parameter(**kwargs)
        else:
            raise ValueError('Unknown stellar parameter, use q, p, or u')

    def add(self, *, name, **kwargs):
        update = False
        if (name == 'q'):
            if (self.q is None):
                update = True
            self.q = Parameter(**kwargs)
        elif (name == 'p'):
            if (self.p is None):
                update = True
            self.p = Parameter(**kwargs)
        elif (name == 'u'):
            if (self.u is None):
                update = True
            self.u = Parameter(**kwargs)
        else:
            raise ValueError('Unknown stellar parameter, use q, p, or u')
        return update

    def validate(self):
        if (self.q is not None and self.p is not None and self.u is not None):
            return True
        else:
            return False


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
