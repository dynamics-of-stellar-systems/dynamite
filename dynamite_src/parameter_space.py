class Parameter(object):

    attributes = []
    def __init__(self,
                 name=None,
                 desc=None,
                 fixed=False,
                 LaTeX=None,
                 sformat="%g",
                 value=None,
                 grid_parspace_settings=None,
                 gpe_parspace_settings=None
                 # lo=None,
                 # hi=None,
                 # step=None,
                 # minstep=None,
                 ):
        self.name = name
        self.desc = desc
        self.fixed = fixed
        self.LaTeX = LaTeX
        self.sformat = sformat
        self.value = value
        self.grid_parspace_settings = grid_parspace_settings
        self.gpe_parspace_settings = gpe_parspace_settings
        # self.lo = lo
        # self.hi = hi
        # self.step = step
        # self.minstep = minstep
        self.__class__.attributes = list(self.__dict__.keys())

    def update(self, **kwargs):
        for k, v in kwargs.items():
            if k not in self.__class__.attributes:
                raise ValueError(f'Invalid parameter key {k}. Allowed keys: {str(tuple(self.__class__.attributes))}')
            setattr(self, k, v)

    def validate(self):
        if sorted(self.__class__.attributes) != sorted(self.__dict__.keys()):
            raise ValueError(f'Parameter attributes can only be {str(tuple(self.__class__.attributes))} , not {str(tuple(self.__dict__.keys()))}')

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.__dict__})')


# class BlackHoleParameters(object):
# #    def __init__(self, *, name, **kwargs):
#     def __init__(self, name=None, **kwargs):
#         self.mass = None
#         self.radius = None
#         if name is not None:
#             self.add(self, name, **kwargs)

# #    def add(self, *, name, **kwargs):
#     def add(self, name, **kwargs):
#         if name == 'mass':
#             self.mass = Parameter(**kwargs)
#         elif name == 'radius':
#             self.radius = Parameter(**kwargs)
#         else:
#             raise ValueError('Unknown black hole parameter ' + name + ', use mass or radius')

#     def validate(self):
#         if self.mass is None or self.radius is None:
#             raise ValueError('Black hole parameters require mass and radius')


# class DarkHaloParameters(object):
#     def __init__(self, name=None, **kwargs):
#         self.dc = None
#         self.f = None
#         if name is not None:
#             self.add(self, name, **kwargs)

# #    def add(self, *, name, **kwargs):
#     def add(self, name, **kwargs):
#         if name == 'dc':
#             self.dc = Parameter(**kwargs)
#         elif name == 'f':
#             self.f = Parameter(**kwargs)
#         else:
#             raise ValueError('Unknown dark matter parameter ' + name + ', use dc or f')

#     def validate(self):
#         if self.dc is None or self.f is None:
#             raise ValueError('Dark matter parameters require dc and f')


# class StellarParameters(object):
# #    def __init__(self, *, name, **kwargs):
#     def __init__(self, name=None, **kwargs):
#         self.q = None
#         self.p = None
#         self.u = None
#         if name is not None:
#             self.add(self, name, **kwargs)
#         # if (name == 'q'):
#         #     self.q = Parameter(**kwargs)
#         # elif (name == 'p'):
#         #     self.p = Parameter(**kwargs)
#         # elif (name == 'u'):
#         #     self.u = Parameter(**kwargs)
#         # else:
#         #     raise ValueError('Unknown stellar parameter ', name, ', use q, p, or u')

# #    def add(self, *, name, **kwargs):
#     def add(self, name, **kwargs):
#         update = False
#         if (name == 'q'):
#             if (self.q is None):
#                 update = True
#             self.q = Parameter(**kwargs)
#         elif (name == 'p'):
#             if (self.p is None):
#                 update = True
#             self.p = Parameter(**kwargs)
#         elif (name == 'u'):
#             if (self.u is None):
#                 update = True
#             self.u = Parameter(**kwargs)
#         else:
#             raise ValueError('Unknown stellar parameter ', name, ', use q, p, or u')
#         return update

#     def validate(self):
#         if (self.q is not None and self.p is not None and self.u is not None):
#             return True
#         else:
#             return False


class ParameterSpace(list):

    def __init__(self, system):
        for cmp in system.cmp_list:
            for par in cmp.parameters:
                self.append(par)
        for par in system.parameters:
            self.append(par)

        self.par_names = []
        for par in self:
            self.par_names.append(par.name)

        self.n_par = len(self)
        self.n_par_fixed = len([p for p in self if p.fixed])
        self.n_par_free = self.n_par - self.n_par_fixed

    def __repr__(self):
        return (f'{self.__class__.__name__}({[p for p in self]}, {self.__dict__})')


class ParameterGenerator(object):

    def __init__(self,
                 par_space=[]):
        self.par_space = par_space

    def generate(self, current_models, n_new):
        # placeholder function to generate a list of "n_new" parameters
        # return new_parameter_list
        # current_models will be an AllModels object
        # ... i.e. all_mod = AllModels(...)
        #          all_mod.tables is the table with params and chi2
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
        # e.g. stop if:
        # i) number of models which have been run > max_n_mods
        # ii) number of iterations > max_n_iter
        # iii) ...
        return stop


class GridSearch(ParameterGenerator):

    def generate(self, current_models=None, n_new=0):
        # actual code to do grid search
        # return new_parameter_list of length n_new
        # if iter == 0 ... or... if there are no models yet:
        #   make a basic grid


        # ...
        #  if min step size has been reached, then fix that param ...
        #
        return []

    def check_specific_stopping_critera(self, current_models):
        stop = True # or false
        # stop if...
        # (i) if iter>1, last iteration did not improve chi2 by min_delta_chi2
        # where we'll set min_delta_chi2 in config file
        # (ii) if step_size < min_step_size for all params
        return stop


class GaussianProcessEmulator(ParameterGenerator):

    def generate(self, current_models=None, n_new=0):
        # actual code to do gaussian process emulation
        # return new_parameter_list of length n_new
        return []

    def check_specific_stopping_critera(self, current_models):
        stop = True # or false
        return stop


# end
