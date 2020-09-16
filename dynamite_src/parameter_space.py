# # some tricks to add the current path to sys.path (so the imports below work)

import os.path
import sys

this_dir = os.path.dirname(__file__)
if not this_dir in sys.path:
    sys.path.append(this_dir)

import numpy as np
import copy
# import schwarzschild
#from astropy import table
import parameter_space as parspace

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
                 gpe_parspace_settings=None,
                 logarithmic=False,
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
        self.logarithmic = logarithmic
        # self.lo = lo
        # self.hi = hi
        # self.step = step
        # self.minstep = minstep
        self.__class__.attributes = list(self.__dict__.keys())

    def update(self, **kwargs):
        for k, v in kwargs.items():
            if k not in self.__class__.attributes:
                raise ValueError(f'Invalid parameter key {k}. Allowed keys: '
                                 f'{str(tuple(self.__class__.attributes))}')
            setattr(self, k, v)

    def validate(self):
        if sorted(self.__class__.attributes) != sorted(self.__dict__.keys()):
            raise ValueError(f'Parameter attributes can only be '
                             f'{str(tuple(self.__class__.attributes))}, '
                             f'not {str(tuple(self.__dict__.keys()))}')

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.__dict__})')

    def get_par_value_from_raw_value(self, raw_value):
        if self.logarithmic is True:
            par_value = 10.**raw_value
        else:
            par_value = raw_value
        return par_value

    def get_raw_value_from_par_value(self, par_value):
        if self.logarithmic is True:
            raw_value = np.log10(par_value)
        else:
            raw_value = par_value
        return raw_value


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
        return (f'{self.__class__.__name__}({[p for p in self]}, '
                f'{self.__dict__})')

    def get_param_value_from_raw_value(self, raw_value):
        par_val = [p.get_par_value_from_raw_value(rv0)
                   for (rv0, p) in zip(raw_value, self)]
        return par_val

    def get_raw_value_from_param_value(self, par_val):
        raw_value = [p.get_raw_value_from_par_value(pv0)
                     for (pv0, p) in zip(par_val, self)]
        return raw_value


class ParameterGenerator(object):

    def __init__(self,
                 par_space=[],
                 parspace_settings=None,
                 name=None):
        self.par_space = par_space
        if not parspace_settings:
            raise ValueError('ParameterGenerator needs parspace_settings')
        self.parspace_settings = parspace_settings
        self.status = {}
        self.name = name

    def generate(self,
                 current_models=None,
                 kw_specific_generate_method={}):
        """Generate new parameter sets. This is a wrapper method around the
        specific_generate_method of child generator classes. This wrapper does
        the following:
        (i) evaluates stopping criteria, and stop if necessary
        (ii) runs the specific_generate_method of the child class, which
        updates self.model_list with a list of propsal models
        (iii) removes previously run models from self.model_list
        (iv) converts parameters from raw_values to par_values
        (v) adds new models to current_models.table
        (vi) update and return the status dictionary

        Parameters
        ----------
        current_models : schwarzschild.AllModels
        kw_specific_generate_method : dict
            keyword arguments passed to the specific_generate_method of the
            child class

        Returns
        -------
        dict
            a dictionary status, with entries
            - stop: bool, whether or not any stopping criteria are met
            - n_new_models: int
            and additional bool entries for the indivdidual stopping criteria
            - last_iter_added_no_new_models
            - n_max_mods_reached
            - n_max_iter_reached
            - plus any criteria specific to the child class

        """
        if current_models is None:
            errormsg = "current_models needs to be a valid " \
                       "schwarzschild.AllModels instance"
            raise ValueError(errormsg)
        else:
            self.current_models = current_models
        self.check_stopping_critera()
        if len(self.current_models.table)==0:
            this_iter = 0
        else:
            this_iter = np.max(self.current_models.table['which_iter']) + 1
        # check whether we need to do anything in the first place...
        if not self.status['stop']:
            self.specific_generate_method(**kw_specific_generate_method)
            newmodels = 0
            # Add new models to current_models.table
            for m in self.model_list:
                if self._is_newmodel(m, eps=1e-10):
                    self.add_model(m, n_iter=this_iter)
                    newmodels += 1
        print(f'{self.name} added {newmodels} new model(s) out of '
              f'{len(self.model_list)}')
        self.status['n_new_models'] = newmodels
        last_iter_check = True if newmodels == 0 else False
        self.status['last_iter_added_no_new_models'] = last_iter_check
        if newmodels == 0:
            self.status['stop'] = True
        return self.status

    def add_model(self, model=None, n_iter=0):
        """
        Adds the model passed as an argument to self.current_models.
        The datetime64 column is populated with the current timestamp
        numpy.datetime64('now', 'ms').
        The 'which_iter' column is populated with the argument value n_iter.

        Parameters
        ----------
        model : List of Parameter objects
        n_iter : integer
            value to write in 'which_iter' column, optional. The default is 0.

        Raises
        ------
        ValueError
            If no or empty model is given.

        Returns
        -------
        None.

        """
        if not model:
            raise ValueError('No or empty model')
        raw_row = [p.value for p in model]
        row = self.par_space.get_param_value_from_raw_value(raw_row)
        # for all columns after parameters, add a entry to this row
        idx_start = self.par_space.n_par
        idx_end = len(self.current_models.table.colnames)
        for i in range(idx_start, idx_end):
            if self.current_models.table.columns[i].dtype == '<M8[ms]':
                # current time
                val = np.datetime64('now', 'ms')
            elif self.current_models.table.columns[i].name == 'which_iter':
                # iteration
                val = n_iter
            else:
                # empty entry for all other columns
                dtype = self.current_models.table.columns[i].dtype
                val = np.dtype(dtype).type(0)
            row.append(val)
        self.current_models.table.add_row(row)

    def check_stopping_critera(self):
        self.status['stop'] = False
        if len(self.current_models.table) > 0:
        # never stop when current_models is empty
            self.check_generic_stopping_critera()
            self.check_specific_stopping_critera()
            for key in [reasons for reasons in self.status \
                        if isinstance(reasons, bool) and reasons != 'stop']:
                if self.status[key]:
                    self.status['stop'] = True
                    break

    def check_generic_stopping_critera(self):
        self.status['n_max_mods_reached'] = \
            len(self.current_models.table) \
                >= self.parspace_settings['stopping_criteria']['n_max_mods']
        self.status['n_max_iter_reached'] = \
            np.max(self.current_models.table['which_iter']) \
                >= self.parspace_settings['stopping_criteria']['n_max_iter']
        # iii) ...

    def _is_newmodel(self, model, eps=1e-6):
        """
        Checks if model is a new model (i.e., its parameter set does not exist
        in self.current_models).

        Parameters
        ----------
        model : A self.model_list element (list of Parameter objects),
                mandatory
        eps : Used for numerical comparison (relative difference w.r.t.
              model values), default is 1e-6

        Returns
        -------
        isnew : True if model is a new model, False otherwise.

        """
        if any(map(lambda t: not isinstance(t, parspace.Parameter), model)):
            raise ValueError('Model arg. must be list of Parameter objects')
        raw_model_values = [p.value for p in model]
        model_values = \
            self.par_space.get_param_value_from_raw_value(raw_model_values)
        if len(self.current_models.table) > 0:
            isnew = True
            for curmod in self.current_models.table[self.par_space.par_names]:
                if np.allclose(list(curmod), model_values, rtol=eps):
                    isnew = False
                    break
        else:
            isnew = True
        return isnew


class LegacyGridSearch(ParameterGenerator):

    def __init__(self, par_space=[], parspace_settings=None):
        super().__init__(par_space=par_space,
                         parspace_settings=parspace_settings,
                         name='LegacyGridSearch')

    def specific_generate_method(self, **kwargs):
        if len(self.current_models.table) == 0:
            # The 'zeroth iteration' results in only one model
            # (all parameters at their .value level)
            self.model_list = [[p for p in self.par_space]]
            return ###########################################################

        thresh = \
            self.parspace_settings['generator_settings']['threshold_del_chi2']
        chi2 = 'chi2' if self.parspace_settings['which_chi2'] == 'chi2' \
            else 'kinchi2'
        min_chi2 = np.min(self.current_models.table[chi2])
        prop_mask = abs(self.current_models.table[chi2] - min_chi2) <= thresh
        prop_list = self.current_models.table[prop_mask]
        self.model_list = []
        for m in prop_list:
            for paridx in range(len(self.par_space)):
                if self.par_space[paridx].fixed: # parameter fixed->do nothing
                    continue
                new_parset = [copy.deepcopy(p) for p in self.par_space]
                lo = self.par_space[paridx].grid_parspace_settings['lo']
                hi = self.par_space[paridx].grid_parspace_settings['hi']
                step = self.par_space[paridx].grid_parspace_settings['step']
                minstep = \
                    self.par_space[paridx].grid_parspace_settings['minstep'] \
                    if 'minstep' in \
                    self.par_space[paridx].grid_parspace_settings else step
                for s in [-1, 1]:
                    value = np.clip(m[paridx].columns[paridx]+s*step, lo, hi)
                    if abs(value-m[paridx].columns[paridx]) >= minstep:
                        new_parset[paridx].value = value
                        self.model_list.append(new_parset)
                    


class GridSearch(ParameterGenerator):

    def __init__(self,
                 par_space=[],
                 parspace_settings=None):
        super().__init__(par_space=par_space,
                         parspace_settings=parspace_settings,
                         name='GridWalk')

    def specific_generate_method(self, **kwargs):
        """
        Generates list of new models self.model_list. Each element of
        self.model_list is a list of Parameter objects. The center of the
        grid search is the parameter set with the smallest chi2+kinchi2 value.

        Parameters
        ----------
        None.

        Raises
        ------
        None.

        Returns
        -------
        None. self.model_list is the list of new models.
        """
        if len(self.current_models.table) == 0:
            # The 'zeroth iteration' results in only one model
            # (all parameters at their .value level)
            self.model_list = [[p for p in self.par_space]]
        else: # Subsequent iterations...
            # center criterion: min(chi2) or min(kinchi2)
            chi2 = 'chi2' if self.parspace_settings['which_chi2'] == 'chi2' \
                else 'kinchi2'
            center_idx = np.argmin(self.current_models.table[chi2])
            n_par = self.par_space.n_par
            center = list(self.current_models.table[center_idx])[:n_par]
            raw_center = self.par_space.get_raw_value_from_param_value(center)
            # print(f'center: {center}')
            # Build model_list by walking the grid
            self.model_list = []
            self.grid_walk(center=raw_center)
            # for m in self.model_list:
            #     print(f'{[(p.name, p.value) for p in m]}')
        return

    def grid_walk(self, center=None, par=None, eps=1e-6):
        """
        Walks the grid defined by self.par_space.grid_parspace_settings
        attributes.
        Clips parameter values to lo/hi attributes. If clipping violates the
        minstep attribute, the resulting model(s) will not be created. If the
        minstep attribute is missing, the step attribute will be used instead.
        Use minstep=0 to eliminate minstep.

        Parameters
        ----------
        center : List of center coordinates. Must be in the same sequence as
                 the parameters in self.par_space. Mandatory argument.
        par : Internal use only. Gives the parameter to start with. Set
              automatically in the recursive process. The default is None.
        eps : Used for numerical comparison (relative tolerance), default 1e-6

        Raises
        ------
        ValueError if center is not specified or fixed parameters != center.

        Returns
        -------
        None. Sets self.model_list to the resulting models.

        """
        if center == None:
            raise ValueError('Need center')
        if not par:
            par = self.par_space[0]
        paridx = self.par_space.index(par)
        # print(f'Call with paridx={paridx}, n_par={self.par_space.n_par}')

        if par.fixed:
            par_values = [par.value]
            if abs(center[paridx] - par.value) > eps:
                raise ValueError('Something is wrong: fixed parameter value '
                                 'not in center')
        else:
            lo = par.grid_parspace_settings['lo']
            hi = par.grid_parspace_settings['hi']
            step = par.grid_parspace_settings['step']
            # up to 3 *distinct* par_values (clipped lo, mid, hi values)
            par_values = []
            # use 'minstep' value if present, otherwise use 'step'
            minstep = par.grid_parspace_settings['minstep'] \
                if 'minstep' in par.grid_parspace_settings else step
            # start with lo...
            delta = center[paridx] - self.clip(center[paridx] - step, lo, hi)
            if abs(delta) >= minstep:
                par_values.append(self.clip(center[paridx] - step, lo, hi))
            # now mid... tol(erance) is necessary in case minstep < eps
            if len(par_values) > 0:
                # check for values differing by more than eps...
                tol = abs(self.clip(center[paridx], lo, hi) - par_values[0])
                if abs(par_values[0]) > eps: # relative tolerance usable=?
                    tol /= abs(par_values[0])
                if tol > eps:
                    par_values.append(self.clip(center[paridx], lo, hi))
            else:
                par_values.append(self.clip(center[paridx], lo, hi))
            # and now hi...
            delta = self.clip(center[paridx] + step, lo, hi) - center[paridx]
            if abs(delta) >= minstep:
                tol = abs(self.clip(center[paridx]+step,lo,hi)-par_values[-1])
                if abs(par_values[-1]) > eps:
                    tol /= abs(par_values[-1])
                if tol > eps:
                    par_values.append(self.clip(center[paridx]+step, lo, hi))

        for value in par_values:
            parcpy = copy.deepcopy(par)
            parcpy.value = value
            if not self.model_list: # add first entry if model_list is empty
                self.model_list = [[parcpy]]
                models_prev = [[]]
                # print(f'new model list, starting w/parameter {parcpy.name}')
            elif parcpy.name in [p.name for p in self.model_list[0]]:
                # in this case, create new (partial) model by copying last
                # models and setting the new parameter value
                for m in models_prev:
                    new_model = m + [parcpy]
                    self.model_list.append(new_model)
                # print(f'{parcpy.name} is in '
                #       f'{[p.name for p in self.model_list[0]]}, '
                #       f'added {parcpy.name}={parcpy.value}')
            else: # new parameter: append it to existing (partial) models
                models_prev = copy.deepcopy(self.model_list)
                for m in self.model_list:
                    m.append(parcpy)
                # print(f'new parameter {parcpy.name}={parcpy.value}')

        # call recursively until all paramaters are done:
        if paridx < self.par_space.n_par - 1:
            self.grid_walk(center=center, par=self.par_space[paridx+1])

    def clip(self, value, mini, maxi):
        """
        Clips value to the interval [mini, maxi]. Similar to the numpy.clip()
        method. If mini==maxi, that value is returned.

        Parameters
        ----------
        value : numeric value
        mini : numeric value
        maxi : numeric value

        Raises
        ------
        ValueError if mini > maxi

        Returns
        -------
        min(max(mini, value), maxi)

        """
        if mini <= maxi:
            return min(max(mini, value), maxi)
        else:
            raise ValueError('Clip error: minimum must be less than '
                             'or equal to maximum')

    def check_specific_stopping_critera(self):
        # stop if...
        # (i) if iter>1, last iteration did not improve chi2 by min_delta_chi2
        self.status['min_delta_chi2_reached'] = False
        last_iter = np.max(self.current_models.table['which_iter'])
        if last_iter > 0:
            mask = self.current_models.table['which_iter'] == last_iter
            models0 = self.current_models.table[mask]
            mask = self.current_models.table['which_iter'] == last_iter-1
            models1 = self.current_models.table[mask]
            chi2 = 'chi2' if self.parspace_settings['which_chi2'] == 'chi2' \
                else 'kinchi2'
            # Don't use abs() so we catch increasing chi2 values, too:
            delta_chi2 = np.min(models1[chi2]) - np.min(models0[chi2])
            if delta_chi2 <= \
                self.parspace_settings['stopping_criteria']['min_delta_chi2']:
                self.status['min_delta_chi2_reached'] = True
        # (ii) if step_size < min_step_size for all params
        #       => dealt with by grid_walk (doesn't create such models)


class GaussianProcessEmulator(ParameterGenerator):

    def generate(self, current_models=None, n_new=0):
        # actual code to do gaussian process emulation
        # return new_parameter_list of length n_new
        return []

    def check_specific_stopping_critera(self, current_models):
        stop = True # or false
        return stop


# end
