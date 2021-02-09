# some tricks to add the current path to sys.path (so the imports below work)

import os.path
import sys

this_dir = os.path.dirname(__file__)
if not this_dir in sys.path:
    sys.path.append(this_dir)

import numpy as np
import copy
import parameter_space as parspace
from astropy.table import Table
import logging

class Parameter(object):

    attributes = []
    def __init__(self,
                 name=None,
                 fixed=False,
                 LaTeX=None,
                 sformat="%g",
                 value=None,
                 par_generator_settings=None,
                 gpe_parspace_settings=None,
                 logarithmic=False,
                 ):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.name = name
        self.fixed = fixed
        self.LaTeX = LaTeX
        self.sformat = sformat
        self.value = value
        self.par_generator_settings = par_generator_settings
        self.gpe_parspace_settings = gpe_parspace_settings
        self.logarithmic = logarithmic
        self.__class__.attributes = list(self.__dict__.keys())

    def update(self, **kwargs):
        for k, v in kwargs.items():
            if k not in self.__class__.attributes:
                text = (f'Invalid parameter key {k}. Allowed keys: '
                        f'{str(tuple(self.__class__.attributes))}')
                self.logger.error(text)
                raise ValueError(text)
            setattr(self, k, v)

    def validate(self):
        if sorted(self.__class__.attributes) != sorted(self.__dict__.keys()):
            text = (f'Parameter attributes can only be '
                    f'{str(tuple(self.__class__.attributes))}, '
                    f'not {str(tuple(self.__dict__.keys()))}')
            self.logger.error(text)
            raise ValueError(text)

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
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        for cmp in system.cmp_list:
            for par in cmp.parameters:
                self.append(par)
                self.system = system
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

    def get_parameter_from_name(self, name):
        name_array = np.array(self.par_names)
        idx = np.where(name_array == name)
        self.logger.debug(f'Checking unique parameter name {name}...')
        error_msg = f"There should be 1 and only 1 parameter named {name}"
        assert len(idx[0]) == 1, error_msg
        self.logger.debug('...check ok.')
        parameter = self[idx[0][0]]
        return parameter

    def get_parset(self):
        t = Table()
        for par in self:
            raw_value = par.value
            par_value = par.get_par_value_from_raw_value(raw_value)
            t[par.name] = [par_value]
        # extract 0th - i.e. the only - row from the table
        parset = t[0]
        return parset

    def validate_parset(self, parset):
        """
        Validates the values of each component's parameters by calling the
        individual components' validate_parameter methods. Does the same for
        system parameters.

        Parameters
        ----------
        parset : list of Parameter objects

        Returns
        -------
        bool
            True if validation was successful, False otherwise

        """
        isvalid = True
        for comp in self.system.cmp_list:
            par = {p.name[:p.name.rfind(f'_{comp.name}')]:p.value \
                   for p in parset if p.name.rfind(f'_{comp.name}')>=0}
            isvalid = isvalid and comp.validate_parset(par)
        par = {p.name:p.value for p in parset \
               if p.name in [n.name for n in self.system.parameters]}
        isvalid = isvalid and self.system.validate_parset(par)
        return isvalid

class ParameterGenerator(object):

    def __init__(self,
                 par_space=[],
                 parspace_settings=None,
                 name=None):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.par_space = par_space
        if not parspace_settings:
            text = 'ParameterGenerator needs parspace_settings'
            self.logger.error(text)
            raise ValueError(text)
        self.parspace_settings = parspace_settings
        which_chi2 = self.parspace_settings.get('which_chi2')
        if which_chi2 not in ['chi2', 'kinchi2']:
          text = 'Unknown or missing which_chi2 setting, use chi2 or kinchi2'
          self.logger.error(text)
          raise ValueError(text)
        self.chi2 = 'chi2' if which_chi2 == 'chi2' else 'kinchi2'
        self.status = {}
        self.name = name
        self.lo = []
        self.hi = []
        try:
            for par in self.par_space:
                settings = par.par_generator_settings
                if par.fixed is False:
                    self.lo.append(settings['lo'])
                    self.hi.append(settings['hi'])
                else:
                    self.lo.append([])
                    self.hi.append([])
        except:
            text = 'ParameterGenerator: non-fixed parameters ' + \
                   'need hi and lo settings'
            self.logger.error(text)
            raise ValueError(text)
        try:
            stop_crit = parspace_settings['stopping_criteria']
        except:
            text = 'ParameterGenerator: need stopping criteria'
            self.logger.error(text)
            raise ValueError(text)
        if not stop_crit.get('n_max_mods') and \
           not stop_crit.get('n_max_iter'):
            text = 'ParameterGenerator: need n_max_mods and ' + \
                   'n_max_iter stopping criteria settings'
            self.logger.error(text)
            raise ValueError(text)

    def generate(self,
                 current_models=None,
                 kw_specific_generate_method={}):
        """Generate new parameter sets. This is a wrapper method around the
        specific_generate_method of child generator classes. This wrapper does
        the following:
        (i) evaluates stopping criteria, and stop if necessary
        (ii) runs the specific_generate_method of the child class, which
        updates self.model_list with a list of propsal models
        (iii) removes previously run and/or invalid models from self.model_list
        (iv) converts parameters from raw_values to par_values
        (v) adds new, valid models to current_models.table
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
            self.logger.error(errormsg)
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
        self.logger.info(f'{self.name} added {newmodels} new model(s) out of '
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
            self.logger.error('No or empty model')
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
        Checks whether model has valid parameter values and it is a new model
        (i.e., its parameter set does not exist in self.current_models).

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
            self.logger.error('Model arg. must be list of Parameter objects')
            raise ValueError('Model arg. must be list of Parameter objects')
        if not self.par_space.validate_parset(model):
            isnew = False
        else:
            isnew = True
            raw_model_values = [p.value for p in model]
            model_values = \
                self.par_space.get_param_value_from_raw_value(raw_model_values)
            if len(self.current_models.table) > 0:
                for mod in self.current_models.table[self.par_space.par_names]:
                    if np.allclose(list(mod), model_values, rtol=eps):
                        isnew = False
                        break
        return isnew


class LegacyGridSearch(ParameterGenerator):

    def __init__(self, par_space=[], parspace_settings=None):
        super().__init__(par_space=par_space,
                         parspace_settings=parspace_settings,
                         name='LegacyGridSearch')
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        # We need a local parameter copy because we don't want to change the
        # minstep in the original par_space:
        self.new_parset = [copy.deepcopy(p) for p in self.par_space]
        self.step = []
        self.minstep = []
        try:
            for par in self.par_space:
                settings = par.par_generator_settings
                if par.fixed is False:
                    self.step.append(settings['step'])
                    # Use 'minstep' value if present, otherwise use 'step'.
                    # Explicitly set minstep=0 to allow arbitrarily
                    # small steps, not recommended.
                    self.minstep.append(settings['minstep'] \
                        if 'minstep' in settings else self.step)
                else:
                    self.step.append(None)
                    self.minstep.append(None)
        except:
            text = 'LegacyGridSearch: non-fixed parameters need step setting'
            self.logger.error(text)
            raise ValueError(text)
        try:
            self.thresh = \
            self.parspace_settings['generator_settings']['threshold_del_chi2']
        except:
            text = 'LegacyGridSearch: need generator_settings - ' + \
                'threshold_del_chi2 (absolute or scaled - see documentation)'
            self.logger.error(text)
            raise ValueError(text)
        stop_crit = parspace_settings['stopping_criteria']
        self.min_delta_chi2_abs = stop_crit.get('min_delta_chi2_abs', False)
        self.min_delta_chi2_rel = stop_crit.get('min_delta_chi2_rel', False)
        if (not self.min_delta_chi2_abs and not self.min_delta_chi2_rel) \
           or \
           (self.min_delta_chi2_abs and self.min_delta_chi2_rel):
            text = 'LegacyGridSearch: specify exactly one of the ' + \
                   'options min_delta_chi2_abs, min_delta_chi2_rel'
            self.logger.error(text)
            raise ValueError(text)

    def specific_generate_method(self, **kwargs):
        """
        Generates list of new models self.model_list. Each element of
        self.model_list is a list of Parameter objects. The grid search is the
        driven by the parameter set with the smallest chi2 value min_chi2.
        Note that the parameter space setting 'which_chi2' determines whether
        chi2 or kinchi2 is used. All models with abs(chi2-min_chi2)<=
        threshold_del_chi2 (given in the config file) are subject to the
        search. New models are generated based on varying one non-fixed
        parameter at a time. Parameter lo and hi attributes as well as
        minstep are considered when creating new models. If no new models can
        be found, the stepsize of all parameters is halved (up to minstepsize).

        Parameters
        ----------
        None.

        Returns
        -------
        None. self.model_list is the list of new models.

        """
        if len(self.current_models.table) == 0:
            # The 'zeroth iteration' results in only one model
            # (all parameters at their .value level)
            self.model_list = [[p for p in self.par_space]]
            return ###########################################################
        min_chi2 = np.min(self.current_models.table[self.chi2])
        prop_mask = \
            abs(self.current_models.table[self.chi2]-min_chi2) <= self.thresh
        prop_list = self.current_models.table[prop_mask]
        self.model_list = []
        step_ok = True
        while step_ok and len(self.model_list) == 0:
            for paridx in range(len(self.new_parset)):
                par = self.new_parset[paridx]
                if par.fixed: # parameter fixed -> do nothing
                    continue
                lo = self.lo[paridx] #par.par_generator_settings['lo']
                hi = self.hi[paridx] #par.par_generator_settings['hi']
                step = self.step[paridx] #par.par_generator_settings['step']
                minstep = self.minstep[paridx]
                for m in prop_list: # for all models within threshold_del_chi2
                    for p in self.new_parset:
                        p.value = p.get_raw_value_from_par_value(m[p.name])
                    center_value = self.new_parset[paridx].value
                    for s in [-1, 1]:
                        new_value = np.clip(center_value + s*step, lo, hi)
                        if abs(new_value-par.value) >= minstep:
                            self.new_parset[paridx].value = new_value
                            if self._is_newmodel(self.new_parset, eps=1e-10):
                                self.model_list.append\
                                  ([copy.deepcopy(p) for p in self.new_parset])
#                                    (copy.deepcopy(self.new_parset))
            # If no new models: cut stepsize in half & try again
            if len(self.model_list) == 0:
                step_ok = False
                for par in [p for p in self.new_parset if not p.fixed]:
                    paridx = self.new_parset.index(par)
                    minstep = self.minstep[paridx]
                    if self.step[paridx]/2 >= minstep:
                        self.step[paridx] /= 2
                        # the following line is just to record the step size
                        # in self.new_parset and can be commented out...
                        par.par_generator_settings['step'] = self.step[paridx]
                        step_ok = True
        return

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
            # Don't use abs() so we stop on increasing chi2 values, too:
            delta_chi2 = np.min(models1[self.chi2])-np.min(models0[self.chi2])
            if self.min_delta_chi2_rel:
                delta_chi2 /= np.min(models1[self.chi2])
                delta_chi2 /= self.min_delta_chi2_rel
            else:
                delta_chi2 /= self.min_delta_chi2_abs
            if delta_chi2 <= 1:
                self.status['min_delta_chi2_reached'] = True
        # (ii) if step_size < min_step_size for all params
        #       => dealt with by grid_walk (doesn't create such models)


class GridWalk(ParameterGenerator):

    def __init__(self,
                 par_space=[],
                 parspace_settings=None):
        super().__init__(par_space=par_space,
                         parspace_settings=parspace_settings,
                         name='GridWalk')
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.step = []
        self.minstep = []
        try:
            for par in self.par_space:
                settings = par.par_generator_settings
                if par.fixed is False:
                    self.step.append(settings['step'])
                    # use 'minstep' value if present, otherwise use 'step'
                    self.minstep.append(settings['minstep'] \
                        if 'minstep' in settings else self.step)
                else:
                    self.step.append([])
                    self.minstep.append([])
        except:
            text = 'GridWalk: non-fixed parameters need step setting'
            self.logger.error(text)
            raise ValueError(text)

        stop_crit = parspace_settings['stopping_criteria']
        self.min_delta_chi2_abs = stop_crit.get('min_delta_chi2_abs', False)
        self.min_delta_chi2_rel = stop_crit.get('min_delta_chi2_rel', False)
        if (not self.min_delta_chi2_abs and not self.min_delta_chi2_rel) \
           or \
           (self.min_delta_chi2_abs and self.min_delta_chi2_rel):
            text = 'GridWalk: specify exactly one of the ' + \
                   'options min_delta_chi2_abs, min_delta_chi2_rel'
            self.logger.error(text)
            raise ValueError(text)

    def specific_generate_method(self, **kwargs):
        """
        Generates list of new models self.model_list. Each element of
        self.model_list is a list of Parameter objects. The center of the
        grid walk is the parameter set with the smallest chi2 value. Note
        that the parameter space setting 'which_chi2' determines whether chi2
        or kinchi2 is used.

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
            # center criterion: min(chi2)
            center_idx = np.argmin(self.current_models.table[self.chi2])
            n_par = self.par_space.n_par
            center = list(self.current_models.table[center_idx])[:n_par]
            raw_center = self.par_space.get_raw_value_from_param_value(center)
            self.logger.debug(f'center: {center}')
            # Build model_list by walking the grid
            self.model_list = []
            self.grid_walk(center=raw_center)
            # for m in self.model_list:
            #     self.logger.debug(f'{[(p.name, p.value) for p in m]}')
        return

    def grid_walk(self, center=None, par=None, eps=1e-6):
        """
        Walks the grid defined by self.par_space.par_generator_settings
        attributes.
        Clips parameter values to lo/hi attributes. If clipping violates the
        minstep attribute, the resulting model(s) will not be created. If the
        minstep attribute is missing, the step attribute will be used instead.
        Explicitly set minstep=0 to allow arbitrarily small steps (not
        recommended).

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
        if center is None:
            text = 'Need center'
            self.logger.error(text)
            raise ValueError(text)
        if not par:
            par = self.par_space[0]
        paridx = self.par_space.index(par)
        self.logger.debug(f'Call with paridx={paridx}, '
                          f'n_par={self.par_space.n_par}')

        if par.fixed:
            par_values = [par.value]
            if abs(center[paridx] - par.value) > eps:
                text='Something is wrong: fixed parameter value not in center'
                self.logger.error(text)
                raise ValueError(text)
        else:
            lo = self.lo[paridx]
            hi = self.hi[paridx]
            step = self.step[paridx]
            minstep = self.minstep[paridx]
            # up to 3 *distinct* par_values (clipped lo, mid, hi values)
            par_values = []
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
                self.logger.debug('new model list, starting w/parameter '
                                  f'{parcpy.name}')
            elif parcpy.name in [p.name for p in self.model_list[0]]:
                # in this case, create new (partial) model by copying last
                # models and setting the new parameter value
                for m in models_prev:
                    new_model = m + [parcpy]
                    self.model_list.append(new_model)
                self.logger.debug(f'{parcpy.name} is in '
                      f'{[p.name for p in self.model_list[0]]}, '
                      f'added {parcpy.name}={parcpy.value}')
            else: # new parameter: append it to existing (partial) models
                models_prev = copy.deepcopy(self.model_list)
                for m in self.model_list:
                    m.append(parcpy)
                self.logger.debug(f'new parameter {parcpy.name}={parcpy.value}')

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
        logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if mini <= maxi:
            return min(max(mini, value), maxi)
        else:
            text = 'Clip error: minimum must be less than or equal to maximum'
            logger.error(text)
            raise ValueError(text)

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
            # Don't use abs() so we stop on increasing chi2 values, too:
            delta_chi2 = np.min(models1[self.chi2])-np.min(models0[self.chi2])
            if self.min_delta_chi2_rel:
                delta_chi2 /= np.min(models1[self.chi2])
                delta_chi2 /= self.min_delta_chi2_rel
            else:
                delta_chi2 /= self.min_delta_chi2_abs
            if delta_chi2 <= 1:
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
