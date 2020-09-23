import model
import parameter_space
import numpy as np


class ModelIterator(object):

    def __init__(self,
                 system=None,
                 all_models=None,
                 settings=None,
                 model_kwargs={},
                 executor=None):
        stopping_crit = settings.parameter_space_settings['stopping_criteria']
        n_max_iter = stopping_crit['n_max_iter']
        # get specified parameter generator
        parspace = parameter_space.ParameterSpace(system)
        par_generator_type = settings.parameter_space_settings['generator_type']
        kwargs = {'parspace_settings':settings.parameter_space_settings}
        par_generator = getattr(parameter_space, par_generator_type)(parspace,
                                                                     **kwargs)
        model_inner_iterator = ModelInnerIterator(
            system=system,
            all_models=all_models,
            settings=settings,
            par_generator=par_generator,
            model_kwargs=model_kwargs)
        for iter in range(n_max_iter):
            print(f'{par_generator_type}: "iteration {iter}"')
            status = model_inner_iterator.run_iteration(iter,
                                                        executor=executor)
            if status['stop'] is True:
                print(f'Stopping after iteration {iter}')
                print(status)
                break


class ModelInnerIterator(object):

    def __init__(self,
                 iter=0,
                 system=None,
                 all_models=None,
                 settings=None,
                 par_generator=None,
                 model_kwargs={}):
        self.system = system
        self.all_models = all_models
        self.settings = settings
        self.parspace = parameter_space.ParameterSpace(system)
        self.par_generator = par_generator
        self.model_kwargs = model_kwargs

    def run_iteration(self, iter, executor=None):
        self.par_generator.generate(current_models=self.all_models)
        # generate parameter sets for this iteration
        if self.par_generator.status['stop'] is False:
            # find models not yet done
            rows_to_do = np.where(self.all_models.table['all_done'] == False)
            rows_to_do = rows_to_do[0]
            n_to_do = len(rows_to_do)
            for i, row in enumerate(rows_to_do):
                print(f'... running model {i+1} out of {n_to_do}')
                # extract the parameter values
                parset0 = self.all_models.table[row]
                parset0 = parset0[self.parspace.par_names]
                # create and run the model
                mod0 = self.create_model(parset0,
                                         model_kwargs=self.model_kwargs,
                                         executor=executor)
                mod0.run()
                self.all_models.table['chi2'][row] = mod0.chi2
                self.all_models.table['kinchi2'][row] = mod0.kinchi2
                # self.all_models.table['which_iter'][row] = iter
                self.all_models.table['all_done'][row] = True
                time_now = np.datetime64('now', 'ms')
                self.all_models.table['time_modified'][row] = time_now
        return self.par_generator.status

    def create_model(self,
                     parset,
                     model_kwargs={},
                     executor=None):
        model_kwargs0 = {'system':self.system,
                         'settings':self.settings,
                         'parspace':self.parspace,
                         'parset':parset,
                         'executor':executor}
        model_kwargs0.update(model_kwargs)
        # create a model object based on choices in settings
        if self.settings.legacy_settings['use_legacy_mode']:
            mod = getattr(model, 'LegacySchwarzschildModel')(**model_kwargs0)
        else:
            # TODO: create other model classes based on a choice of:
            # (i) orbit library generator
            # (i) weight solver
            # (iii) colour solver
            # mod = getattr(model, '...')(**model_kwargs0)
            raise ValueError("""
                             Only Legacy Mode currently implemented. Set
                                 legacy_settings:
                                     use_legacy_mode: True
                             in the config file
                             """)
        return mod







# end
