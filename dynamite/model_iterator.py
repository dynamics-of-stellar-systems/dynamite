import model
import parameter_space
import numpy as np


class ModelIterator(object):

    def __init__(self,
                 system=None,
                 all_models=None,
                 settings=None,
                 model_kwargs={},
                 executor=None,
                 do_dummy_run=None,
                 dummy_chi2_function=None):
        stopping_crit = settings.parameter_space_settings['stopping_criteria']
        n_max_iter = stopping_crit['n_max_iter']
        self.n_max_mods = stopping_crit['n_max_mods']
        self.all_models = all_models
        # get specified parameter generator
        parspace = parameter_space.ParameterSpace(system)
        par_generator_type = settings.parameter_space_settings['generator_type']
        kwargs = {'parspace_settings':settings.parameter_space_settings}
        par_generator = getattr(parameter_space, par_generator_type)(parspace,
                                                                     **kwargs)
        model_inner_iterator = ModelInnerIterator(
            system=system,
            all_models=self.all_models,
            settings=settings,
            par_generator=par_generator,
            do_dummy_run=do_dummy_run,
            dummy_chi2_function=dummy_chi2_function)
        if len(self.all_models.table)>0:
            previous_iter = np.max(self.all_models.table['which_iter'])+1
        else:
            previous_iter = 0
        status = {}
        status['stop'] = False
        for iter in range(n_max_iter):
            total_iter_count = previous_iter + iter
            n_models_done = np.sum(self.all_models.table['all_done'])
            if n_models_done>=self.n_max_mods:
                status['n_max_mods_reached'] = True
                status['stop'] = True
            if status['stop'] is True:
                print(f'Stopping after iteration {total_iter_count}')
                print(status)
                break
            print(f'{par_generator_type}: "iteration {total_iter_count}"')
            status = model_inner_iterator.run_iteration(iter, executor=executor)


class ModelInnerIterator(object):

    def __init__(self,
                 iter=0,
                 system=None,
                 all_models=None,
                 settings=None,
                 par_generator=None,
                 do_dummy_run=False,
                 dummy_chi2_function=None):
        self.system = system
        self.all_models = all_models
        self.settings = settings
        self.parspace = parameter_space.ParameterSpace(system)
        self.par_generator = par_generator
        self.do_dummy_run = do_dummy_run
        if self.do_dummy_run:
            assert dummy_chi2_function is not None
            # TODO: assert dummy_chi2_function is a valid function of parset
        self.dummy_chi2_function = dummy_chi2_function

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
                mod0 = self.create_model(parset0, executor=executor)
                if self.do_dummy_run:
                    mod0.chi2 = self.dummy_chi2_function(parset0)
                    mod0.kinchi2 = 0.
                else:
                    mod0.setup_directories()
                    mod0.get_orblib()
                    self.all_models.table['orblib_done'][row] = True
                    mod0.get_weights()
                    self.all_models.table['weights_done'][row] = True
                # store results
                self.all_models.table['chi2'][row] = mod0.chi2
                self.all_models.table['kinchi2'][row] = mod0.kinchi2
                # 'which_iter' column is filled by ParameterGenerator.add_model
                self.all_models.table['all_done'][row] = True
                time_now = np.datetime64('now', 'ms')
                self.all_models.table['time_modified'][row] = time_now
        return self.par_generator.status

    def create_model(self,
                     parset,
                     executor=None):
        model_kwargs = {'system':self.system,
                        'settings':self.settings,
                        'parspace':self.parspace,
                        'parset':parset,
                        'executor':executor}
        # create a model object based on choices in settings
        if self.settings.legacy_settings['use_legacy_mode']:
            mod = getattr(model, 'LegacySchwarzschildModel')(**model_kwargs)
        else:
            # TODO: create other model classes based on a choice of:
            # (i) orbit library generator
            # (i) weight solver
            # (iii) colour solver
            # mod = getattr(model, '...')(**model_kwargs)
            raise ValueError("""
                             Only Legacy Mode currently implemented. Set
                                 legacy_settings:
                                     use_legacy_mode: True
                             in the config file
                             """)
        return mod







# end
