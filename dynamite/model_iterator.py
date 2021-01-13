import model
import parameter_space
import numpy as np
import pathos
from pathos.multiprocessing import Pool

class ModelIterator(object):

    def __init__(self,
                 system=None,
                 all_models=None,
                 settings=None,
                 model_kwargs={},
                 do_dummy_run=None,
                 dummy_chi2_function=None,
                 ncpus=1):
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
            dummy_chi2_function=dummy_chi2_function,
            ncpus=ncpus)
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
                print(f'Saving all_models table')
                break
            print(f'{par_generator_type}: "iteration {total_iter_count}"')
            status = model_inner_iterator.run_iteration(iter)
        self.all_models.save()


class ModelInnerIterator(object):

    def __init__(self,
                 iter=0,
                 system=None,
                 all_models=None,
                 settings=None,
                 par_generator=None,
                 do_dummy_run=False,
                 dummy_chi2_function=None,
                 ncpus=1):
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
        self.ncpus = ncpus

    def run_iteration(self, iter):
        self.par_generator.generate(current_models=self.all_models)
        # generate parameter sets for this iteration
        if self.par_generator.status['stop'] is False:
            # find models not yet done
            rows_to_do = np.where(self.all_models.table['all_done'] == False)
            rows_to_do = rows_to_do[0]
            self.n_to_do = len(rows_to_do)
            input_list = []
            for i, row in enumerate(rows_to_do):
                input_list += [(i, row)]
            with Pool(self.ncpus) as p:
                output = p.map(self.create_and_run_model, input_list)
            # save the output
            self.write_output_to_all_models_table(rows_to_do, output)
        return self.par_generator.status

    def create_model(self,
                     parset):
        model_kwargs = {'system':self.system,
                        'settings':self.settings,
                        'parspace':self.parspace,
                        'parset':parset}
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

    def create_and_run_model(self, input):
        i, row = input
        print(f'... running model {i+1} out of {self.n_to_do}')
        # extract the parameter values
        parset0 = self.all_models.table[row]
        parset0 = parset0[self.parspace.par_names]
        # create and run the model
        mod0 = self.create_model(parset0)
        if self.do_dummy_run:
            mod0.chi2 = self.dummy_chi2_function(parset0)
            mod0.kinchi2 = 0.
        else:
            mod0.setup_directories()
            orblib = mod0.get_orblib()
            orb_done = True
            mod0.get_weights(orblib)
            wts_done = True
        all_done = True
        time = np.datetime64('now', 'ms')
        output = orb_done, wts_done, mod0.chi2, mod0.kinchi2, all_done, time
        return output

    def write_output_to_all_models_table(self, rows_to_do, output):
        for i, row in enumerate(rows_to_do):
            orb_done, wts_done, chi2, kinchi2, all_done, time = output[i]
            self.all_models.table['orblib_done'][row] = orb_done
            self.all_models.table['weights_done'][row] = wts_done
            self.all_models.table['chi2'][row] = chi2
            self.all_models.table['kinchi2'][row] = kinchi2
            self.all_models.table['all_done'][row] = all_done
            self.all_models.table['time_modified'][row] = time




# end
