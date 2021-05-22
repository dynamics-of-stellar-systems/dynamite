import os
import numpy as np
import logging
from pathos.multiprocessing import Pool
import matplotlib.pyplot as plt

from dynamite import model
from dynamite import parameter_space
from dynamite import plotter

class ModelIterator(object):

    def __init__(self,
                 system=None,
                 all_models=None,
                 settings=None,
                 model_kwargs={},
                 do_dummy_run=None,
                 dummy_chi2_function=None,
                 ncpus=1,
                 plots=True):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')

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

        if plots:
            self.the_plotter = plotter.Plotter(system = system,
                                               settings = settings,
                                               parspace = parspace,
                                               all_models = all_models)

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
        for iteration in range(n_max_iter):
            total_iter_count = previous_iter + iteration
            n_models_done = np.sum(self.all_models.table['all_done'])
            if n_models_done>=self.n_max_mods:
                status['n_max_mods_reached'] = True
                status['stop'] = True
            if status['stop'] is True:
                self.logger.info(f'Stopping after iteration {total_iter_count}')
                self.logger.debug(status)
                break
            self.logger.info(f'{par_generator_type}: "iteration '
                        f'{total_iter_count}"')
            status = model_inner_iterator.run_iteration()
            if plots:
                self.the_plotter.make_chi2_vs_model_id_plot()
                self.the_plotter.make_chi2_plot()
                self.the_plotter.plot_kinematic_maps(kin_set='all',
                                                     cbar_lims='data')
                plt.close('all') # just to make sure...


class ModelInnerIterator(object):

    def __init__(self,
                 system=None,
                 all_models=None,
                 settings=None,
                 par_generator=None,
                 do_dummy_run=False,
                 dummy_chi2_function=None,
                 ncpus=1):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.system = system
        self.all_models = all_models
        self.settings = settings
        self.parspace = parameter_space.ParameterSpace(system)
        self.orblib_parameters = self.parspace.par_names[:]
        ml = 'ml'
        try:
            self.orblib_parameters.remove(ml)
        except:
            self.logger.error(f"Parameter '{ml}' not found - check "
                              "implementation")
            raise
        self.logger.debug(f'orblib_parameters: {self.orblib_parameters}')
        self.par_generator = par_generator
        self.do_dummy_run = do_dummy_run
        if self.do_dummy_run:
            assert dummy_chi2_function is not None
            # TODO: assert dummy_chi2_function is a valid function of parset
        self.dummy_chi2_function = dummy_chi2_function
        self.ncpus = ncpus

    def run_iteration(self):
        """
        Executes one iteration step: run all models in self.all_models.table
        for which all_done == False. The model runs (1) build the orbit
        library and (2) execute weight_solver. The models are run in parallel
        threads as defined in the ncpus parameter in the configuration file.
        In case multiple models comprise the same (new) orbit library it is
        ensured that they are calculated only once, avoiding conflicting
        threads.

        Returns
        -------
        dict
            ParameterGenerator.status.

        """

        self.par_generator.generate(current_models=self.all_models)
        self.all_models.save() # save all_models table once parameters are added
        # generate parameter sets for this iteration
        if self.par_generator.status['stop'] is False:
            # find models not yet done
            rows_to_do = np.where(self.all_models.table['all_done'] == False)
            rows_to_do = rows_to_do[0]
            self.logger.debug(f'rows_to_do: {rows_to_do}.')
            self.n_to_do = len(rows_to_do)
            # rows_to_do_orblib are the rows that need orblib and weight_solver
            rows_to_do_orblib=[i for i in rows_to_do if self.is_new_orblib(i)]
            n_orblib = len(rows_to_do_orblib)
            # rows_to_do_ml are the rows that need weight_solver only
            rows_to_do_ml=[i for i in rows_to_do if i not in rows_to_do_orblib]
            input_list_orblib = [i+(True,) for i in enumerate(rows_to_do_orblib)]
            input_list_ml=[i+(False,) for i in enumerate(rows_to_do_ml, start=n_orblib)]
            self.logger.debug(f'input_list_orblib: {input_list_orblib}, '
                              f'input_list_ml: {input_list_ml}.')
            # input_list = []
            # for i, row in enumerate(rows_to_do):
            #     input_list += [(i, row)]
            # self.logger.debug(f'input_list: {input_list}')
            with Pool(self.ncpus) as p:
                # output = p.map(self.create_and_run_model, input_list)
                output_orblib = \
                    p.map(self.create_and_run_model, input_list_orblib)
                output_ml = p.map(self.create_and_run_model, input_list_ml)
            # save the output
            # self.write_output_to_all_models_table(rows_to_do, output)
            self.write_output_to_all_models_table(rows_to_do_orblib,
                                                  output_orblib)
            self.write_output_to_all_models_table(rows_to_do_ml, output_ml)
            self.all_models.save() # save all_models table once models are run
        return self.par_generator.status

    def is_new_orblib(self, row_idx):
        """
        Checks whether the orbit library characterized by the parameters in
        row number row_idx exists in earlier rows of self.all_models.table.

        Parameters
        ----------
        row_idx : int
            Row index of the model entry to be checked.

        Returns
        -------
        is_new : bool
            True if no earlier row contains the orbit library, False otherwise.

        """
        all_data = self.all_models.table[self.orblib_parameters]
        row_data = all_data[row_idx]
        previous_data = all_data[:row_idx]
        if any(np.allclose(tuple(row_data), tuple(r)) for r in previous_data):
            self.logger.debug('Orblib exists above in table: '
                              f'{row_data}.')
            is_new = False
        else:
            self.logger.debug(f'New orblib: {row_data}.')
            is_new = True
        return is_new

    def create_and_run_model(self, input):
        i, row, new_orblib = input
        self.logger.info(f'... running model {i+1} out of {self.n_to_do}')
        # extract the parameter values
        parset0 = self.all_models.table[row]
        parset0 = parset0[self.parspace.par_names]
        # create and run the model
        mod0 = model.Model(system=self.system,
                           settings=self.settings,
                           parspace=self.parspace,
                           parset=parset0)
        orb_done = False
        wts_done = False
        if self.do_dummy_run:
            mod0.chi2 = self.dummy_chi2_function(parset0)
            mod0.kinchi2 = 0.
        else:
            # If new orblib: check for duplicate directory conflict
            if new_orblib:
                orblib_dir = mod0.get_model_directory()
                orblib_dir = orblib_dir[:orblib_dir.rindex('/', 0, -1)]
                orblib_dir = orblib_dir[orblib_dir.rindex('/')+1:]
                self.logger.debug(f'orblib_dir: {orblib_dir}')
                model_dir = \
                    self.settings.io_settings['output_directory'] + 'models/'
                if os.path.isdir(model_dir):
                    _, orblib_dirs, _ = next(os.walk(model_dir))
                    self.logger.debug(f'orblib_dirs: {orblib_dirs}')
                    if orblib_dir in orblib_dirs:
                        t = 'Cannot create orbit library directory ' \
                            f'{model_dir}{orblib_dir} because it already ' \
                            'exists. Caused by model with parameter set ' \
                            f'{parset0}. ' \
                            'Hint: check the parameter values, their ' \
                            'stepsize, and the parameter string formats in ' \
                            'the Component classes in physical_system.py.'
                        self.logger.error(t)
                        raise RuntimeError(t)
                        # mod0.chi2 = float('nan')
                        # mod0.kinchi2 = float('nan')
                        # all_done = orb_done and wts_done
                        # time = np.datetime64('now', 'ms')
                        # output = (orb_done, wts_done, mod0.chi2,
                        #           mod0.kinchi2, all_done, time)
                        # return output
                else:
                    self.logger.debug('...is a new orblib directory.')
            # Carry on, there is no orblib directory conflict...
            mod0.setup_directories()
            orblib = mod0.get_orblib()
            orb_done = True
            weight_solver = mod0.get_weights(orblib)
            wts_done = True
        all_done = orb_done and wts_done
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
