import os
import numpy as np
import logging
from pathos.multiprocessing import Pool
import matplotlib.pyplot as plt

from dynamite import model
from dynamite import parameter_space
from dynamite import plotter

class ModelIterator(object):
    """Iterator for models

    Creating this ``ModelIterator`` object will (i) generate parameters sets,
    (ii) run models for those parameters, (iii) check stopping criteria, and
    iterate this procedure till a stopping criterion is met. This is implemented
    by created a ``ModelInnerIterator`` object whose ``run_iteration`` method is
    called a number of times.

    Parameters
    ----------
    system : a ``dyn.physical_system.System`` object
    all_models : a ``dyn.model.AllModels`` object
    settings : a ``dyn.config_reader.Settings`` object
    model_kwargs : dict
        other kewyord argument required for this model
    do_dummy_run : Bool
        whether this is a dummy run - if so, dummy_chi2_funciton is executed
        instead of the model (for testing!)
    dummy_chi2_function : function
        a function of model parameters to be executed instead of the real model
    ncpus : int
        number of cpus for multiprocessing
    plots : bool
        whether or not to make plots

    """
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
                self.logger.info(f'Stopping at iteration {total_iter_count}')
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
    """Class to run all models in a single iteration.

    Uses ``pathos.multiprocessing.Pool`` to execute the models

    Parameters
    ----------
    system : a DYNMAITE system object
    all_models : a DYNMAITE all_models object
    settings : a DYNMAITE settings object
    model_kwargs : type
        Description of parameter `model_kwargs`.
    do_dummy_run : Bool
        whether this is a dummy run - if so, dummy_chi2_funciton is executed
        instead of the model (for testing!)
    dummy_chi2_function : function
        a function of model parameters to be executed instead of the real model
    ncpus : int
        number of cpus for multiprocessing

    """
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
        """run one iteration step

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
            input_list_orblib = [i for i in enumerate(rows_to_do_orblib)]
            input_list_ml=[i for i in enumerate(rows_to_do_ml, start=n_orblib)]
            self.logger.debug(f'input_list_orblib: {input_list_orblib}, '
                              f'input_list_ml: {input_list_ml}.')
            self.assign_model_directories(rows_to_do_orblib, rows_to_do_ml)

            with Pool(self.ncpus) as p:
                output_orblib = \
                    p.map(self.create_and_run_model, input_list_orblib)
                output_ml = p.map(self.create_and_run_model, input_list_ml)
            # save the output
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

    def assign_model_directories(self, rows_orblib=None, rows_ml=None):
        """
        Assigns model directories in all_models.table.

        Models indexed by rows_orblib:
        The model directories follow the pattern orblib_xxx_yyy/mlz.zz where
        xxx is the iteration number, yyy a consecutive number of that
        iteration's orbit library, and z.zz is the value of the models'
        ml parameter in the 01.2f format (the sformat set in the System class).

        Models indexed by rows_ml:
        These models re-use an existing orbit library. Hence, their directory
        strings re-use an existing orblib_xxx_yyy part and get augmented with
        the appropriate /mlz.zz.

        Parameters
        ----------
        rows_orblib : list, optional
            Indices of models with new orbit libraries. The default is None.
        rows_ml : list, optional
            Indices of models with existing orbit libraries.
            The default is None.

        Raises
        ------
        ValueError
            If the orbit library of a model in rows_ml cannot be found in
            all_models.table.

        Returns
        -------
        None.

        """
        iteration = self.all_models.table['which_iter'][-1]
        # new orblib model directories
        for row in rows_orblib:
            n=np.sum(self.all_models.table[:row]['which_iter']==iteration)
            orblib_dir = f'orblib_{iteration:03d}_{n:03d}'
            self.all_models.table[row]['directory'] = orblib_dir
        # existing orblib directories
        orblib_data = self.all_models.table[self.orblib_parameters]
        for row in rows_ml:
            row_data = orblib_data[row]
            for idx, orblib in enumerate(orblib_data[:row]):
                if np.allclose(tuple(row_data), tuple(orblib)):
                    orblib_dir = self.all_models.table[idx]['directory']
                    orblib_dir = orblib_dir[:orblib_dir[:-1].rindex('/')]
                    break
            else:
                text = f'Unexpected: cannot find orblib {dict(row_data)}.'
                self.logger.error(text)
                raise ValueError(text)
            self.all_models.table[row]['directory'] = orblib_dir
        # ml directories
        sformat = self.system.parameters[0].sformat # this is ml's format
        for row in rows_orblib+rows_ml:
            ml_dir = f"/ml{self.all_models.table['ml'][row]:{sformat}}/"
            self.all_models.table[row]['directory'] += ml_dir
            self.logger.debug(f"New model directory "
                f"{self.all_models.table[row]['directory']} assigned.")

    def create_and_run_model(self, which_model):
        """main method to create and run a model

        Parameters
        ----------
        input : tuple
            (i, row) where i is the index of a model in this iteration, and row
            is the row index of the all_models table for this model

        Returns
        -------
        tuple
            all the output for this model, bundles up in a tuple

        """
        i, row = which_model
        self.logger.info(f'... running model {i+1} out of {self.n_to_do}')
        mod = self.all_models.get_model_from_row(row)
        orb_done = False
        wts_done = False
        if self.do_dummy_run:
            parset = self.all_models.get_parset_from_row(row)
            mod.chi2 = self.dummy_chi2_function(parset)
            mod.kinchi2 = 0.
        else:
            mod.setup_directories()
            orblib = mod.get_orblib()
            orb_done = True
            weight_solver = mod.get_weights(orblib)
            wts_done = True
        all_done = orb_done and wts_done
        time = np.datetime64('now', 'ms')
        output = orb_done, wts_done, mod.chi2, mod.kinchi2, all_done, time
        return output

    def write_output_to_all_models_table(self, rows_to_do, output):
        """write_output_to_all_models_table

        Parameters
        ----------
        rows_to_do : list of ints
            which rows of all models table to fill with output
        output : list
            output returned by Pool.map

        Returns
        -------
        Fills output into the all_models table

        """
        for i, row in enumerate(rows_to_do):
            orb_done, wts_done, chi2, kinchi2, all_done, time = output[i]
            self.all_models.table['orblib_done'][row] = orb_done
            self.all_models.table['weights_done'][row] = wts_done
            self.all_models.table['chi2'][row] = chi2
            self.all_models.table['kinchi2'][row] = kinchi2
            self.all_models.table['all_done'][row] = all_done
            self.all_models.table['time_modified'][row] = time
# end
