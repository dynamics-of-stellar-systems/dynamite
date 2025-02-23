import os
import logging
import numpy as np
from astropy import table
from pathos.multiprocessing import Pool
import matplotlib.pyplot as plt

from dynamite import parameter_space
from dynamite import plotter

class ModelIterator(object):
    """Iterator for models

    Creating this ``ModelIterator`` object will (i) generate parameters sets,
    (ii) run models for those parameters, (iii) check stopping criteria, and
    iterate this procedure till a stopping criterion is met. This is implemented
    by creating a ``ModelInnerIterator`` object whose ``run_iteration`` method
    is called a number of times.

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    model_kwargs : dict
        other kewyord argument required for this model
    do_dummy_run : Bool
        whether this is a dummy run - if so, dummy_chi2_funciton is executed
        instead of the model (for testing!)
    dummy_chi2_function : function
        a function of model parameters to be executed instead of the real model
    plots : bool
        whether or not to make plots

    """
    def __init__(self,
                 config=None,
                 model_kwargs={},
                 do_dummy_run=None,
                 dummy_chi2_function=None,
                 plots=True):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None:
            text = f'{__class__.__name__} needs configuration object, ' \
                   'None provided.'
            self.logger.error(text)
            raise ValueError(text)
        self.config = config
        parameter_space_settings = config.settings.parameter_space_settings
        stopping_crit = parameter_space_settings['stopping_criteria']
        # get specified parameter generator
        par_generator_type = parameter_space_settings['generator_type']
        kwargs = {'parspace_settings':parameter_space_settings}
        par_generator = getattr(parameter_space,
                                par_generator_type)(config.parspace, **kwargs)

        if plots:
            the_plotter = plotter.Plotter(config)

        iterator = config.settings.multiprocessing_settings['modeliterator']
        model_inner_iterator = globals()[iterator](
            config=config,
            par_generator=par_generator,
            do_dummy_run=do_dummy_run,
            dummy_chi2_function=dummy_chi2_function)
        if len(config.all_models.table)>0:
            previous_iter = np.max(config.all_models.table['which_iter'])
        else:
            previous_iter = -1
        status = {}
        status['stop'] = False
        # if configured, re-calculate weights for past models where weight
        # calculation failed
        if config.settings.weight_solver_settings['reattempt_failures']:
            self.reattempt_failed_weights()
        iteration = 1
        while iteration <= stopping_crit['n_max_iter']:
            total_iter = previous_iter + iteration
            n_models_done = np.sum(config.all_models.table['all_done'])
            if n_models_done >= stopping_crit['n_max_mods']:
                status['n_max_mods_reached'] = True
                status['stop'] = True
            if status['stop'] is True:
                self.logger.info(f'Stopping at iteration {total_iter}')
                self.logger.debug(status)
                break
            if total_iter > 0:
                self.logger.info(f'{par_generator_type}: iteration '
                                 f'{total_iter}')
                iteration += 1
            else:
                self.logger.info(f'{par_generator_type}: iterations 0 and 1')
                iteration += 2
            status = model_inner_iterator.run_iteration()
            if plots and not status['last_iter_added_no_new_models']:
                try:
                    self.chi2_vs_model_id_plot = \
                        the_plotter.make_chi2_vs_model_id_plot()
                    self.chi2_plot = the_plotter.make_chi2_plot()
                    self.kinematic_maps = \
                        the_plotter.plot_kinematic_maps(kin_set='all',
                                                        cbar_lims='default')
                    plt.close('all')  # just to make sure...
                except:
                    self.logger.warning(f'Iteration {total_iter}: '
                                        'plotting failed!')

    def get_plots(self):
        """
        Returns the latest iteration's plots as figure objects.

        Returns
        -------
        tuple of matplotlib.pyplot.figure:
            matplotlib.pyplot.figure: chi2 vs. model id plot
            matplotlib.pyplot.figure: chisquare plot
            (matplotlib.pyplot.figure, str): kinematic maps of best model so
            far, kinematics name

        """
        chi2_vs_model_id_plot = self.chi2_vs_model_id_plot \
            if hasattr(self, 'chi2_vs_model_id_plot') else None
        chi2_plot = self.chi2_plot if hasattr(self, 'chi2_plot') else None
        kinematic_maps = self.kinematic_maps \
            if hasattr(self, 'kinematic_maps') else None
        return chi2_vs_model_id_plot, chi2_plot, kinematic_maps

    def reattempt_failed_weights(self):
        config = self.config
        rows_with_orbits_but_no_weights = \
            [i for i,t in enumerate(config.all_models.table) \
                if t['orblib_done'] \
                    and not t['weights_done'] \
                    and not t['all_done']]
        if len(rows_with_orbits_but_no_weights) > 0:
            to_do = config.all_models.table[rows_with_orbits_but_no_weights]
            to_do = to_do['directory']
            self.logger.info(f'Reattempting weight solving: models {to_do}.')
            n_proc = config.settings.multiprocessing_settings['ncpus_weights']
            with Pool(n_proc) as p:
                output = p.map(self.get_missing_weights,
                               rows_with_orbits_but_no_weights)
            for i, row in enumerate(rows_with_orbits_but_no_weights):
                chi2, kinchi2, kinmapchi2, time = output[i]
                config.all_models.table[row]['chi2'] = chi2
                config.all_models.table[row]['kinchi2'] = kinchi2
                config.all_models.table[row]['kinmapchi2'] = kinmapchi2
                config.all_models.table[row]['time_modified'] = time
                directory = config.all_models.table[row]['directory']
                if not all(np.isnan([chi2, kinchi2, kinmapchi2])):
                    config.all_models.table[row]['weights_done'] = True
                    config.all_models.table[row]['all_done'] = True
                    self.logger.info('Reattempting weight solving for model '
                                     f'in row {i}, {directory=} successful.')
                else:
                    self.logger.info('Reattempting weight solving for model '
                                     f'in row {i}, {directory=} failed.')
            config.all_models.save()

    def get_missing_weights(self, row):
        mod = self.config.all_models.get_model_from_row(row)
        mod.setup_directories()
        self.logger.debug(f'Reattempting weight solving for model {row}, '
                          f'directory={mod.directory}.')
        orblib = mod.get_orblib()
        weight_solver = mod.get_weights(orblib)
        time = str(np.datetime64('now', 'ms'))
        return mod.chi2, mod.kinchi2, mod.kinmapchi2, time

class ModelInnerIterator(object):
    """Class to run all models in a single iteration.

    Uses ``pathos.multiprocessing.Pool`` to execute the models

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    par_generator : a ``dyn.parameter_space.ParameterGenerator`` child object
    do_dummy_run : Bool
        whether this is a dummy run - if so, dummy_chi2_function is executed
        instead of the model (for testing!)
    dummy_chi2_function : function
        a function of model parameters to be executed instead of the real model

    """
    def __init__(self,
                 config=None,
                 par_generator=None,
                 do_dummy_run=False,
                 dummy_chi2_function=None):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.config = config
        self.system = config.system
        self.all_models = config.all_models
        self.orblib_parameters = config.parspace.par_names[:]
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
        self.ncpus = config.settings.multiprocessing_settings['ncpus']
        self.n_to_do = 0

    def run_iteration(self, split_orblib_weights=False):
        """Run one iteration step

        Executes one iteration step: run all models in self.all_models.table
        for which all_done == False. The model runs (1) build the orbit
        library and (2) execute weight_solver. The models are run in parallel
        threads as defined in the ncpus parameter in the configuration file.
        In case multiple models comprise the same (new) orbit library it is
        ensured that they are calculated only once, avoiding conflicting
        threads.

        Parameters
        ----------
        split_orblib_weights : bool, optional
            If True, first the orbit libraries are calculated in self.ncpus
            parallel pools, then the weights in self.ncpus_weights parallel
            pools. The default is False.

        Returns
        -------
        dict
            ParameterGenerator.status.

        """

        self.par_generator.generate(current_models=self.all_models)
        self.all_models.save() # save all_models table once parameters are added
        if not self.par_generator.status['stop']:
            # find new models which are those with an empty directory string
            # and always correspond to either iteration 0+1 or
            # the current iteration if iteration > 1
            rows_to_do = np.where(self.all_models.table['directory'] == '')
            rows_to_do = rows_to_do[0]
            self.logger.debug(f'rows_to_do: {rows_to_do}.')
            self.n_to_do = len(rows_to_do)
            # rows_to_do_orblib are the rows that need orblib and weight_solver
            rows_to_do_orblib=[i for i in rows_to_do if self.is_new_orblib(i)]
            n_orblib = len(rows_to_do_orblib)
            # rows_to_do_ml are the rows that need weight_solver only
            rows_to_do_ml=[i for i in rows_to_do if i not in rows_to_do_orblib]
            self.assign_model_directories(rows_to_do_orblib, rows_to_do_ml)
            # save all_models here - as it is useful to have directories saved
            # even if the run fails and we don't reach the next save
            self.all_models.save()
            if split_orblib_weights:  # first the orblibs, then all weights
                do_orblib, do_weights = True, False  # orblibs only
                input_list_orblib = [i + (do_orblib, do_weights)
                                     for i in enumerate(rows_to_do_orblib)]
                if len(input_list_orblib) > 0:
                    with Pool(self.ncpus) as p:
                        output = p.map(self.create_and_run_model,
                                       input_list_orblib)
                    self.write_output_to_all_models_table(rows_to_do_orblib,
                                                          output)
                    self.all_models.save()
                do_orblib, do_weights = False, True  # all the weights
                input_list_ml = [i + (do_orblib, do_weights)
                                 for i in enumerate(rows_to_do)]
                if len(rows_to_do) > 0:
                    with Pool(self.ncpus_weights) as p:
                        output = p.map(self.create_and_run_model, input_list_ml)
                    self.write_output_to_all_models_table(rows_to_do, output)
            else:  # first the orblibs + their weights, then remaining weights
                input_list_orblib = list(enumerate(rows_to_do_orblib,
                                                   start=n_orblib))
                input_list_ml = list(enumerate(rows_to_do_ml))
                if len(input_list_orblib) + len(input_list_ml) > 0:
                    with Pool(self.ncpus) as p:
                        if len(input_list_orblib) > 0:
                            output_orblib = p.map(self.create_and_run_model,
                                                  input_list_orblib)
                        if len(input_list_ml) > 0:
                            output_ml = p.map(self.create_and_run_model,
                                              input_list_ml)
                    if len(input_list_orblib) > 0:
                        self.write_output_to_all_models_table(rows_to_do_orblib,
                                                              output_orblib)
                    if len(input_list_ml) > 0:
                        self.write_output_to_all_models_table(rows_to_do_ml,
                                                              output_ml)
            self.all_models.save()  # save all_models table once models are run
            self.logger.info('Iteration done, '
                             f'{self.n_to_do} model(s) calculated.')
            self.delete_staging_files(rows_to_do) # delete all staging files
        return self.par_generator.status

    def delete_staging_files(self, rows):
        """
        Deletes staging files.

        Parameters
        ----------
        rows : iterable of ints
            The all_models table rows indicating models whose staging files
            are to be deleted.

        Returns
        -------
        n_files : int
            Number of staging files deleted.

        """
        for row in rows:
            f_name = self.all_models.get_model_from_row(row).directory + \
                'model_done_staging.ecsv'
            if os.path.isfile(f_name):
                os.remove(f_name)
            else:
                self.logger.warning(f'Strange: {f_name} does not exist.')
        n_files = len(rows)
        self.logger.info(f'{n_files} staging file(s) deleted.')
        return n_files

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
            # self.logger.debug('Orblib exists above in table: '
            #                   f'{row_data}.')
            is_new = False
        else:
            # self.logger.debug(f'New orblib: {row_data}.')
            is_new = True
        return is_new

    def assign_model_directories(self, rows_orblib=[], rows_ml=[]):
        """
        Assigns model directories in all_models.table.

        Models indexed by rows_orblib:
        The model directories follow the pattern orblib_xxx_yyy/mlzz.zz/ where
        xxx is the iteration number, yyy a consecutive number of that
        iteration's orbit library, and zz.zz is the value of the models'
        ml parameter in the format given in its sformat attribute.

        Models indexed by rows_ml:
        These models re-use an existing orbit library. Hence, their directory
        strings re-use an existing orblib_xxx_yyy part and get augmented with
        the appropriate /mlzz.zz/.

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
        # new orblib model directories
        nodir = ''
        for row in rows_orblib:
            iteration = self.all_models.table[row]['which_iter']
            t = self.all_models.table[:row]
            n = np.sum((t['which_iter']==iteration) & (t['directory']!=nodir))
            orblib_dir = f'orblib_{iteration:03d}_{n:03d}'
            self.all_models.table[row]['directory'] = orblib_dir
        # existing orblib directories
        orblib_data = self.all_models.table[self.orblib_parameters]
        for row in rows_ml:
            row_data = orblib_data[row]
            for idx, orblib in enumerate(orblib_data[:row]):
                if np.allclose(tuple(row_data), tuple(orblib)):
                    orblib_dir = self.all_models.table[idx]['directory']
                    if orblib_dir[-1] == '/': # need to strip ml subdirectory?
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

    def create_and_run_model(self, data_input):
        """Main method to create and run a model

        Parameters
        ----------
        data_input : tuple of length 2 or 4
            len(input)==2:
            (i, row) where i is the index of a model in this iteration, and row
            is the row index of the all_models table for this model
            Both the orblib and the weights will be computed.
            len(input==4):
            (i, row, get_orblib, get_weights) where i and row are as above,
            get_orblib==True if the orblib needs to be computed
            get_weights==True if the weights need to be computed

        Returns
        -------
        tuple
            all the output for this model, bundled up in a tuple

        """
        if len(data_input) == 2:
            i, row = data_input
            get_orblib, get_weights = True, True
        elif len(data_input) == 4:
            i, row, get_orblib, get_weights = data_input
        else:
            msg = 'Unexpected input, need a tuple of length 2 or 4.'
            self.logger.error(msg)
            raise ValueError(msg)
        mod = self.all_models.get_model_from_row(row)
        self.logger.info(f'... running model {i+1} out of {self.n_to_do}: '
                         f'{mod.directory}.')
        orb_done = False
        wts_done = False
        if self.do_dummy_run:
            parset = self.all_models.get_parset_from_row(row)
            mod.chi2 = self.dummy_chi2_function(parset)
            mod.kinchi2 = np.nan
            mod.kinmapchi2 = np.nan
        else:
            if not (get_orblib or get_weights):
                msg = 'Nothing to run, specify get_orblib and/or get_weights'
                self.logger.error(msg)
                raise ValueError(msg)
            mod.setup_directories()
            if not get_orblib and self.is_new_orblib(row) and \
                              not self.all_models.table['orblib_done'][row]:
                msg = f'Unexpected: orbit library in row {row}, ' \
                      f'directory {mod.directory} not existing! ' \
                      'Will calculate it (beware of multiprocessing - use ' \
                      'the restart feature or ncpus=1 in case of problems)...'
                self.logger.warning(msg)
            cwd = os.getcwd()
            try:
                orblib = mod.get_orblib()
                orb_done = True
                if get_weights:
                    weight_solver = mod.get_weights(orblib)
                    if not np.isnan(mod.weights[0]):
                        wts_done = True
                else:
                    mod.chi2, mod.kinchi2, mod.kinmapchi2 \
                        = np.nan, np.nan, np.nan
            except RuntimeError:
                os.chdir(cwd)
                mod.chi2, mod.kinchi2, mod.kinmapchi2 = np.nan, np.nan, np.nan
                w_txt = f'Model {i+1} (row {row}, ' \
                        f'directory {mod.directory}): get_orblib ' \
                        + ('or get_weights ' if get_weights else '')+'failed.'\
                        + (' all chi2 values set to nan!' \
                           if get_weights else '')
                self.logger.warning(w_txt)
        all_done = orb_done and wts_done
        time = str(np.datetime64('now', 'ms'))
        # Build and write model_done_staging.ecsv
        current_model_row = table.Table(self.all_models.table[row])
        for name, value in zip(
                ['orblib_done','weights_done','chi2',
                 'kinchi2','kinmapchi2','all_done','time_modified'],
                [orb_done, wts_done, mod.chi2,
                 mod.kinchi2, mod.kinmapchi2, all_done, time]):
            current_model_row[name][0] = value
        file_name = mod.directory + 'model_done_staging.ecsv'
        current_model_row.write(file_name, format='ascii.ecsv', overwrite=True)
        self.logger.info(f'Model {i+1}: {file_name} written.')
        output = orb_done, wts_done, mod.chi2, \
                 mod.kinchi2, mod.kinmapchi2, all_done, time
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
            orb_done, wts_done, \
                chi2, kinchi2, kinmapchi2, all_done, time = output[i]
            self.all_models.table['orblib_done'][row] = orb_done
            self.all_models.table['weights_done'][row] = wts_done
            self.all_models.table['chi2'][row] = chi2
            self.all_models.table['kinchi2'][row] = kinchi2
            self.all_models.table['kinmapchi2'][row] = kinmapchi2
            self.all_models.table['all_done'][row] = all_done
            self.all_models.table['time_modified'][row] = time


class SplitModelIterator(ModelInnerIterator):
    """Class to run all models in a single iteration.

    First calculates the orbit libraries and then runs weight solving.
    Uses ``pathos.multiprocessing.Pool`` to execute the models.
    Orbit integration uses a pool of ncpus parallel processes, weight
    solving uses ncpus_weights from the Configuration object.

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    par_generator : a ``dyn.parameter_space.ParameterGenerator`` child object
    do_dummy_run : Bool
        whether this is a dummy run - if so, dummy_chi2_funciton is executed
        instead of the model (for testing!)
    dummy_chi2_function : function
        a function of model parameters to be executed instead of the real model

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.ncpus_weights = \
            self.config.settings.multiprocessing_settings['ncpus_weights']

    def run_iteration(self):
        """Execute one iteration step

        Calls the parameter generator and (a) calculates all new orbit
        libraries and consecutively (b) does the weight solving.
        (a) and (b) are run in their respective parallel pools
        as defined by the ncpus and ncpus_weights parameters, respectively.

        Returns
        -------
        dict
            ParameterGenerator.status.

        """
        return super().run_iteration(split_orblib_weights=True)

# end
