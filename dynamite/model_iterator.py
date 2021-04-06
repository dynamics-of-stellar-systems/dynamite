import model
import parameter_space
import plotter
import os
import numpy as np
import logging
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
        for iter in range(n_max_iter):
            total_iter_count = previous_iter + iter
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
            status = model_inner_iterator.run_iteration(iter)
            self.make_in_progress_plots(settings, iter)

    def make_in_progress_plots(self, settings, iteration=None,
                               chi2_progress=None,
                               chi2_plot=None,
                               kin_map=None):
        """
        Creates three plots: (kin)chi2 vs. model id, (kin)chi2 and non-fixed
        parameters ("chi2 plot"), kinematic map of best fit model so-far.
        The parameter space settings in the config file determine whether
        chi2 or kinchi2 is used. Will choose file names automatically and
        append the iteration counter to avoid duplicate file names.

        Parameters
        ----------
        settings : Settings object
            Needed for plot directory and which_chi2 setting.
        iteration : int, optional
            Iteration counter; mandatory for automatic file names.
            The default is None.
        chi2_progress : str, optional
            File name of the (kin)chi2 vs. model id plot. If None, the
            file name will be created automatically. The default is None.
        chi2_plot : str, optional
            File name of the "chi2 plot". If None, the
            file name will be created automatically. The default is None.
        kin_map : str, optional
            File name of the kinematic map. If None, the
            file name will be created automatically. The default is None.

        Raises
        ------
        ValueError
            Will be raised if at least one of the file names is None
            and iteration is not an integer.

        Returns
        -------
        None.

        """
        if (chi2_progress is None or chi2_plot is None or kin_map is None) \
            and type(iteration) is not int:
                text = 'iteration must be an integer when automatic file ' + \
                       'names are used.'
                self.logger.error(text)
                raise ValueError(text)
        plot_dir = settings.io_settings['output_directory'] + '/plots/'
        which_chi2 = settings.parameter_space_settings['which_chi2']
        if chi2_progress is None:
            chi2_progress = f'{which_chi2}_progress_plot_{iteration}'
        chi2_progress = plot_dir + chi2_progress
        if chi2_plot is None:
            chi2_plot = f'{which_chi2}_plot_{iteration}'
        chi2_plot = plot_dir + chi2_plot
        if kin_map is None:
            kin_map = f'kinematics_map_{iteration}'
        kin_map = plot_dir + kin_map
        self.delete_if_exists([chi2_progress, chi2_plot, kin_map])

        self.the_plotter.make_chi2_vs_model_id_plot().savefig(chi2_progress)
        self.the_plotter.make_chi2_plot().savefig(chi2_plot)
        self.the_plotter.plot_kinematic_maps(cbar_lims='data').savefig(kin_map)

    def delete_if_exists(self, files):
        """
        Given a file name or a list or tuple of file names, this method
        will check if the file(s) exist and if so, remove it/them.

        Parameters
        ----------
        files : str or list or tuple
            File name including the path or list or tuple of file names.

        Raises
        ------
        ValueError
            Will be raised if files is neither a string nor a list nor
            a tuple.

        Returns
        -------
        None.

        """
        if type(files) == list or type(files) == tuple:
            for f in files:
                if os.path.isfile(f):
                    os.remove(f)
        elif type(files) == str:
            if os.path.isfile(f):
                os.remove(f)
        else:
            text = 'files must be of type str, list, or tuple.'
            self.logger.error(text)
            raise ValueError(text)


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
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
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
        self.all_models.save() # save all_models table once parameters are added
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
            self.all_models.save() # save all_models table once models are run
        return self.par_generator.status

    def create_and_run_model(self, input):
        i, row = input
        self.logger.info(f'... running model {i+1} out of {self.n_to_do}')
        # extract the parameter values
        parset0 = self.all_models.table[row]
        parset0 = parset0[self.parspace.par_names]
        # create and run the model
        mod0 = model.Model(system=self.system,
                           settings=self.settings,
                           parspace=self.parspace,
                           parset=parset0)
        if self.do_dummy_run:
            mod0.chi2 = self.dummy_chi2_function(parset0)
            mod0.kinchi2 = 0.
        else:
            mod0.setup_directories()
            orblib = mod0.get_orblib()
            orb_done = True
            weight_solver = mod0.get_weights(orblib)
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
