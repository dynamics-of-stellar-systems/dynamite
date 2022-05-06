import os
import copy
import logging
import shutil
import numpy as np
from astropy import table
from astropy.io import ascii

from dynamite import weight_solvers as ws
from dynamite import orblib as dyn_orblib

class AllModels(object):
    """All models which have been run so far

    The main attribute ``self.table`` is an Astropy table holding all models run
    so far.

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    from_file : bool
        whether to create this ojbect from a saved `all_models.ecsv` file

    """
    def __init__(self, config=None, from_file=True):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None:
            text = f'{__class__.__name__} needs configuration object, ' \
                   'None provided.'
            self.logger.error(text)
            raise ValueError(text)
        self.config = config
        self.system = config.system
        self.set_filename(config.settings.io_settings['all_models_file'])
        if from_file and os.path.isfile(self.filename):
            self.logger.info('Previous models have been found: '
                        f'Reading {self.filename} into '
                        f'{__class__.__name__}.table')
            self.read_completed_model_file()
        else:
            self.logger.info(f'No previous models (file {self.filename}) '
                        'have been found: '
                        f'Making an empty table in {__class__.__name__}.table')
            self.make_empty_table()

    def set_filename(self, filename):
        """Set the name (including path) for this model

        Parameters
        ----------
        filename : string
            name for this model

        Returns
        -------
        None
            sets ``self.filename``

        """
        outdir = self.config.settings.io_settings['output_directory']
        filename = f'{outdir}{filename}'
        self.filename = filename

    def make_empty_table(self):
        """Make an empty Astropy table

        Returns
        -------
        None
            sets ``self.table``

        """
        names = self.config.parspace.par_names.copy()
        dtype = [np.float64 for n in names]
        # add the columns from legacy version
        names += ['chi2', 'kinchi2', 'time_modified']
        dtype += [np.float64, np.float64, str]
        # add extra columns
        names += ['orblib_done', 'weights_done', 'all_done']
        dtype += [bool, bool, bool]
        # which_iter will record which iteration of parameters a model came from
        names.append('which_iter')
        dtype.append(int)
        # directory will be the model directory name in the models/ directory
        names.append('directory')
        dtype.append('<S256') # little-endian string of max. 256 characters
        self.table = table.Table(names=names, dtype=dtype)

    def read_completed_model_file(self):
        """Read table from file ``self.filename``

        Dealing with incomplete models:
        Models with all_done==False but an existing model_done_staging.ecsv
        will be updated in the table and the staging file will be deleted.
        Models with all_done==False and no existing orblib will be deleted
        from the table and their model directory will be deleted, too.
        The configuration setting reattempt_failures determines how partially
        completed models with all_done==False but existing orblibs are treated:
        If reattempt_failures==True, their orblib_done will be set to True
        and weight solving will be done based on the existing orblibs.
        If reattempt_failures==False, the model and its directory will be
        deleted.
        Note that orbit libraries on disk will not be deleted as they
        may be in use by other models.

        Returns
        -------
        None
            sets ``self.table``

        """
        self.table = ascii.read(self.filename)
        self.logger.debug(f'{len(self.table)} models read '
                          f'from file {self.filename}')

        table_modified = False
        for i, row in enumerate(self.table):
            if not row['all_done']:
                table_modified = True
                mod = self.get_model_from_row(i)
                staging_filename = mod.directory+'model_done_staging.ecsv'
                check1 = os.path.isfile(
                    mod.directory_noml+'datfil/orblib.dat.bz2'
                    )
                check2 = os.path.isfile(
                    mod.directory_noml+'datfil/orblibbox.dat.bz2'
                    )
                check_if_orblibs_present = check1 and check2
                if os.path.isfile(staging_filename):
                    # the model has completed but was not entered in the table
                    staging_file = ascii.read(staging_filename)
                    self.table[i] = staging_file[0]
                    self.logger.info(f'Staging file {staging_filename} '
                                f'used to update {__class__.__name__}.table.')
                    os.remove(staging_filename)
                    self.logger.debug(
                        f'Staging file {staging_filename} deleted.')
                elif check_if_orblibs_present:
                    self.logger.debug(f'Row {i}: orblibs were computed '
                                      'but not weights.')
                    self.table[i]['orblib_done'] = True
                else:
                    self.logger.debug(f'Row {i}: neither orblibs nor '
                                      'weights were completed.')
        # collect failed models to delete (both their directory and table entry)
        to_delete = []
        # if we will reattempt weight solving, only delete models with no orblib
        if self.config.settings.weight_solver_settings['reattempt_failures']:
            for i, row in enumerate(self.table):
                if not row['orblib_done']:
                    to_delete.append(i)
                    self.logger.info('No orblibs calculated for model in '
                                     f'{row["directory"]} - removing row {i}.')
        # otherwise delete any model which is not `all_done`
        else:
            for i, row in enumerate(self.table):
                if not row['all_done']:
                    to_delete.append(i)
                    self.logger.info('No finished model found in '
                                     f'{row["directory"]} - removing row {i}.')
        # do the deletion
        cwd = os.getcwd()
        os.chdir(self.config.settings.io_settings['model_directory'])
        for row in to_delete:
            shutil.rmtree(self.table[row]['directory'])
            self.logger.info(f"Model {row}'s directory "
                             f"{self.table[row]['directory']} removed.")
        os.chdir(cwd)
        self.table.remove_rows(to_delete)
        if table_modified:
            self.save()

    def read_legacy_chi2_file(self, legacy_filename):
        """
        Read the `legacy` AKA schwpy format of chi2 files

        Taken from schw_basics.py, reads in legacy files named similar to
        griddata/_chi2.cat

        Parameters
        -----------
        legacy_filename: string
            the legacy_filename (probably griddata/_chi2.cat)
        """
        ### read the header
        head1 = np.genfromtxt(legacy_filename, max_rows=1)
        Nf = int(head1[0])
        npar = int(head1[1])
        ### read the main matrix
        mpar = np.genfromtxt(legacy_filename,
                             max_rows=Nf,
                             skip_header=1)
        mtest = np.genfromtxt(legacy_filename,
                              max_rows=1,
                              skip_header=Nf + 1)
        len_mtest = len(np.atleast_1d(mtest))
        ### read the last modification time
        if np.mod(Nf, len_mtest) > 0:
            mlast = 1  # add 1 to line counts for the next step, reading the fls
            mtime1 = np.genfromtxt(legacy_filename,
                                   max_rows=int(Nf / len_mtest),
                                   skip_header=Nf + 1)
            mtime2 = np.genfromtxt(legacy_filename,
                                   max_rows=1,
                                   skip_header=Nf + int(Nf / len(mtest)) + 1)
            mtime = np.hstack((np.ravel(mtime1), np.ravel(mtime2)))
        else:
            mlast = 0
            mtime = np.ravel(np.genfromtxt(legacy_filename,
                                           max_rows=int(Nf / len_mtest),
                                           skip_header=Nf + 1))
        ### read the file paths
        fls = np.genfromtxt(legacy_filename,
                            max_rows=Nf,
                            skip_header=Nf + int(Nf / len_mtest) + 1 + mlast,
                            dtype=str)
        self.logger.debug(f'Legacy chi2 file {legacy_filename} read')
        return Nf, npar, mpar.T, mtime.T, fls.T

    def convert_legacy_chi2_file(self, legacy_filename=None):
        """
        Convert `legacy` format of chi2 files

        `legacy` AKA schwpy format were likely called ```griddata/_chi2.cat``.

        Parameters
        -----------
        legacy_filename: string
            the legacy_filename (probably griddata/_chi2.cat)
        """
        # TODO: (maybe...?)
        # make more general if legacy parameters have different names
        Nf, npar, mpar, mtime, fls = self.read_legacy_chi2_file(legacy_filename)
        # legacy parameter names are fixed: bh, dc, etc...
        mods = table.Table()
        mods['bh'] = mpar[0,:]
        mods['dc'] = mpar[1,:]
        mods['f'] = mpar[2,:]
        mods['q'] = mpar[3,:]
        mods['p'] = mpar[4,:]
        mods['u'] = mpar[5,:]
        mods['ml'] = mpar[6,:]
        mods['chi2'] = mpar[7,:]
        mods['kinchi2'] = mpar[8,:]
        mods['time_modified'] = mtime
        # currently fls are strings "{param_directory}/{ml_directory}/nn"
        # cleaner to have just the directry name instead...?
        direcs = [file0[:-2] for file0 in fls]
        mods['directory'] = table.Column(direcs, dtype='<U100')
        # add any extra columns we want to have in the table
        mods['ics_done'] = True
        mods['orblib_done'] = True
        mods['weights_done'] = True
        mods['all_done'] = True
        # which_iter will record which iteration of parameters a model came from
        mods['which_iter'] = 0  # was not recorded for schwpy so set to 0
        mods.write(self.filename, format='ascii.ecsv')
        self.logger.debug(f'Legacy chi2 file {legacy_filename} converted ' + \
                          f'to {self.filename}')

    def get_parset_from_row(self, row_id):
        """Get a parset set given a table row

        Parameters
        ----------
        row_id : int
            which row

        Returns
        -------
        list
            a list of ``dyn.parspace.Parameter`` objects

        """
        parset = self.table[row_id][self.config.parspace.par_names]
        return parset

    def get_model_from_parset(self, parset):
        """Get the ``Model`` from a parset

        Parameters
        ----------
        parset : list
            a list of ``dyn.parspace.Parameter`` objects

        Returns
        -------
        a ``dyn.model.Model`` object

        """
        for idx, row in enumerate(self.table[self.config.parspace.par_names]):
            if np.allclose(tuple(parset), tuple(row)):
                mod = self.get_model_from_row(idx)
                break
        else:
            text = f'parset not in all_models table. parset={parset}, ' \
                   f'all_models table: {self.table}'
            self.logger.error(text)
            raise ValueError(text)
        return mod

    def get_model_from_row(self, row_id):
        """Get a ``Model`` given a table row

        Parameters
        ----------
        row_id : int
            which row

        Returns
        -------
        a ``dyn.model.Model`` object

        """
        parset = self.get_parset_from_row(row_id)
        mod = Model(config=self.config,
                    parset=parset,
                    directory=self.table['directory'][row_id])
        return mod

    def get_row_from_model(self, model=None):
        """Get the table row for a ``Model``

        Parameters
        ----------
        model : a ``Model`` object

        Raises
        ------
        ValueError
            If the model is not found in ``self.table``.

        Returns
        -------
        row_id : int
            The ``self.table`` row index of model.

        """
        row_comp = tuple(model.parset[self.config.parspace.par_names])
        for row_id, row in enumerate(
                                self.table[self.config.parspace.par_names]):
            if np.allclose(row_comp, tuple(row)):
                break
        else:
            text = 'Cannot find model in all_models table.'
            self.logger.error(text)
            raise ValueError(text)
        return row_id

    def get_ml_of_original_orblib(self, model_id):
        """Get ``ml`` of model number model_id's original orblib

        The original ``ml`` is required to rescale orbit libraries for rescaled
        potentials. This method searches ``self.table``, the all_models table.

        Parameters
        ----------
        model_id : int
            The ``self.table`` row index of the model.

        Raises
        ------
        ValueError
            If the ``ml`` parameter is not in the parameter space or the
            model's orblib cannot be found.

        Returns
        -------
        ml_orblib : float
            the original ``ml``

        """
        orblib_parameters = self.config.parspace.par_names[:]
        ml = 'ml'
        try:
            orblib_parameters.remove(ml)
        except:
            self.logger.error(f"Parameter '{ml}' not found - check "
                              "implementation")
            raise
        row_comp = tuple(self.table[orblib_parameters][model_id])
        for row_id, row in enumerate( \
                                 self.table[orblib_parameters][:model_id+1]):
            if np.allclose(row_comp, tuple(row)):
                ml_orblib = self.table['ml'][row_id]
                self.logger.debug(f'Orblib of model #{model_id} has original '
                                  f'ml value of {ml_orblib} '
                                  f'(model #{row_id}).')
                break
        else:
            text = f'Cannot find orblib for model #{model_id} in ' \
                   'all_models table.'
            self.logger.error(text)
            raise ValueError(text)
        return ml_orblib

    def get_model_velocity_scaling_factor(self, model_id=None, model=None):
        """Get the model's velocity scaling factor

        Returns sqrt(model_ml/original_orblib_ml).
        The model can be either given by its row id in ``self.table`` or
        as a ``Model`` object. Note that the parameters model_id and model
        are mutually exclusive.

        Parameters
        ----------
        model_id : int compatible
            The model's row id in ``self.table``.
        model : a ``Model`` object

        Raises
        ------
        ValueError
            If not exactly one of model_id and model are supplied.

        Returns
        -------
        scaling_factor : float
            The model's velocity scaling factor
            sqrt(model_ml/original_orblib_ml).

        """
        if model_id is None and isinstance(model, Model):
            model_id = self.get_row_from_model(model)
        elif not (model_id == int(model_id) and model is None):
            text = 'Need to pass either model_id (int) or model (Model).'
            self.logger.error(text)
            raise ValueError(text)
        ml_orblib = self.get_ml_of_original_orblib(model_id)
        scaling_factor = np.sqrt(self.table['ml'][model_id]/ml_orblib)
        return scaling_factor

    def save(self):
        """Save the all_models table

        """
        self.table.write(self.filename, format='ascii.ecsv', overwrite=True)
        self.logger.debug(f'Model table written to file {self.filename}')

    def get_best_n_models(self, n=10, which_chi2=None):
        """Get the best n models so far

        Parameters
        ----------
        n : int, optional
            How many models to get. The default is 10.
        which_chi2 : str, optional
            Which chi2 is used for determining the best models. Must be
            None, chi2, or kinchi2. If None, the setting from the
            configuration file will be used. The default is None.

        Raises
        ------
        ValueError
            If which_chi2 is neither None, chi2, nor kinchi2.

        Returns
        -------
        a new ``astropy.table`` object holding the best n models

        """
        if which_chi2 is None:
            which_chi2 = \
                self.config.settings.parameter_space_settings['which_chi2']
        if which_chi2 not in ('chi2', 'kinchi2'):
            text = 'which_chi2 needs to be chi2 or kinchi2, ' \
                   f'but it is {which_chi2}'
            self.logger.error(text)
            raise ValueError(text)
        table = copy.deepcopy(self.table)
        table.sort(which_chi2)
        table = table[:n]
        return table

    def get_mods_within_chi2_thresh(self, which_chi2=None, delta=None):
        """Get models within delta threshold of best

        Parameters
        ----------
        which_chi2 : str, optional
            Which chi2 is used for determining the best models. Must be
            None, chi2, or kinchi2. If None, the setting from the
            configuration file will be used. The default is None.
        delta : float, optional
            The threshold value. Models with (kin)chi2 values differing
            from the opimum by at most delta will be returned. If none,
            models within 10% of the optimal value will be returned.
            The default is None.

        Raises
        ------
        ValueError
            If which_chi2 is neither None, chi2, nor kinchi2.

        Returns
        -------
        a new ``astropy.table`` object holding the ''delta-best'' models

        """
        if which_chi2 is None:
            which_chi2 = \
                self.config.settings.parameter_space_settings['which_chi2']
        if which_chi2 not in ('chi2', 'kinchi2'):
            text = 'which_chi2 needs to be chi2 or kinchi2, ' \
                   f'but it is {which_chi2}'
            self.logger.error(text)
            raise ValueError(text)
        chi2_min = min(self.table[which_chi2])
        if delta is None:
            delta = chi2_min * 0.1
        models = self.table[self.table[which_chi2] <= chi2_min+delta]
        return models


class Model(object):
    """A DYNAMITE model.

    The model can be run by running the methods (i) get_orblib, (ii) get_weights
    and (iii) (in the future) do_orbit_colouring. Running each of these methods
    will return the appropriate object, e.g. model.get_orblib() --> returns an
    OrbitLibrary object model.get_weights(...) --> returns a WeightSolver object

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    parset : row of an Astropy Table
        contains the values of the potential parameters for this model
    directory : str
        The model directory name (without path). If None or not specified,
        the all_models_file will be searched for the directory name. If the
        all_models file does not exist, the model directory will be set to
        ``orblib_000_000/ml{ml}``.

    Returns
    -------
    Nothing returned. Attributes holding outputs are are added to the
    object when methods are run.

    """
    def __init__(self, config=None, parset=None, directory=None):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None or parset is None:
            text = f'{__class__.__name__} needs configuration object and ' \
                   'parameter set.'
            self.logger.error(text)
            raise ValueError(text)
        self.config = config
        self.check_parset(config.parspace, parset)
        self.parset = parset
        # directory of the input kinematics
        if directory is None:
           self.directory = self.get_model_directory()
        else:
           self.directory=self.config.settings.io_settings['output_directory']\
                          + 'models/' + directory
        self.logger.debug(f'Model directory string: {self.directory}')
        self.directory_noml=self.directory[:self.directory[:-1].rindex('/')+1]
        self.logger.debug('Model directory string up to ml: '
                          f'{self.directory_noml}')

    def get_model_directory(self):
        """get the name of this model's output directory

        Returns
        -------
        string

        """
        output_directory = self.config.settings.io_settings['output_directory']
        directory = output_directory + 'models/'
        models_file = output_directory \
                      + self.config.settings.io_settings['all_models_file']
        try:
            all_models = ascii.read(models_file)
            self.logger.debug(f'Setting model dir from file {models_file}...')
        except FileNotFoundError:
            sformat = self.config.system.parameters[0].sformat # ml's format
            ml_dir = f"/ml{self.parset['ml']:{sformat}}/"
            directory += f'orblib_000_000{ml_dir}'
            self.logger.info(f'The all_models file {models_file} does not '
                             f'exist - model directory set to {directory}.')
            return directory #######################################
        except:
            self.logger.error('Error reading all_models file. '
                              'Cannot set model directory.')
            raise
        for idx, parset in enumerate(all_models[self.config.parspace.par_names]):
            if np.allclose(tuple(parset),tuple(self.parset)):
                dir_string = all_models['directory'][idx]
                if dir_string == 'None': # yes, really...
                    text = f'Something went wrong reading the model # {idx} ' \
                           'directory. Model not yet computed?'
                    self.logger.error(text)
                    raise ValueError(text)
                directory += dir_string
                break
        else:
            text = f'Cannot set model directory: parset {self.parset} ' \
                   f'not found in {models_file}.'
            self.logger.error(text)
            raise ValueError(text)
        self.logger.debug(f'...model directory {directory} read from file.')
        return directory

    def create_model_directory(self, path):
        """make a directory if it does not yet exist

        Parameters
        ----------
        path : string
            directory path to make
        """
        if not os.path.exists(path):
            os.makedirs(path)
            self.logger.debug(f'Created directory {path}')
        else:
            self.logger.debug(f'Using existing directory {path}')

    def setup_directories(self):
        """setup directories
        """
        # create self.directory if it doesn't exist
        self.create_model_directory(self.directory)
        # and also the model directories
        self.create_model_directory(self.directory_noml+'infil/')
        self.create_model_directory(self.directory_noml+'datfil/')

    def get_orblib(self):
        """Make the orbit library

        Returns
        -------
        a ``dyn.orblib.OrbitLibrary`` object

        """
        orblib = dyn_orblib.LegacyOrbitLibrary(
                config=self.config,
                mod_dir=self.directory_noml,
                parset=self.parset)
        orblib.get_orblib()
        return orblib

    def get_weights(self, orblib=None):
        """Get the orbital weights

        Parameters
        ----------
        orblib : a ``dyn.orblib.OrbitLibrary`` object

        Returns
        -------
        a ``dyn.weight_solver.WeightSolver`` object

        """
        ws_type = self.config.settings.weight_solver_settings['type']
        if ws_type=='LegacyWeightSolver':
            weight_solver = ws.LegacyWeightSolver(
                    config=self.config,
                    directory_with_ml=self.directory)
        elif ws_type=='NNLS':
            weight_solver = ws.NNLS(
                    config=self.config,
                    directory_with_ml=self.directory)
        else:
            raise ValueError('Unknown WeightSolver type')
        weights, chi2_tot, chi2_kin = weight_solver.solve(orblib)
        self.chi2 = chi2_tot # instrinsic/projected mass + GH coeeficients 1-Ngh
        self.kinchi2 = chi2_kin # GH coeeficients 1-Ngh
        self.weights = weights
        return weight_solver

    def check_parset(self, parspace, parset):
        """
        Validate a parset

        Given parameter values in parset, the validate_parspace method of
        the parameter space is executed. If a parameter exists in parspace
        but not in parset, a warning will be issued and the parameter
        will remain unchanged. If parset tries to set the value of
        a parameter not existing in parspace, an exception will be raised.
        Validating relies on exceptions raised by validate_parspace.

        Parameters
        ----------
        parspace : ``dyn.parameter_space.ParameterSpace``
            A list of parameter objects.
        parset : row of an Astropy Table
            Contains parameter values to be checked against the settings in
            parspace.

        Raises
        ------
        ValueError
            If at least one parameter in parset is unknown to parspace.

        Returns
        -------
        None.

        """
        parspace_copy = copy.deepcopy(parspace)
        parspace_par_names = [p.name for p in parspace_copy]
        unknown_pars = [p for p in parset.colnames
                          if p not in parspace_par_names]
        if len(unknown_pars) > 0:
            text = f"Parset parameters {unknown_pars} don't exist in parset."
            self.logger.error(text)
            raise ValueError(text)
        for par_idx, par_name in enumerate(parspace_par_names):
            if par_name not in parset.colnames:
                self.logger.warning(f'Parspace parameter {par_name} '
                                    'unchanged (not in parset).')
            else:
                par = parspace_copy[par_idx]
                par.par_value = parset[par_name]
        parspace_copy.validate_parspace()
        self.logger.debug('parset validated against parspace.')

# end
