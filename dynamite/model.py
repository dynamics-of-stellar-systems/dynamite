import os
import copy
import logging
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
    system : a ``dyn.physical_system.System`` object
    from_file : bool
        whether to create this ojbect from a saved `all_models.ecsv` file
    settings : a ``dyn.config_reader.Settings`` object
    parspace : a ``dyn.parameter_space.parspace`` object

    """
    def __init__(self,
                 system=None,
                 from_file=True,
                 settings=None,
                 parspace=None):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.system = system
        self.settings = settings
        self.set_filename(settings.io_settings['all_models_file'])
        self.parspace = parspace
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
        outdir = self.settings.io_settings['output_directory']
        filename = f'{outdir}{filename}'
        self.filename = filename

    def make_empty_table(self):
        """Make an empty Astropy table

        Returns
        -------
        None
            sets ``self.table``

        """
        names = self.parspace.par_names.copy()
        dtype = [np.float64 for n in names]
        # add the columns from legacy version
        names += ['chi2', 'kinchi2', 'time_modified']
        dtype += [np.float64, np.float64, '<M8[ms]']
        # add extra columns
        names += ['orblib_done', 'weights_done', 'all_done']
        dtype += [bool, bool, bool]
        # which_iter will record which iteration of parameters a model came from
        names.append('which_iter')
        dtype.append(int)
        # directory will be the model directory name in the models/ directory
        names.append('directory')
        dtype.append(np.object)
        self.table = table.Table(names=names, dtype=dtype)

    def read_completed_model_file(self):
        """read table from file ``self.self.filename``

        Returns
        -------
        None
            sets ``self.table``

        """
        self.table = ascii.read(self.filename)
        self.logger.debug(f'Models read from file {self.filename}')

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
        Nf = int(head1[0]);
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
        parset = self.table[row_id][self.parspace.par_names]
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
        for idx, row in enumerate(self.table[self.parspace.par_names]):
            if np.allclose(tuple(parset), tuple(row)):
                mod = self.get_model_from_row(idx)
                break
        else:
            text = f'parset not in all_models table. parset={parset}, ' \
                   f'all_models table: {self.table}'
            self.logger.error(text)
            raise ValueError(text)
        # if parset not in [row[self.parspace.par_names] for row in self.table]:
        #     text = f'parset not in all_models table. parset={parset}, ' \
        #            f'all_models table: {self.table}'
        #     self.logging.error(text)
        #     raise ValueError(text)
        # mod = Model(system=self.system,
        #             settings=self.settings,
        #             parspace=self.parspace,
        #             parset=parset)
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
        mod = Model(system=self.system,
                      settings=self.settings,
                      parspace=self.parspace,
                      parset=parset,
                      directory=self.table['directory'][row_id])
        return mod

    def save(self):
        """Save the all_models table

        """
        self.table.write(self.filename, format='ascii.ecsv', overwrite=True)
        self.logger.debug(f'Model table written to file {self.filename}')


class Model(object):
    """A DYNAMITE model.

    The model can be run by running the methods (i) get_orblib, (ii) get_weights
    and (iii) (in the future) do_orbit_colouring. Running each of these methods
    will return the appropriate object, e.g. model.get_orblib() --> returns an
    OrbitLibrary object model.get_weights(...) --> returns a WeightSolver object

    Parameters
    ----------
    system : dyn.physical_system.System
        Object holding information about the physical system being modelled.
    settings : dyn.config_reader.Settings
        Object holding other settings
    parspace : dyn.parameter_space.ParameterSpace
        A list of parameter objects for this model
    parset : row of an Astropy Table
        contains the values of the potential parameters for this model
    directory : str
        The model directory name (without path). If None or not specified,
        the all_models_file will be searched for the directory name.

    Returns
    -------
    Nothing returned. Attributes holding outputs are are added to the
    object when methods are run.

    """
    def __init__(self,
                 system=None,
                 settings=None,
                 parspace=None,
                 parset=None,
                 directory=None):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.check_parset(parspace, parset)
        self.system = system
        self.settings = settings
        self.parset = parset
        self.parspace = parspace
        # directory of the legacy fortran files
        self.legacy_directory = self.settings.legacy_settings['directory']
        # directory of the input kinematics
        self.in_dir = self.settings.io_settings['input_directory']
        if directory is None:
            self.directory = self.get_model_directory()
        else:
            self.directory = self.settings.io_settings['output_directory'] + \
                             'models/' + directory
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
        directory = self.settings.io_settings['output_directory'] + 'models/'
        models_file = directory + self.settings.io_settings['all_models_file']
        try:
            all_models = ascii.read(models_file)
            self.logger.debug(f'Setting model dir from file {models_file}...')
        except FileNotFoundError:
            sformat = self.system.parameters[0].sformat # this is ml's format
            ml_dir = f"/ml{self.parset['ml']:{sformat}}/"
            directory += f'orblib_000_000{ml_dir}'
            self.logger.info(f'The all_models file {models_file} does not '
                                'exist - setting model '
                                f'directory to {directory}.')
            return directory #######################################
        except:
            self.logger.error('Error reading all_models file. '
                              'Cannot set model directory.')
            raise
        for idx, parset in enumerate(all_models[self.parspace.par_names]):
            if np.allclose(tuple(parset),tuple(self.parset)):
                directory += 'models/' + all_models['directory'][idx]
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
                system=self.system,
                mod_dir=self.directory_noml,
                settings=self.settings.orblib_settings,
                legacy_directory=self.legacy_directory,
                input_directory=self.settings.io_settings['input_directory'],
                parset=self.parset)
        orblib.get_orblib()
        orblib.read_losvd_histograms()
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
        ws_type = self.settings.weight_solver_settings['type']
        if ws_type=='LegacyWeightSolver':
            weight_solver = ws.LegacyWeightSolver(
                    system=self.system,
                    mod_dir=self.directory_noml,
                    settings=self.settings.weight_solver_settings,
                    legacy_directory=self.legacy_directory,
                    ml=self.parset['ml'])
        elif ws_type=='NNLS':
            weight_solver = ws.NNLS(
                    system=self.system,
                    settings=self.settings.weight_solver_settings,
                    directory_with_ml=self.directory)
        else:
            raise ValueError('Unknown WeightSolver type')
        weights, chi2_tot, chi2_kin = weight_solver.solve(orblib)
        self.chi2 = chi2_tot # instrinsic/projected mass + GH coeeficients 1-Ngh
        self.kinchi2 = chi2_kin # GH coeeficients 1-Ngh
        self.weights = weights
        return weight_solver

    def get_ml_of_original_orblib(self):
        """Get ``ml`` of original orblib with shared parameters

        The original ``ml`` is required to rescale orbit libraries for rescaled
        potentials. This method calls the model's orblib's method of the
        same name.

        Returns
        -------
        float
            the original ``ml``

        """
        orblib = dyn_orblib.LegacyOrbitLibrary(
                system=self.system,
                mod_dir=self.directory_noml,
                settings=self.settings.orblib_settings,
                legacy_directory=self.legacy_directory,
                input_directory=self.settings.io_settings['input_directory'],
                parset=self.parset)
        ml_original = orblib.get_ml_of_original_orblib()
        return ml_original

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
                par.value = par.get_raw_value_from_par_value(parset[par_name])
        parspace_copy.validate_parspace()
        self.logger.debug('parset validated against parspace.')

# end
