import os
# import shutil #used to easily copy files
import numpy as np
from astropy import table
from astropy.io import ascii
import logging

import weight_solvers as ws
import orblib as dyn_orblib

class AllModels(object):

    def __init__(self,
                 system=None,
                 from_file=True,
                 filename='all_models.ecsv',
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
        outdir = self.settings.io_settings['output_directory']
        filename = f'{outdir}{filename}'
        self.filename = filename

    def make_empty_table(self):
        names = self.parspace.par_names.copy()
        dtype = [np.float64 for n in names]
        # add the columns from legacy version
        names += ['chi2', 'kinchi2', 'time_modified']
        dtype += [np.float64, np.float64, '<M8[ms]']
        # add extra columns
        names += ['orblib_done', 'weights_done', 'all_done']
        dtype += [bool, bool, bool]
        # which_iter will record which iteration of parameters a model came from
        names += ['which_iter']
        dtype += [int]
        # ncols = len(names)
        # data = np.zeros((0, ncols))
        # self.table = table.Table(data,
        #                          names=names,
        #                          dtype=dtype)
        # self.logger.debug(names, dtype)
        self.table = table.Table(names=names, dtype=dtype)
        return

    def read_completed_model_file(self):
        self.table = ascii.read(self.filename)
        self.logger.debug(f'Models read from file {self.filename}')
        return

    def read_legacy_chi2_file(self, legacy_filename):
        """
        Taken from schw_basics.py
        reads in Chi2 of all the models finished, ./griddata/_chi2.cat
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
        return

    def get_parset_from_row(self, row_id):
        parset = self.table[row_id][self.parspace.par_names]
        return parset

    def get_model_from_parset(self, parset):
        if parset not in [row[self.parspace.par_names] for row in self.table]:
            text = f'parset not in all_models table. parset={parset}, ' \
                   f'all_models table: {self.table}'
            self.logging.error(text)
            raise ValueError(text)
        mod = Model(system=self.system,
                    settings=self.settings,
                    parspace=self.parspace,
                    parset=parset)
        return mod

    def get_model_from_row(self, row_id):
        parset0 = self.get_parset_from_row(row_id)
        mod = self.get_model_from_parset(parset0)
        # mod0 = Model(system=self.system,
        #              settings=self.settings,
        #              parspace=self.parspace,
        #              parset=parset0)
        return mod

    def save(self):
        self.table.write(self.filename, format='ascii.ecsv', overwrite=True)
        self.logger.debug(f'Model table written to file {self.filename}')


class Model(object):
    '''
    A DYNAMITE model.

    The model can be run by running the methods (i) get_orblib, (ii) get_weights
    and (iii) (in the future) do_orbit_colouring. Running each of these methods
    will return the appropriate object, e.g. model.get_orblib() --> returns an
    OrbitLibrary object model.get_weights(...) --> returns a WeightSolver object
    '''
    def __init__(self,
                 system=None,
                 settings=None,
                 parspace=None,
                 parset=None):
        """
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

        Returns
        -------
        Nothing returned. Attributes holidng outputs are are added to the object
        when methods are run.

        """
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.system = system
        self.settings = settings
        self.parset = parset
        self.parspace = parspace
        # directory of the Schwarzschild fortran files
        self.legacy_directory = self.settings.legacy_settings['directory']
        # directory of the input kinematics
        self.in_dir = self.settings.io_settings['input_directory']
        self.directory = self.get_model_directory()
        self.logger.debug(f'Model directory: {self.directory}')
        self.directory_noml=self.directory[:self.directory[:-1].rindex('/')+1]
        self.logger.debug(f'Model directory up to ml: {self.directory_noml}')

    def get_model_directory(self):
        out_dir = self.settings.io_settings['output_directory']
        #out_dir += self.system.name
        out_dir += 'models/'
        # add all parameters to directory name except ml
        for par0, pval0 in zip(self.parspace, self.parset):
            pname0 = par0.name
            psfmt0  = par0.sformat
            if pname0 != 'ml':
                out_dir += f'{pname0}'
                out_dir += format(pval0, psfmt0)+'_'
        out_dir = out_dir[:-1]
        # add ml to directory name
        out_dir += '/'
        for par0, pval0 in zip(self.parspace, self.parset):
            pname0 = par0.name
            psfmt0  = par0.sformat
            if par0.name == 'ml':
                out_dir += f'{pname0}'
                out_dir += format(pval0, psfmt0)
        out_dir += '/'
        # remove all whitespace
        out_dir = out_dir.replace(" ", "")
        return out_dir

    def create_model_directory(self, path):
        if not os.path.exists(path):
            os.makedirs(path)
            self.logger.debug(f'Created directory {path}')
        else:
            self.logger.debug(f'Using existing directory {path}')

    def setup_directories(self):
        # create self.directory if it doesn't exist
        self.create_model_directory(self.directory)
        # and also the model directories
        self.create_model_directory(self.directory_noml+'infil/')
        self.create_model_directory(self.directory_noml+'datfil/')

    def get_orblib(self):
        # make orbit libary object
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
        # create the weight solver object
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



# end
