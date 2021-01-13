import os
import shutil #used to easily copy files
import numpy as np
from astropy import table
from astropy.io import ascii

import weight_solvers as ws
import orblib as dyn_orblib

class AllModels(object):

    def __init__(self,
                 from_file=True,
                 filename='all_models.ecsv',
                 settings=None,
                 parspace=None):
        # self.settings = settings
        outdir = settings.io_settings['output_directory']
        filename = settings.io_settings['all_models_file']
        filename = f'{outdir}{filename}'
        self.filename = filename
        self.parspace = parspace
        if from_file and os.path.isfile(filename):
            print('Previous models have been found:')
            print(f'Reading {filename} into AllModels.table')
            self.read_completed_model_file()
        else:
            print('No previous models have been found:')
            print(f'Making an empty table in AllModels.table')
            self.make_empty_table()

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
        # print(names, dtype)
        self.table = table.Table(names=names, dtype=dtype)
        return

    def read_completed_model_file(self):
        self.table = ascii.read(self.filename)
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
        return

    def get_parset_from_row(self, row_id):
        parset = self.table[row_id][self.parspace.par_names]
        return parset

    def save(self):
        self.table.write(self.filename, format='ascii.ecsv', overwrite=True)


class Model(object):
    '''
    A DYNAMITE model. The Model can be run by running the methods (i)
    get_orblib, (ii) get_weights, (iii) (in the future) do_orbit_colouring.
    Ruuning each of these methods will adds a new attribute to the model, e.g.
    model.get_orblib(...) --> creates an attribute --> model.orblib
    model.get_weights(...) --> creates an attribute --> model.weight_solver
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
        self.system = system
        self.settings = settings
        self.parset = parset
        self.parspace = parspace

    def get_model_directory(self):
        out_dir = self.settings.io_settings['output_directory']
        #out_dir += self.system.name
        out_dir += 'models/'
        # add all parameters to directory name except ml
        for par0, pval0 in zip(self.parspace, self.parset):
            pname0 = par0.name
            psfmt0  = par0.sformat
            if pname0 != 'ml':
                out_dir += pname0
                out_dir += format(pval0, psfmt0)
        # add ml to directory name
        out_dir += '/'
        for par0, pval0 in zip(self.parspace, self.parset):
            pname0 = par0.name
            psfmt0  = par0.sformat
            if par0.name == 'ml':
                out_dir += pname0
                out_dir += format(pval0, psfmt0)
        out_dir += '/'
        # remove all whitespace
        out_dir = out_dir.replace(" ", "")
        return out_dir

    def create_model_directory(self, path):
        if not os.path.exists(path):
            os.makedirs(path)


class LegacySchwarzschildModel(Model):

    def __init__(self,
                 **kwargs):
        super().__init__(**kwargs)
        # directory of the Schwarzschild fortran files
        self.legacy_directory = self.settings.legacy_settings['directory']
        # directory of the input kinematics
        self.in_dir = self.settings.io_settings['input_directory']
        self.directory = self.get_model_directory()
        # TODO: replace following with use string.split() or /../
        self.directory_noml=self.directory[:-7]

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

    def get_weights(self, orblib):
        # create the weight solver object
        weight_solver = ws.LegacyWeightSolver(
                system=self.system,
                mod_dir=self.directory_noml,
                settings=self.settings.weight_solver_settings,
                legacy_directory=self.legacy_directory,
                ml=self.parset['ml'])

        # weight_solver = ws.PrashsCoolNewWeightSolver(
        #     system=self.system,
        #     settings=self.settings.weight_solver_settings)

        # TODO: extract other outputs e.g. orbital weights
        chi2, kinchi2 = weight_solver.solve(orblib)
        # store chi2 to the model
        self.chi2 = chi2
        self.kinchi2 = kinchi2



# end
