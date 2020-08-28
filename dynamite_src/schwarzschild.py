import os
import numpy as np
from astropy import table
from astropy.io import ascii

from . import weight_solvers as ws
from . import dynamics as dyn


class AllModels(object):

    def __init__(self,
                 from_file=True,
                 filename='all_models.ecsv',
                 settings=None,
                 parspace=None,
                 *args,
                 **kwargs):
        self.settings = settings
        outdir = settings.output_settings['directory']
        filename = f'{outdir}{filename}'
        self.filename = filename
        print(filename)
        if from_file and os.path.isfile(filename):
            print(f'reading {filename} into table attribute')
            self.read_completed_model_file(*args, **kwargs)
        else:
            print(f'making empty table attribute')
            self.make_empty_table(parspace, *args, **kwargs)

    def make_empty_table(self, parspace, *args, **kwargs):
        names = parspace.par_names.copy()
        dtype = [np.float64 for n in names]
        # add the columns from legacy version
        names += ['chi2', 'kinchi2', 'time_modified', 'directory']
        dtype += [np.float64, np.float64, np.float64, '<U100']
        # add extra columns
        names += ['ics_done', 'orblib_done', 'weights_done']
        dtype += [bool, bool, bool]
        ncols = len(names)
        data = np.zeros((0, ncols))
        self.table = table.Table(data,
                                 names=names,
                                 dtype=dtype,
                                 *args, **kwargs)
        return

    def read_completed_model_file(self, *args, **kwargs):
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
        mods.write(self.filename, format='ascii.ecsv')
        return


class SchwarzschildModel(object):

    def __init__(self,
                 system=None,
                 settings=None,
                 parset=None,
                 parspace=None):
        self.system = system
        self.settings = settings
        self.parset = parset
        self.parspace = parspace

        # orb_lib = dyn.OrbitLibrary(
        #     system=system,
        #     settings=settings.orblib_settings)

        # weight_solver = ... instantiate the selected solver
        # weight_solver.set_orb_lib(orb_lib)
        # orb_wts, chi2 = weight_solver.solve()
        #
        # # do colouring
        # orb_labels = []
        # for cmp in system.cmp_list:
        #     for pop_data0 in cmp.population_data:
        #         orb_labels0 = pop_data0.colouring_recipe(orb_lib,
        #                                                   orb_wts)
        #         orb_labels += [orb_labels0]
        #
        # # store resuult
        # self.parset = parset
        # self.orb_lib = orb_lib
        # self.orb_wts = orb_wts
        # self.chi2 = chi2
        # self.orb_labels = orb_labels

    def get_model_directory(self):
        out_dir = self.settings.output_settings['directory']
        out_dir += self.system.name
        out_dir += '/'
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


class LegacySchwarzschildModel(SchwarzschildModel):

    def __init__(self,
                 system=None,
                 settings=None,
                 parset=None,
                 parspace=None):

        self.system = system
        self.settings = settings
        self.parset = parset
        self.parspace = parspace

        # In general, arguments for:
        #    - OrbitLibrary = system, settings.orblib_settings
        #    - WeightSolver = orb_lib, system, settings.weight_solver_settings
        # Legacy versions however access these quantities through files which
        # are saved in mod_dir, so we can just pass the directory name instead

        self.directory = self.get_model_directory()

        # TODO:
        # here we can make self.directory if it doesn't exist, and fill it with
        # all the necessary schwpy/fortran input files

        orb_lib = dyn.LegacyOrbitLibrary(
            mod_dir=self.directory,
            settings=settings.orblib_settings)

        weight_solver = ws.LegacyWeightSolver(
            mod_dir=self.directory,
            settings=settings.weight_solver_settings)

        chi2, kinchi2 = weight_solver.solve()
        # TODO: extract other outputs e.g. orbital weights

        # store result
        self.chi2 = chi2
        self.kinchi2 = kinchi2


class SchwarzschildModelLoop(object):

    def __init__(self,
                 system=None,
                 parset_list=None,
                 weight_solver=None,
                 settings=None
                 ):
        pass
        self.models = []
        # for parset0 in parset_list:
        #     if settings.legacy_mode:
        #         mod0 = LegacySchwarzschildModel(
        #             system=system,
        #             settings=settings
        #             parset=parset0,
        #             parspace=parspace)
        #     else:
        #         mod0 = SchwarzschildModel(
        #             system=system,
        #             settings=settings,
        #             parset=parset0,
        #             parspace=parspace)
        #     self.models += [mod0]








# end
