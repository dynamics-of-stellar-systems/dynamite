import os
import glob
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
                 executor=None,
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
        executor : dyn.executor.Executor
            Handles differences between execting models on your local machines
            vs submission on clusters, HPC modes etc
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
        self.executor = executor

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
                executor=self.executor)
        self.orblib = orblib
        # check if orbit library was calculated already
        check1 = os.path.isfile(self.directory_noml+'datfil/orblib.dat.bz2')
        check2 = os.path.isfile(self.directory_noml+'datfil/orblibbox.dat.bz2')
        if not check1 or not check2:
            # prepare the fortran input files for orblib
            self.create_fortran_input_orblib(self.directory_noml+'infil/')
            kinematics = self.system.cmp_list[2].kinematic_data[0]
            old_filename = self.directory_noml+'infil/kin_data.dat'
            kinematics.convert_to_old_format(old_filename)
            aperture_file = self.in_dir + kinematics.aperturefile
            shutil.copyfile(aperture_file,
                            self.directory_noml+'infil/aperture.dat')
            binfile = self.in_dir + kinematics.binfile
            shutil.copyfile(binfile,
                            self.directory_noml+'infil/bins.dat')
            # calculate orbit libary
            self.orblib.get_orbit_library()

    def get_weights(self):
        # prepare fortran input file for nnls
        self.create_fortran_input_nnls(self.directory_noml)
        # create the weight solver object
        self.weight_solver = ws.LegacyWeightSolver(
                system=self.system,
                mod_dir=self.directory_noml,
                settings=self.settings.weight_solver_settings,
                legacy_directory=self.legacy_directory,
                ml=self.parset['ml'],
                executor=self.executor)
        # TODO: extract other outputs e.g. orbital weights
        chi2, kinchi2 = self.weight_solver.solve()
        # store chi2 to the model
        self.chi2 = chi2
        self.kinchi2 = kinchi2

    def create_fortran_input_orblib(self,path):

        #-------------------
        #write parameters.in
        #-------------------

        stars=self.system.cmp_list[2]

        #used to derive the viewing angles
        q=self.parset['q']
        p=self.parset['p']
        u=self.parset['u']

        #the minimal flattening from stellar mge
        qobs=np.amin(stars.mge.data['q'])

        #TODO: add softening length somewhere
        # r_BH='1d-03'
        # Done! This is now paramater "a" of the Plummer

        #TODO: which dark matter profile
        dm_specs='1 2'

        theta,psi,phi=self.system.cmp_list[2].triax_pqu2tpp(p,q,qobs,u)

        #header
        len_mge=len(stars.mge.data)
        #footer (#double check the order of theta, phi, psi) and dm properties
        text=str(self.system.distMPc)+'\n'+ \
             '{:06.9f}'.format(theta)+' '+ '{:06.9f}'.format(phi)+' '+ '{:06.9f}'.format(psi) + '\n' + \
             str(self.parset['ml'])+'\n' + \
             str(self.parset['mass'])+'\n' + \
             str(self.parset['a'])+'\n' + \
             str(self.settings.orblib_settings['nE']) +' ' +str(self.settings.orblib_settings['logrmin']) +' ' +str(self.settings.orblib_settings['logrmax'])+ '\n' + \
             str(self.settings.orblib_settings['nI2']) +'\n' + \
             str(self.settings.orblib_settings['nI3']) +'\n' + \
             str(self.settings.orblib_settings['dithering']) +'\n' + \
             dm_specs +'\n' + \
             str(self.parset['dc']) +' ' + str(self.parset['f'])

        #parameters.in
        np.savetxt(path+'parameters.in',stars.mge.data,header=str(len_mge),footer=text,comments='',fmt=['%10.2f','%10.5f','%10.5f','%10.2f'])

        #parmsb.in (assumed to be the same as paramters.in)
        np.savetxt(path+'paramsb.in',stars.mge.data,header=str(len_mge),footer=text,comments='',fmt=['%10.2f','%10.5f','%10.5f','%10.2f'])

        #-------------------
        #write orbstart.in
        #-------------------

        text='infil/parameters.in' +'\n' + \
        'datfil/orbstart.dat' +'\n' + \
        'datfil/begin.dat' +'\n' + \
        'datfil/beginbox.dat'

        orbstart_file= open(path+'orbstart.in',"w")
        orbstart_file.write(text)
        orbstart_file.close()

        #-------------------
        #write orblib.in
        #-------------------

        i=0
        psf_weight=(stars.kinematic_data[0].PSF['weight'])[i]
        psf_sigma=(stars.kinematic_data[0].PSF['sigma'])[i]
        n_psf=[[1]]   #len(stars.kinematic_data) #needs to be revised

        #TODO:needs to be slightly changed for more psfs, loop

        text1='#counterrotation_setupfile_version_1' +'\n' + \
            'infil/parameters.in' +'\n' + \
            'datfil/begin.dat' +'\n' + \
            str(self.settings.orblib_settings['orbital_periods']) + '                            [orbital periods to intergrate orbits]' +'\n' + \
            str(self.settings.orblib_settings['sampling']) + '                          [points to sample for each orbit in the merid. plane]' +'\n' + \
            str(self.settings.orblib_settings['starting_orbit']) + '                              [starting orbit]' +'\n' + \
            str(self.settings.orblib_settings['number_orbits']) + '                             [orbits  to intergrate; -1 --> all orbits]' +'\n' + \
            str(self.settings.orblib_settings['accuracy']) + '                         [accuracy]' +'\n' + \
            str(len(stars.kinematic_data)) + '                              [number of psfs of the kinematic cubes]' +'\n'

        psf= str(len(stars.kinematic_data[0].PSF['sigma'])) + '                              [# of gaussians components]'  +'\n' + \
             str(psf_weight) + '   ' + str(psf_sigma) + '                    [weight, sigma]' +  '\n'


        text2=str(len(stars.kinematic_data)) + '                              [apertures]' +'\n'  + \
              '"infil/' + stars.kinematic_data[0].aperturefile +'"' +'\n'  + \
              '1                              [use psf 1] ' +'\n'  + \
              self.settings.orblib_settings['hist_vel'] + '  ' + self.settings.orblib_settings['hist_sigma'] + '  ' + self.settings.orblib_settings['hist_bins'] +'   [histogram]' +'\n'  + \
              '1                              [use binning for aperture 1] ' +'\n'  + \
              '"infil/' + stars.kinematic_data[0].binfile +'"' +'\n'  + \
              'datfil/orblib.dat '

        orblib_file= open(path+'orblib.in',"w")
        orblib_file.write(text1)
        orblib_file.write(psf)
        orblib_file.write(text2)
        orblib_file.close()

        #-------------------
        #write orblibbox.in
        #-------------------

        #TODO:why not paramsb.in?
        text1='#counterrotation_setupfile_version_1' +'\n' + \
            'infil/parameters.in' +'\n' + \
            'datfil/beginbox.dat' +'\n' + \
            str(self.settings.orblib_settings['orbital_periods']) + '                            [orbital periods to intergrate orbits]' +'\n' + \
            str(self.settings.orblib_settings['sampling']) + '                          [points to sample for each orbit in the merid. plane]' +'\n' + \
            str(self.settings.orblib_settings['starting_orbit']) + '                              [starting orbit]' +'\n' + \
            str(self.settings.orblib_settings['number_orbits']) + '                             [orbits  to intergrate; -1 --> all orbits]' +'\n' + \
            str(self.settings.orblib_settings['accuracy']) + '                         [accuracy]' +'\n' + \
            str(len(stars.kinematic_data)) + '                              [number of psfs of the kinematic cubes]' +'\n'

        text2=str(len(stars.kinematic_data)) + '                              [apertures]' +'\n'  + \
              '"infil/' + stars.kinematic_data[0].aperturefile +'"' +'\n'  + \
              '1                              [use psf 1] ' +'\n'  + \
              self.settings.orblib_settings['hist_vel'] + '  ' + self.settings.orblib_settings['hist_sigma'] + '  ' + self.settings.orblib_settings['hist_bins'] +'   [histogram]' +'\n'  + \
              '1                              [use binning for aperture 1] ' +'\n'  + \
              '"infil/' + stars.kinematic_data[0].binfile +'"' +'\n'  + \
              'datfil/orblibbox.dat '

        orblibbox_file= open(path+'orblibbox.in',"w")
        orblibbox_file.write(text1)
        orblibbox_file.write(psf) #this is the same as for orblib.in
        orblibbox_file.write(text2)
        orblibbox_file.close()

        #-------------------
        #write triaxmass.in
        #-------------------

        text='infil/paramsb.in' +'\n' + \
        'datfil/orblib.dat' +'\n' + \
        'datfil/mass_radmass.dat' +'\n' + \
        'datfil/mass_qgrid.dat'

        triaxmass_file= open(path+'triaxmass.in',"w")
        triaxmass_file.write(text)
        triaxmass_file.close()

        #-------------------
        #write triaxmassbin.in
        #-------------------

        text='infil/paramsb.in' +'\n' + \
              str(int(np.max(n_psf))) + '                              [# of apertures]'  +'\n'  + \
              '"infil/' + stars.kinematic_data[0].aperturefile +'"' + '\n' + \
              str(len(stars.kinematic_data[0].PSF['sigma'])) + '                              [# of gaussians components]'  +'\n' + \
              str(psf_weight) + '   ' + str(psf_sigma) + '                     [weight sigma]' +  '\n'  + \
              '"infil/' + stars.kinematic_data[0].binfile +'"' +'\n'  + \
              '"datfil/mass_aper.dat"'

        triaxmassbin_file= open(path+'triaxmassbin.in',"w")
        triaxmassbin_file.write(text)
        triaxmassbin_file.close()

    def create_fortran_input_nnls(self,path):

        #for the ml the model is only scaled. We therefore need to know what is the ml that was used for the orbit library
        infile=path+'infil/parameters.in'
        lines = [line.rstrip('\n').split() for line in open(infile)]
        ml_orblib=float((lines[-9])[0])

        #-------------------
        #write nn.in
        #-------------------

        text='infil/parameters.in' +'\n' + \
        str(self.settings.weight_solver_settings['regularisation'])   + '                                  [ regularization strength, 0 = no regularization ]' +'\n'  + \
        'ml'+ '{:01.2f}'.format(self.parset['ml']) + '/nn' +'\n' + \
        'datfil/mass_qgrid.dat' +'\n' + \
        'datfil/mass_aper.dat' +'\n' + \
        str(self.settings.weight_solver_settings['number_GH']) + '	                           [ # of GH moments to constrain the model]' +'\n' + \
        'infil/kin_data.dat' +'\n' + \
        str(self.settings.weight_solver_settings['GH_sys_err']) + '    [ systemic error of v, sigma, h3, h4... ]' + '\n' + \
        str(self.settings.weight_solver_settings['lum_intr_rel_err']) + '                               [ relative error for intrinsic luminosity ]' +'\n' + \
        str(self.settings.weight_solver_settings['sb_proj_rel_err']) + '                               [ relative error for projected SB ]' + '\n' + \
        str(np.sqrt(self.parset['ml']/ml_orblib))  + '                                [ scale factor related to M/L, sqrt( (M/L)_k / (M/L)_ref ) ]' + '\n' + \
        'datfil/orblib.dat' +'\n' + \
        'datfil/orblibbox.dat' +'\n' + \
        str(self.settings.weight_solver_settings['nnls_solver']) + '                                  [ nnls solver ]'

        nn_file= open(path+'ml'+'{:01.2f}'.format(self.parset['ml'])+'/nn.in',"w")
        nn_file.write(text)
        nn_file.close()


















# end
