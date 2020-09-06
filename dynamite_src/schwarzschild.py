import os
import glob
import shutil #used to easily copy files
import numpy as np
from astropy import table
from astropy.io import ascii

import weight_solvers as ws
import dynamics as dyn

class AllModels(object):

    def __init__(self,
                 from_file=True,
                 filename='all_models.ecsv',
                 settings=None,
                 parspace=None):
        self.settings = settings
        outdir = settings.io_settings['output_directory']
        filename = f'{outdir}{filename}'
        self.filename = filename
        print(filename)
        self.parspace = parspace
        if from_file and os.path.isfile(filename):
            print(f'reading {filename} into table attribute')
            self.read_completed_model_file()
        else:
            print(f'making empty table attribute')
            self.make_empty_table()

    def make_empty_table(self):
        names = self.parspace.par_names.copy()
        dtype = [np.float64 for n in names]
        # add the columns from legacy version
        names += ['chi2', 'kinchi2', 'time_modified', 'directory']
        dtype += [np.float64, np.float64, np.float64, '<U100']
        # add extra columns
        names += ['ics_done', 'orblib_done', 'weights_done', 'all_done']
        dtype += [bool, bool, bool, bool]
        # which_iter will record which iteration of parameters a model came from
        names += ['which_iter']
        dtype += [int]
        ncols = len(names)
        data = np.zeros((0, ncols))
        self.table = table.Table(data,
                                 names=names,
                                 dtype=dtype)
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


    def create_model_directory(self,path):
        if not os.path.exists(path):
            os.makedirs(path)




class LegacySchwarzschildModel(SchwarzschildModel):

    def __init__(self,
                 system=None,
                 settings=None,
                 parset=None,
                 parspace=None,
                 execute_run=True):

        self.system = system
        self.settings = settings
        self.parset = parset
        self.parspace = parspace

        #directory of the Schwarzschild fortran files
        self.legacy_directory = self.settings.legacy_settings['directory']

        #directory of the input kinematics
        self.in_dir = self.settings.io_settings['input_directory']

        # In general, arguments for:
        #    - OrbitLibrary = system, settings.orblib_settings
        #    - WeightSolver = orb_lib, system, settings.weight_solver_settings
        # Legacy versions however access these quantities through files which
        # are saved in mod_dir, so we can just pass the directory name instead

        self.directory = self.get_model_directory()

        #might be removed later, use string.split() or /../
        self.directory_noml=self.directory[:-7]

        if execute_run:

            # create self.directory if it doesn't exist
            self.create_model_directory(self.directory)

            #and also the model directories
            self.create_model_directory(self.directory_noml+'infil/')
            self.create_model_directory(self.directory_noml+'datfil/')


            #add communication to the user
            #check if orbit library was calculated already
            if not os.path.isfile(self.directory_noml+'datfil/orblib.dat.bz2') or not os.path.isfile(self.directory_noml+'datfil/orblibbox.dat.bz2') :

                #prepare the fortran input files for orblib, possibly as template files
                self.create_fortran_input_orblib(self.directory_noml+'infil/')

                #FOR NOW: copy the *.dat-files from datafiles, TODO: will be created via prepare_data, links hardcoded not so good
                shutil.copyfile(self.in_dir+'kin_data.dat', self.directory_noml+'infil/kin_data.dat')
                shutil.copyfile(self.in_dir+'aperture.dat', self.directory_noml+'infil/aperture.dat')
                shutil.copyfile(self.in_dir+'bins.dat', self.directory_noml+'infil/bins.dat')



                #calculate orbit libary
                orb_lib = dyn.LegacyOrbitLibrary(
                        system=self.system,
                        mod_dir=self.directory_noml,
                        settings=self.settings.orblib_settings,
                        legacy_directory=self.legacy_directory)


            #prepare fortran input file for nnls
            self.create_fortran_input_nnls(self.directory)

            #apply the nnls
            weight_solver = ws.LegacyWeightSolver(
                    system=self.system,
                    mod_dir=self.directory_noml,
                    settings=self.settings.weight_solver_settings,
                    legacy_directory=self.legacy_directory,
                    ml=self.parset['ml'])

            chi2, kinchi2 = weight_solver.solve()

            # TODO: extract other outputs e.g. orbital weights

            # store result
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
             str(10**self.parset['mass'])+'\n' + \
             str(self.parset['a'])+'\n' + \
             str(self.settings.orblib_settings['nE']) +' ' +str(self.settings.orblib_settings['logrmin']) +' ' +str(self.settings.orblib_settings['logrmax'])+ '\n' + \
             str(self.settings.orblib_settings['nI2']) +'\n' + \
             str(self.settings.orblib_settings['nI3']) +'\n' + \
             str(self.settings.orblib_settings['dithering']) +'\n' + \
             dm_specs +'\n' + \
             str(self.parset['dc']) +' ' + str(self.parset['f']) +'\n'

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
              '"' + stars.kinematic_data[0].aperturefile +'"' +'\n'  + \
              '1                              [use psf 1] ' +'\n'  + \
              self.settings.orblib_settings['hist_vel'] + '  ' + self.settings.orblib_settings['hist_sigma'] + '  ' + self.settings.orblib_settings['hist_bins'] +'   [histogram]' +'\n'  + \
              '1                              [use binning for aperture 1] ' +'\n'  + \
              '"' + stars.kinematic_data[0].binfile +'"' +'\n'  + \
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
              '"' + stars.kinematic_data[0].aperturefile +'"' +'\n'  + \
              '1                              [use psf 1] ' +'\n'  + \
              self.settings.orblib_settings['hist_vel'] + '  ' + self.settings.orblib_settings['hist_sigma'] + '  ' + self.settings.orblib_settings['hist_bins'] +'   [histogram]' +'\n'  + \
              '1                              [use binning for aperture 1] ' +'\n'  + \
              '"' + stars.kinematic_data[0].binfile +'"' +'\n'  + \
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
              '"' + stars.kinematic_data[0].aperturefile +'"' +  '\n'  + \
              str(len(stars.kinematic_data[0].PSF['sigma'])) + '                              [# of gaussians components]'  +'\n' + \
              str(psf_weight) + '   ' + str(psf_sigma) + '                     [weight sigma]' +  '\n'  + \
              '"' + stars.kinematic_data[0].binfile +'"' +'\n'  + \
              '"datfil/mass_aper.dat"'

        triaxmassbin_file= open(path+'triaxmassbin.in',"w")
        triaxmassbin_file.write(text)
        triaxmassbin_file.close()


    def create_fortran_input_nnls(self,path):

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
        str(self.settings.weight_solver_settings['ml_scale_factor']) + '                                [ scale factor related to M/L, sqrt( (M/L)_k / (M/L)_ref ) ]' + '\n' + \
        'datfil/orblib.dat' +'\n' + \
        'datfil/orblibbox.dat' +'\n' + \
        str(self.settings.weight_solver_settings['nnls_solver']) + '                                  [ nnls solver ]'

        nn_file= open(path+'nn.in',"w")
        nn_file.write(text)
        nn_file.close()



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
