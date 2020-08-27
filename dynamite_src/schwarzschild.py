import os
import glob
import numpy as np
from astropy import table
from astropy.io import ascii

from . import weight_solvers as ws
from . import dynamics as dyn


class AllModels(object):

    def __init__(self,
                 from_file=True,
                 filename='all_models.ecsv',
                 config=None,
                 parspace=None,
                 *args,
                 **kwargs):
        outdir = config.output_settings['directory']
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
                 config=None,
                 parset=None,
                 parspace=None):
        self.system = system
        self.config = config
        self.parset = parset
        self.parspace = parspace

        # orb_lib = dyn.OrbitLibrary(
        #     system=system,
        #     settings=config.orblib_settings)

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
        out_dir = self.config.output_settings['directory']
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
    
    
    def create_model_directory(self,path):
        if not os.path.exists(path):
            os.makedirs(path)
            
        


class LegacySchwarzschildModel(SchwarzschildModel):

    def __init__(self,
                 system=None,
                 config=None,
                 parset=None,
                 parspace=None):

        self.system = system
        self.config = config
        self.parset = parset
        self.parspace = parspace

        #we might add this to config_legacy_settings
        self.legacy_directory = '/Users/sabine/Triaxschwarzschild/triaxschwarzschild'

        # In general, arguments for:
        #    - OrbitLibrary = system, config.orblib_settings
        #    - WeightSolver = orb_lib, system, config.weight_solver_settings
        # Legacy versions however access these quantities through files which
        # are saved in mod_dir, so we can just pass the directory name instead

        self.directory = self.get_model_directory()
        
        #might be removed later, use string.split() or /../
        self.directory_noml=self.directory[:-7]

        # create self.directory if it doesn't exist
        self.create_model_directory(self.directory)
        
        #and also the model directories
        self.create_model_directory(self.directory_noml+'infil/')
        self.create_model_directory(self.directory_noml+'datfil/')        
        
        #and fill it with all the necessary schwpy/fortran input files, possibly as template files
        self.create_fortran_input(self.directory_noml+'infil/')

        ##orb_lib = dyn.LegacyOrbitLibrary(
        ##    mod_dir=self.directory,
        ##    settings=config.orblib_settings)

        ##weight_solver = ws.LegacyWeightSolver(
        ##    mod_dir=self.directory,
        ##    settings=config.weight_solver_settings)

        ##chi2, kinchi2 = weight_solver.solve()
        # TODO: extract other outputs e.g. orbital weights

        # store result
        ##self.chi2 = chi2
        ##self.kinchi2 = kinchi2


    def create_fortran_input(self,path):
        
        
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
        r_BH='1d-03'
        
        #TODO: which dark matter profile
        dm_specs='1 2'
        
        #TODO:
        #Dm might be calculated similar to schadjustinput
        
        theta,psi,phi=self.system.cmp_list[2].triax_pqu2tpp(p,q,qobs,u)
        
        #header
        len_mge=len(stars.mge.data)   
        #footer (#double check the order of theta, phi, psi) and dm properties
        text=str(self.system.distMPc)+'\n'+ \
             '{:06.9f}'.format(theta)+' '+ '{:06.9f}'.format(phi)+' '+ '{:06.9f}'.format(psi) + '\n' + \
             str(self.parset['ml'])+'\n' + \
             str(10**self.parset['bh'])+'\n' + \
             r_BH                           +'\n' + \
             str(self.config.orblib_settings['nE']) +' ' +str(self.config.orblib_settings['logrmin']) +' ' +str(self.config.orblib_settings['logrmax'])+ '\n' + \
             str(self.config.orblib_settings['nI2']) +'\n' + \
             str(self.config.orblib_settings['nI3']) +'\n' + \
             str(self.config.orblib_settings['dithering']) +'\n' + \
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
        
        #TODO?read kinematic psf from kinematics file
        n_psf,psf_weight,psf_sigma=np.loadtxt('datafiles/kinpsffile.dat',skiprows=2)
        
        #TODO the psfs need to be defined in an array before
        
        #DODO:needs to be slightly changed for more psfs
        #TODO aperture and bins.dat need to be clearly defined, 
        text1='#counterrotation_setupfile_version_1' +'\n' + \
            'infil/parameters.in' +'\n' + \
            'datfil/begin.dat' +'\n' + \
            str(self.config.orblib_settings['orbital_periods']) + '                            [orbital periods to intergrate orbits]' +'\n' + \
            str(self.config.orblib_settings['sampling']) + '                          [points to sample for each orbit in the merid. plane]' +'\n' + \
            str(self.config.orblib_settings['starting_orbit']) + '                              [starting orbit]' +'\n' + \
            str(self.config.orblib_settings['number_orbits']) + '                             [orbits  to intergrate; -1 --> all orbits]' +'\n' + \
            str(self.config.orblib_settings['accuracy']) + '                         [accuracy]' +'\n' + \
            str(int(np.max(n_psf))) + '                              [number of psfs of the kinemtic cubes]' +'\n' 

        
        psf= str(len(n_psf[n_psf==1])) + '                              [gaussians components]'  +'\n' + \
             str(psf_weight) + '   ' + str(psf_sigma) + '                    [weight sigma]' +  '\n' 
             
             
        text2=str(int(np.max(n_psf))) + '                              [apertures]' +'\n'  + \
              '"infil/aperture.dat"' +'\n'  + \
              '1                              [use psf 1] ' +'\n'  + \
              self.config.orblib_settings['hist_vel'] + '  ' + self.config.orblib_settings['hist_sigma'] + '  ' + self.config.orblib_settings['hist_bins'] +'   [histogram]' +'\n'  + \
              '1                              [use binning for aperture 1] ' +'\n'  + \
              '"infil/bins.dat"  ' +'\n'  + \
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
            str(self.config.orblib_settings['orbital_periods']) + '                            [orbital periods to intergrate orbits]' +'\n' + \
            str(self.config.orblib_settings['sampling']) + '                          [points to sample for each orbit in the merid. plane]' +'\n' + \
            str(self.config.orblib_settings['starting_orbit']) + '                              [starting orbit]' +'\n' + \
            str(self.config.orblib_settings['number_orbits']) + '                             [orbits  to intergrate; -1 --> all orbits]' +'\n' + \
            str(self.config.orblib_settings['accuracy']) + '                         [accuracy]' +'\n' + \
            str(int(np.max(n_psf))) + '                              [number of psfs of the kinemtic cubes]' +'\n' 

             
             
        text2=str(int(np.max(n_psf))) + '                              [apertures]' +'\n'  + \
              '"infil/aperture.dat"' +'\n'  + \
              '1                              [use psf 1] ' +'\n'  + \
              self.config.orblib_settings['hist_vel'] + '  ' + self.config.orblib_settings['hist_sigma'] + '  ' + self.config.orblib_settings['hist_bins'] +'   [histogram]' +'\n'  + \
              '1                              [use binning for aperture 1] ' +'\n'  + \
              '"infil/bins.dat"  ' +'\n'  + \
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
              '"infil/aperture.dat"' +  '\n'  + \
              str(len(n_psf[n_psf==1])) + '                              [# of gaussians components]'  +'\n' + \
              str(psf_weight) + '   ' + str(psf_sigma) + '                     [weight sigma]' +  '\n'  + \
              '"infil/bins.dat"  ' +'\n'  + \
              '"datfil/mass_aper.dat"' 
        
        triaxmassbin_file= open(path+'triaxmassbin.in',"w")
        triaxmassbin_file.write(text)
        triaxmassbin_file.close()
        
        #-------------------        
        #write nn.in
        #-------------------
        
        text='infil/parameters.in' +'\n' + \
        str(self.config.weight_solver_settings['regularization'])   + '                                  [ regularization strength, 0 = no regularization ]' +'\n'  + \
        'ml'+ str(self.parset['ml']) + '/nn' +'\n' + \
        'datfil/mass_qgrid.dat' +'\n' + \
        'datfil/mass_aper.dat' +'\n' + \
        str(self.config.weight_solver_settings['number_GH']) + '	                           [ # of GH moments to constrain the model]' +'\n' + \
        'infil/kin_data.dat' +'\n' + \
        str(self.config.weight_solver_settings['GH_sys_err']) + '    [ systemic error of v, sigma, h3, h4... ]' + '\n' + \
        str(self.config.weight_solver_settings['lum_intr_rel_err']) + '                               [ relative error for intrinsic luminosity ]' +'\n' + \
        str(self.config.weight_solver_settings['sb_proj_rel_err']) + '                               [ relative error for projected SB ]' + '\n' + \
        str(self.config.weight_solver_settings['ml_scale_factor']) + '                                [ scale factor related to M/L, sqrt( (M/L)_k / (M/L)_ref ) ]' + '\n' + \
        'datfil/orblib.dat' +'\n' + \
        'datfil/orblibbox.dat' +'\n' + \
        str(self.config.weight_solver_settings['nnls_solver']) + '                                  [ nnls solver ]'
        
        nn_file= open(path+'nn.in',"w")
        nn_file.write(text)
        nn_file.close()

        

class SchwarzschildModelLoop(object):

    def __init__(self,
                 system=None,
                 parset_list=None,
                 weight_solver=None,
                 config=None
                 ):
        pass
        # self.models = []
        # for parset0 in parset_list:
        #     if config.legacy_mode:
        #         mod0 = LegacySchwarzschildModel(
        #             system=system,
        #             parset=parset0,
        #             weight_solver=weight_solver,
        #             orblib_pars = orblib_pars)
        #     else:
        #         mod0 = SchwarzschildModel(
        #             system=system,
        #             parset=parset0,
        #             weight_solver=weight_solver,
        #             orblib_pars = orblib_pars)
        #     self.models += [mod0]








# end
