import os
import numpy as np
import subprocess
import shutil
from scipy.io import FortranFile
from astropy.table import Table, vstack
import logging

import sys
this_dir = os.path.dirname(__file__)
if not this_dir in sys.path:
    sys.path.append(this_dir)
import physical_system as physys
import kinematics as dyn_kin

class OrbitLibrary(object):

    def __init__(self,
                 system=None,
                 settings=None):
        """A class for orbit libraries.

        Parameters
        ----------
        system : type
            Description of parameter `system`.
        settings : type
            Description of parameter `settings`.

        Returns
        -------
        type
            Description of returned object.

        """
        self.system = system
        self.settings = settings
        self.generate_ics()
        self.integrate_loop(timesteps)

    def generate_ics(self):
        #ics: initial conditions
        pass

    def integrate_orbits(self):
        pass


class LegacyOrbitLibrary(OrbitLibrary):

    def __init__(self,
                 system=None,
                 mod_dir=None,
                 settings=None,
                 legacy_directory=None,
                 input_directory=None,
                 parset=None):

        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')

        self.system = system
        self.mod_dir = mod_dir
        self.settings = settings
        self.legacy_directory = legacy_directory
        self.in_dir = input_directory
        self.parset = parset

    def get_orblib(self):
        # check if orbit library was calculated already
        check1 = os.path.isfile(self.mod_dir+'datfil/orblib.dat.bz2')
        check2 = os.path.isfile(self.mod_dir+'datfil/orblibbox.dat.bz2')
        if not check1 or not check2:
            # prepare the fortran input files for orblib
            self.create_fortran_input_orblib(self.mod_dir+'infil/')
            stars = self.system.get_component_from_class( \
                                            physys.TriaxialVisibleComponent)
            kinematics = stars.kinematic_data
            #create the kinematic input files for each kinematic dataset
            for i in np.arange(len(kinematics)):
                if len(kinematics)==1:
                    old_filename = self.mod_dir+'infil/kin_data.dat'
                else:
                    old_filename = self.mod_dir+'infil/kin_data_'+str(i)+'.dat'
                kinematics[i].convert_to_old_format(old_filename)

                aperture_file = self.in_dir + kinematics[i].aperturefile
                shutil.copyfile(aperture_file,
                            self.mod_dir+'infil/'+ kinematics[i].aperturefile)
                binfile = self.in_dir + kinematics[i].binfile
                shutil.copyfile(binfile,
                            self.mod_dir+'infil/'+ kinematics[i].binfile)

            #combined kinematics for legacy dynamite
            if len(kinematics)>1:
                kinematics_combined=kinematics[0]

                kinematics_combined.data=vstack((kinematics[0].data, kinematics[1].data))
                old_filename = self.mod_dir+'infil/kin_data_combined.dat'
                kinematics_combined.convert_to_old_format(old_filename)

            # calculate orbit libary
            self.get_orbit_ics()
            self.get_orbit_library()

    def create_fortran_input_orblib(self, path):
        #---------------------------------------------
        #write parameters_pot.in and parameters_lum.in
        #---------------------------------------------
        stars = \
          self.system.get_component_from_class(physys.TriaxialVisibleComponent)
        bh = self.system.get_component_from_class(physys.Plummer)
        dh = self.system.get_component_from_class(physys.NFW)
        # used to derive the viewing angles
        q=self.parset[f'q-{stars.name}']
        p=self.parset[f'p-{stars.name}']
        u=self.parset[f'u-{stars.name}']
        # the minimal flattening from stellar mge
        # qobs=np.amin(stars.mge.data['q'])
        # TODO: which dark matter profile
        dm_specs='1 2'
        theta, psi, phi = stars.triax_pqu2tpp(p,q,u)
        # header
        len_mge=len(stars.mge_pot.data) # must be the same as for mge_lum
        # footer (#double check the order of theta, phi, psi) and dm properties
        text=str(self.system.distMPc)+'\n'+ \
             '{:06.9f}'.format(theta)+' '+ '{:06.9f}'.format(phi)+' '+ '{:06.9f}'.format(psi) + '\n' + \
             str(self.parset['ml'])+'\n' + \
             str(self.parset[f'm-{bh.name}'])+'\n' + \
             str(self.parset[f'a-{bh.name}'])+'\n' + \
             str(self.settings['nE']) +' ' +str(self.settings['logrmin']) +' ' +str(self.settings['logrmax'])+ '\n' + \
             str(self.settings['nI2']) +'\n' + \
             str(self.settings['nI3']) +'\n' + \
             str(self.settings['dithering']) +'\n' + \
             dm_specs +'\n' + \
             str(self.parset[f'c-{dh.name}']) +' ' + str(self.parset[f'f-{dh.name}'])

        #parameters_pot.in
        np.savetxt(path+'parameters_pot.in',stars.mge_pot.data,header=str(len_mge),footer=text,comments='',fmt=['%10.2f','%10.5f','%10.5f','%10.2f'])

        #parameters_lum.in
        np.savetxt(path+'parameters_lum.in',stars.mge_lum.data,header=str(len_mge),footer=text,comments='',fmt=['%10.2f','%10.5f','%10.5f','%10.2f'])

        #-------------------
        #write orbstart.in
        #-------------------

        text = f"{self.settings['random_seed']}\n"
        text += 'infil/parameters_pot.in' +'\n' + \
        'datfil/orbstart.dat' +'\n' + \
        'datfil/begin.dat' +'\n' + \
        'datfil/beginbox.dat'

        orbstart_file= open(path+'orbstart.in',"w")
        orbstart_file.write(text)
        orbstart_file.close()

        #-------------------
        #write orblib.in
        #-------------------

        n_psf=len(stars.kinematic_data)

        text0 = f"{self.settings['random_seed']}\n"

        text1='#counterrotation_setupfile_version_1' +'\n' + \
            'infil/parameters_pot.in' +'\n' + \
            'datfil/begin.dat' +'\n' + \
            str(self.settings['orbital_periods']) + '                            [orbital periods to intergrate orbits]' +'\n' + \
            str(self.settings['sampling']) + '                          [points to sample for each orbit in the merid. plane]' +'\n' + \
            str(self.settings['starting_orbit']) + '                              [starting orbit]' +'\n' + \
            str(self.settings['number_orbits']) + '                             [orbits  to intergrate; -1 --> all orbits]' +'\n' + \
            str(self.settings['accuracy']) + '                         [accuracy]' +'\n' + \
            str(n_psf) + '                              [number of psfs of the kinematic data]' +'\n'

        psf=''

        for i in np.arange(n_psf):
            psf= psf+str(len(stars.kinematic_data[i].PSF['sigma'])) + '                              [# of gaussians components for psf '+str(i+1)+']'  +'\n'

        #loop over psf for each kinematic set. A psf can be described by multiple Gaussians

        for i in np.arange(n_psf):
            n_psf_comp=len((stars.kinematic_data[i].PSF['weight']))
            for k in np.arange(n_psf_comp):
                psf_weight=(stars.kinematic_data[i].PSF['weight'])[k]
                psf_sigma=(stars.kinematic_data[i].PSF['sigma'])[k]
                psf=psf+(str(psf_weight) + '   ' + str(psf_sigma) + '                     [weight, sigma]' +  '\n')


        apertures=str(n_psf) + '                              [# of apertures]' +'\n'


        for i in np.arange(n_psf):
            apertures+='"infil/' + stars.kinematic_data[i].aperturefile +'"\n' + \
            str(i+1) + ' '*30 + \
            f'[use psf {i+1} {stars.kinematic_data[i].name}] \n'
        for i in np.arange(n_psf):
            # apertures+= self.settings['hist_width'] + '  ' + self.settings['hist_center'] + '  ' + self.settings['hist_bins'] +'             [histogram]' +'\n'
            apertures+= stars.kinematic_data[i].hist_width + '  ' + \
                        stars.kinematic_data[i].hist_center + '  ' + \
                        stars.kinematic_data[i].hist_bins + \
                        f'             [histogram {stars.kinematic_data[i].name}]\n'

        for i in np.arange(n_psf):
            apertures+=f'1                              [use binning for aperture {1+i}]\n'

        for i in np.arange(n_psf):
              apertures+='"infil/' + stars.kinematic_data[i].binfile + f'"          [binning for aperture {i+1}]\n'

        text2 = 'datfil/orblib.dat'

        orblib_file= open(path+'orblib.in',"w")
        orblib_file.write(text0)
        orblib_file.write(text1)
        orblib_file.write(psf)
        orblib_file.write(apertures)
        orblib_file.write(text2)
        orblib_file.close()

        #-------------------
        #write orblibbox.in
        #-------------------

        text0 = f"{self.settings['random_seed']}\n"
        text1='#counterrotation_setupfile_version_1' +'\n' + \
            'infil/parameters_pot.in' +'\n' + \
            'datfil/beginbox.dat' +'\n' + \
            str(self.settings['orbital_periods']) + '                            [orbital periods to intergrate orbits]' +'\n' + \
            str(self.settings['sampling']) + '                          [points to sample for each orbit in the merid. plane]' +'\n' + \
            str(self.settings['starting_orbit']) + '                              [starting orbit]' +'\n' + \
            str(self.settings['number_orbits']) + '                             [orbits  to intergrate; -1 --> all orbits]' +'\n' + \
            str(self.settings['accuracy']) + '                         [accuracy]' +'\n' + \
            str(len(stars.kinematic_data)) + '                              [number of psfs of the kinematic cubes]' +'\n'

        text2='datfil/orblibbox.dat'

        orblibbox_file= open(path+'orblibbox.in',"w")
        orblibbox_file.write(text0)
        orblibbox_file.write(text1)
        orblibbox_file.write(psf) #this is the same as for orblib.in
        orblibbox_file.write(apertures)
        orblibbox_file.write(text2)
        orblibbox_file.close()

        #-------------------
        #write triaxmass.in
        #-------------------

        text='infil/parameters_lum.in' +'\n' + \
        'datfil/orblib.dat' +'\n' + \
        'datfil/mass_radmass.dat' +'\n' + \
        'datfil/mass_qgrid.dat'

        triaxmass_file= open(path+'triaxmass.in',"w")
        triaxmass_file.write(text)
        triaxmass_file.close()

        #-------------------
        #write triaxmassbin.in
        #-------------------

        text='infil/parameters_lum.in' +'\n' + \
              str(int(np.max(n_psf))) + '                              [# of apertures]'  +'\n'
        apertures=''
        for i in np.arange(n_psf):
            psf_weight=(stars.kinematic_data[i].PSF['weight'])[k]
            psf_sigma=(stars.kinematic_data[i].PSF['sigma'])[k]
            
            apertures+= \
              '"infil/' + stars.kinematic_data[i].aperturefile +'"' + '\n' + \
              str(len(stars.kinematic_data[i].PSF['sigma'])) + '                              [# of gaussians components]'  +'\n' + \
              str(psf_weight) + '   ' + str(psf_sigma) + '                     [weight sigma]' +  '\n'  + \
              '"infil/' + stars.kinematic_data[i].binfile +'"' +'\n'

        text2='"datfil/mass_aper.dat"'

        triaxmassbin_file= open(path+'triaxmassbin.in',"w")
        triaxmassbin_file.write(text)
        triaxmassbin_file.write(apertures)
        triaxmassbin_file.write(text2)
        triaxmassbin_file.close()

    def get_orbit_ics(self):
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        cmdstr = self.write_executable_for_ics()
        self.logger.info('Calculating initial conditions')
        p = subprocess.call('bash '+cmdstr, shell=True)
        self.logger.debug('...done. ' + \
                          f'Logfile: {self.mod_dir}datfil/orbstart.log')
        os.chdir(cur_dir)

    def write_executable_for_ics(self):
        cmdstr = 'cmd_orb_start'
        #create the fortran executable
        txt_file = open(cmdstr, "w")
        txt_file.write('#!/bin/bash' + '\n')
        txt_file.write(
        #    'grep finished datfil/orbstart.dat || ' + self.legacy_directory +'/orbitstart < infil/orbstart.in >> datfil/orbstart.log' + '\n')
             self.legacy_directory +'/orbitstart < infil/orbstart.in 2>&1 >> datfil/orbstart.log' + '\n')
        txt_file.close()
        #returns the name of the executable
        return cmdstr

    def get_orbit_library(self):
        # move to model directory
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        cmdstr_tube, cmdstr_box = self.write_executable_for_integrate_orbits()
        self.logger.info('Integrating orbit library tube orbits')
        p = subprocess.call('bash '+cmdstr_tube, shell=True)
        self.logger.debug('...done. ' + \
                          f'Logfiles: {self.mod_dir}datfil/orblib.log, ' + \
                          f'{self.mod_dir}datfil/triaxmass.log, ' + \
                          f'{self.mod_dir}datfil/triaxmassbin.log')
        self.logger.info(f'Integrating orbit library box orbits')
        p = subprocess.call('bash '+cmdstr_box, shell=True)
        self.logger.debug('...done. ' + \
                          f'Logfile: {self.mod_dir}datfil/orblibbox.log')
        # move back to original directory
        os.chdir(cur_dir)

    def write_executable_for_integrate_orbits(self):
        #tubeorbits
        cmdstr_tube = 'cmd_tube_orbs'
        txt_file = open(cmdstr_tube, "w")
        txt_file.write('#!/bin/bash' + '\n')
        #txt_file.write('grep Writing datfil/orblib.dat.tmp && rm -f datfil/orblib.dat.tmp datfil/orblib.dat' + '\n')
        txt_file.write('touch datfil/orblib.dat.tmp datfil/orblib.dat' + '\n')
        txt_file.write('rm -f datfil/orblib.dat.tmp datfil/orblib.dat' + '\n')
        txt_file.write( self.legacy_directory +'/orblib < infil/orblib.in >> datfil/orblib.log' + '\n')
        txt_file.write('touch datfil/mass_qgrid.dat datfil/mass_radmass.dat datfil/mass_aper.dat' + '\n')
        txt_file.write('rm datfil/mass_qgrid.dat datfil/mass_radmass.dat datfil/mass_aper.dat' + '\n')
        txt_file.write( self.legacy_directory + '/triaxmass       < infil/triaxmass.in >> datfil/triaxmass.log' + '\n')
        txt_file.write( self.legacy_directory + '/triaxmassbin    < infil/triaxmassbin.in >> datfil/triaxmassbin.log' + '\n')
        txt_file.write('# if the gzipped orbit library does not exist zip it' + '\n')
        txt_file.write('test -e datfil/orblib.dat.bz2 || bzip2 -k datfil/orblib.dat' + '\n')
        txt_file.write('rm datfil/orblib.dat' + '\n')
        txt_file.close()
        #boxorbits
        cmdstr_box = 'cmd_box_orbs'
        txt_file = open(cmdstr_box, "w")
        txt_file.write('#!/bin/bash' + '\n')
        #txt_file.write(
        #    'grep Writing datfil/orblibbox.dat.tmp && rm -f datfil/orblibbox.dat.tmp datfil/orblibbox.dat' + '\n')
        txt_file.write('touch datfil/orblibbox.dat.tmp datfil/orblibbox.dat' + '\n')
        txt_file.write('rm -f datfil/orblibbox.dat.tmp datfil/orblibbox.dat' + '\n')
        txt_file.write(self.legacy_directory + '/orblib < infil/orblibbox.in >> datfil/orblibbox.log' + '\n')
        txt_file.write('# if the gzipped orbit library does not exist zip it' + '\n')
        txt_file.write('test -e datfil/orblibbox.dat.bz2 || bzip2 -k datfil/orblibbox.dat' + '\n')
        txt_file.write('rm datfil/orblibbox.dat' + '\n')
        txt_file.close()
        #returns the name of the executables
        return cmdstr_tube, cmdstr_box

    def read_ics(self):
        # ...
        pass

    def read_orbit_base(self, fileroot):
        """Read a zipped Fortran orbit library from the file
            datfil/{fileroot}.dat.bz2'
        relative to the model directory.

        Parameters
        ----------
        fileroot : string
            this will probably be either 'orblib' or 'orblibbox'

        Returns
        -------
        Histogram
            the orbit library stored in a Histogram object

        """
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        # unzip orblib to a temproary file with ml value attached
        # ml value is needed to prevent different processes clashing
        ml = self.parset['ml']
        aaa = subprocess.run(['bunzip2', '-c', f'datfil/{fileroot}.dat.bz2'],
                             stdout=subprocess.PIPE)
        tmpfname = f'datfil/{fileroot}_{ml}.dat'
        tmpf = open(tmpfname, "wb")
        tmpf.write(aaa.stdout)
        tmpf.close()
        # read the fortran file
        orblibf = FortranFile(tmpfname, 'r')
        # remove temproary file
        subprocess.call(['rm', tmpfname])
        # read size of orbit library
        # from integrator_setup_write, lines 506 - 5129:
        tmp = orblibf.read_ints(np.int32)
        norb, t2, t3, t4, ndith = tmp
        # from qgrid_setup_write, lines 2339-1350:
        quad_light_grid_sizes = orblibf.read_ints(np.int32)
        size_ql, size_qph, size_qth, size_qlr = quad_light_grid_sizes
        quad_lr = orblibf.read_reals(float)
        quad_lth = orblibf.read_reals(float)
        quad_lph = orblibf.read_reals(float)
        # from histogram_setup_write, lines 1917-1926:
        tmp = orblibf.read_record(np.int32, np.int32, float)
        nconstr = tmp[0][0]
        nvhist = tmp[1][0]
        dvhist = tmp[2][0]
        # Next read the histograms themselves.
        orbtypes = np.zeros((norb, ndith**3), dtype=int)
        nbins_vhist = 2*nvhist + 1
        velhist = np.zeros((norb, nbins_vhist, nconstr))
        density_3D = np.zeros((norb, size_qlr, size_qth, size_qph))
        for j in range(norb):
            t1,t2,t3,t4,t5 = orblibf.read_ints(np.int32)
            orbtypes[j, :] = orblibf.read_ints(np.int32)
            quad_light = orblibf.read_reals(float)
            quad_light = np.reshape(quad_light, quad_light_grid_sizes[::-1])
            # quad_light stores orbit features in 3D (r,th,phi) bins.
            # Quad_light[ir,it,ip,XXX] stores
            # - the zeroth moment i.e. density for XXX=0,
            # - the first moments x,y,z,vx,vy,vz for XXX=slice(1,7)
            # - 2nd moments vx^2,vy^2,vz^2,vx*vy,vy*vz,vz*vx for XXX=slice(7,13)
            # - an averaged orbit classification for XXX=slice(13,None)
            # in the bin indexed by (ir,it,ip).
            # We need to extract 3D density for use in weight solving.
            density_3D[j] = quad_light[:,:,:,0]
            for k in range(nconstr):
                ivmin, ivmax = orblibf.read_ints(np.int32)
                if (ivmin <= ivmax):
                    tmp = orblibf.read_reals(float)
                    velhist[j, ivmin+nvhist:ivmax+nvhist+1, k] = tmp
        orblibf.close()
        os.chdir(cur_dir)
        vedg_pos = np.arange(1, nbins_vhist+1, 2) * dvhist/2.
        vedg_neg = -vedg_pos[::-1]
        vedg = np.concatenate((vedg_neg, vedg_pos))
        velhist = dyn_kin.Histogram(xedg=vedg,
                                    y=velhist,
                                    normalise=False)
        return velhist, density_3D

    def duplicate_flip_and_interlace_orblib(self, orblib):
        """ Take an orbit library, create a duplicate library with the velocity
        signs flipped, then interlace the two i.e. so that resulting library
        alternates between flipped/unflipped.

        This creates an orbit library consistent with the Fortran output,
        enforcing the ordering created by the for loops in lines 157-178 of
        triaxnnls_CRcut.f90

        Parameters
        ----------
        orblib : Histogram

        Returns
        -------
        Histogram
            the duplicated, flipped and interlaced orblib

        """
        self.logger.debug('Checking for symmetric velocity array...')
        error_msg = 'velocity array must be symmetric'
        assert np.all(orblib.xedg == -orblib.xedg[::-1]), error_msg
        self.logger.debug('...check ok.')
        losvd = orblib.y
        n_orbs, n_vel_bins, n_spatial_bins = losvd.shape
        reveresed_losvd = losvd[:, ::-1, :]
        new_losvd = np.zeros((2*n_orbs, n_vel_bins, n_spatial_bins))
        new_losvd[0::2] = losvd
        new_losvd[1::2, :] = reveresed_losvd
        new_orblib = dyn_kin.Histogram(xedg=orblib.xedg,
                                       y=new_losvd,
                                       normalise=False)
        return new_orblib

    def combine_orblibs(self, orblib1, orblib2):
        """Combine two histogrammed orbit libraries into one.

        Parameters
        ----------
        orblib1 : Histogram
        orblib2 : Histogram

        Returns
        -------
        Histogram
            the combined orbit libraries

        """
        # check orblibs are compatible
        n_orbs1, n_vel_bins1, n_spatial_bins1 = orblib1.y.shape
        n_orbs2, n_vel_bins2, n_spatial_bins2 = orblib2.y.shape
        self.logger.debug('Checking number of velocity bins...')
        error_msg = 'orblibs have different number of velocity bins'
        assert n_vel_bins1==n_vel_bins2, error_msg
        self.logger.debug('Checking velocity arrays...')
        error_msg = 'orblibs have different velocity arrays'
        assert np.array_equal(orblib1.x, orblib2.x), error_msg
        self.logger.debug('Checking number of spatial bins...')
        error_msg = 'orblibs have different number of spatial bins'
        assert n_spatial_bins1==n_spatial_bins2, error_msg
        self.logger.debug('...checks ok.')
        new_losvd = np.zeros((n_orbs1 + n_orbs2,
                              n_vel_bins1,
                              n_spatial_bins1))
        new_losvd[:n_orbs1] = orblib1.y
        new_losvd[n_orbs1:] = orblib2.y
        new_orblib = dyn_kin.Histogram(xedg=orblib1.xedg,
                                       y=new_losvd,
                                       normalise=False)
        return new_orblib

    def read_losvd_histograms(self):
        '''
        Reads the LOSVD histograms: reads box orbits and tube orbits, flips
        the latter, and combines. Sets the 'losvd_histograms' attribute which
        is a Histogram of the combined orbit libraries.
        '''
        # TODO: check if this ordering is compatible with weights read in by
        # LegacyWeightSolver.read_weights
        tube_orblib, tube_density_3D = self.read_orbit_base('orblib')
        tube_orblib = self.duplicate_flip_and_interlace_orblib(tube_orblib)
        # duplicate and interlace tube_density_3D array
        tube_density_3D = np.repeat(tube_density_3D, 2, axis=0)
        box_orblib, box_density_3D = self.read_orbit_base('orblibbox')
        orblib = self.combine_orblibs(tube_orblib, box_orblib)
        # combine the two density_3D arrays
        density_3D = np.vstack((tube_density_3D, box_density_3D))
        self.losvd_histograms = orblib
        ml_current = self.parset['ml']
        ml_original = self.get_ml_of_original_orblib()
        scale_factor = np.sqrt(ml_current/ml_original)
        self.losvd_histograms.scale_x_values(scale_factor)
        self.intrinsic_masses = density_3D
        self.n_orbs = self.losvd_histograms.y.shape[0]
        self.projected_masses = np.sum(self.losvd_histograms.y, 1)

    def get_ml_of_original_orblib(self):
        infile = self.mod_dir + 'infil/parameters_pot.in'
        lines = [line.rstrip('\n').split() for line in open(infile)]
        ml_original = float((lines[-9])[0])
        return ml_original

# end
