import os
import numpy as np
import subprocess
from scipy.io import FortranFile

import sys
this_dir = os.path.dirname(__file__)
if not this_dir in sys.path:
    sys.path.append(this_dir)
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
                 legacy_directory=None):
        self.system = system
        self.mod_dir = mod_dir
        self.settings = settings
        self.legacy_directory = legacy_directory

    def get_orbit_ics(self):
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        print("Calculating the initial conditions for this potential")
        cmdstr = self.write_executable_for_ics()
        p = subprocess.call('bash '+cmdstr, shell=True)
        os.chdir(cur_dir)
        print("Orbit integration is finished.")

    def write_executable_for_ics(self):
        cmdstr = 'cmd_orb_start'
        #create the fortran executable
        txt_file = open(cmdstr, "w")
        txt_file.write('#!/bin/bash' + '\n')
        txt_file.write(
        #    'grep finished datfil/orbstart.dat || ' + self.legacy_directory +'/orbitstart < infil/orbstart.in >> datfil/orbstart.log' + '\n')
             self.legacy_directory +'/orbitstart < infil/orbstart.in >> datfil/orbstart.log' + '\n')
        txt_file.close()
        #returns the name of the executable
        return cmdstr

    def read_ics(self):
        # ...
        pass

    def get_orbit_library(self):
        # move to model directory
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        print("Calculating the orbit library for the proposed potential.")
        cmdstr_tube, cmdstr_box = self.write_executable_for_integrate_orbits()
        p = subprocess.call('bash '+cmdstr_tube, shell=True)
        p = subprocess.call('bash '+cmdstr_box, shell=True)
        print("Orbit integration is finished.")
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
        txt_file.write( self.legacy_directory + '/triaxmass       < infil/triaxmass.in ' + '\n')
        txt_file.write( self.legacy_directory + '/triaxmassbin    < infil/triaxmassbin.in ' + '\n')
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
        subprocess.call(['bunzip2', '-k', f'datfil/{fileroot}.dat.bz2'])
        orblibf = FortranFile(f'datfil/{fileroot}.dat', 'r')
        # read size of orbit library
        # from integrator_setup_write, lines 506 - 5129:
        tmp = orblibf.read_ints(np.int32)
        norb, t2, t3, t4, ndith = tmp
        # from qgrid_setup_write, lines 2339-1350:
        tmp = orblibf.read_ints(np.int32)
        size_ql_1, size_qph_minus1, size_qth_minus1, size_qlr_minus1 = tmp
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
        for j in range(norb):
            t1,t2,t3,t4,t5 = orblibf.read_ints(np.int32)
            orbtypes[j, :] = orblibf.read_ints(np.int32)
            quad_light = orblibf.read_reals(float)
            for k in range(nconstr):
                ivmin, ivmax = orblibf.read_ints(np.int32)
                if (ivmin <= ivmax):
                    tmp = orblibf.read_reals(float)
                    velhist[j, ivmin+nvhist:ivmax+nvhist+1, k] = tmp
        subprocess.call(['rm', f'datfil/{fileroot}.dat'])
        os.chdir(cur_dir)
        vedg_pos = np.arange(1, nbins_vhist+1, 2) * dvhist/2.
        vedg_neg = -vedg_pos[::-1]
        vedg = np.concatenate((vedg_neg, vedg_pos))
        velhist = dyn_kin.Histogram(xedg=vedg,
                                    y=velhist,
                                    normalise=False)
        return velhist

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
        error_msg = 'velocity array must be symmetric'
        assert np.all(orblib.xedg == -orblib.xedg[::-1]), error_msg

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
        error_msg = 'orblibs have different number of velocity bins'
        assert n_vel_bins1==n_vel_bins2, error_msg
        error_msg = 'orblibs have different velocity arrays'
        assert np.array_equal(orblib1.x, orblib2.x), error_msg
        error_msg = 'orblibs have different number of spatial bins'
        assert n_spatial_bins1==n_spatial_bins2, error_msg
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
        Reads the LOSVD histograms: reads box orbits and non-box orbits, flips
        the latter, and combines. Sets the 'losvd_histograms' attribute which
        is a Histogram of the combined orbit libraries. (I thnk) this ordering
        is compatible with the weights read in by LegacyWeightSolver.read_weights
        '''
        tube_orblib = self.read_orbit_base('orblib')
        tube_orblib = self.duplicate_flip_and_interlace_orblib(tube_orblib)
        box_orblib = self.read_orbit_base('orblibbox')
        orblib = self.combine_orblibs(tube_orblib, box_orblib)
        self.losvd_histograms = orblib

# end
