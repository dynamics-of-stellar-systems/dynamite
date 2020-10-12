import numpy as np
import subprocess

class Executor(object):

    def __init__(self,
                 system=None,
                 legacy_directory=None,
                 executor_settings=None):
        """Class to handle different modes of executions: e.g. Local, Slurm, etc
        Any method which varies for different execution mode should go here:
        This is a dummy class - specific implementations are given by child
        classes below.
        """
        self.system = system
        self.legacy_directory = legacy_directory
        self.executor_settings = executor_settings

    def write_executable_for_ics(self):
        """called by OrbitLibrary.generate_ics()

        Returns
        -------
        string
            name of the executable file which has been created

        """
        #
        executable_filename = ''
        return executable_filename

    def write_executable_for_integrate_orbits(self):
        """called by OrbitLibrary.integrate_orbits()

        Returns
        -------
        string, string
            name of the executable files created for box and tube orbits
        """
        #
        executable_filename_tube = ''
        executable_filename_box = ''
        return executable_filename_tube, executable_filename_box

    def write_executable_for_weight_solver(self, ml):
        """called by WeightSolver.solve()

        Returns
        -------
        string
            name of the executable file which has been created

        """
        #
        executable_filename = ''
        return executable_filename

    def execute(self, executable_filename):
        """Execute a given file

        Parameters
        ----------
        executable_filename : string

        Returns
        -------
        None
        """
        return


class Local(Executor):

    def __init__(self,
                 **kwargs):
        super().__init__(**kwargs)

    def write_executable_for_ics(self):
        cmdstr = 'cmda' + str(self.system.name) + '_' + str(int(np.random.uniform(0, 1) * 100000.0))
        #create the fortran executable
        txt_file = open(cmdstr, "w")
        txt_file.write('#!/bin/bash' + '\n')
        txt_file.write(
        #    'grep finished datfil/orbstart.dat || ' + self.legacy_directory +'/orbitstart < infil/orbstart.in >> datfil/orbstart.log' + '\n')
             self.legacy_directory +'/orbitstart < infil/orbstart.in >> datfil/orbstart.log' + '\n')
        txt_file.close()
        #returns the name of the executable
        return cmdstr

    def write_executable_for_integrate_orbits(self):
        #tubeorbits
        cmdstr_tube = 'cmdb' + str(self.system.name) + '_' + str(int(np.random.uniform(0, 1) * 100000.0))
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
        cmdstr_box = 'cmdc' + str(self.system.name) + '_' + str(int(np.random.uniform(0, 1) * 100000.0))
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

    def write_executable_for_weight_solver(self, ml):
        '{:06.2f}'.format(ml)
        nn='ml'+'{:01.2f}'.format(ml)+'/nn'
        cmdstr = 'cmdd' + str(self.system.name) + '_' + str(int(np.random.uniform(0, 1) * 100000.0))
        txt_file = open(cmdstr, "w")
        txt_file.write('#!/bin/bash' + '\n')
        txt_file.write('# if the gzipped orbit library exist unzip it' + '\n')
        txt_file.write('test -e ../datfil/orblib.dat || bunzip2 -k  datfil/orblib.dat.bz2 ' + '\n')
        txt_file.write('test -e ../datfil/orblibbox.dat || bunzip2 -k  datfil/orblibbox.dat.bz2' + '\n')
        txt_file.write('test -e ' + str(nn) + '_kinem.out || ' +
                           self.legacy_directory +'/triaxnnls_CRcut < ' + str(nn) + '.in >>' +str(nn) + 'ls.log' + '\n') #TODO: specify which nnls to use
        txt_file.write('rm datfil/orblib.dat' + '\n')
        txt_file.write('rm datfil/orblibbox.dat' + '\n')
        txt_file.close()
        return cmdstr

    # def execute_ics(self, cmdstr):
    #     p1 = subprocess.call('bash '+cmdstr, shell=True)
    #
    # def execute_integrate_orbits(self, cmdstr):
    #     p2 = subprocess.call('bash '+cmdstr, shell=True)
    #
    # def execute_weight_solver(self, cmdstr):
    #     p4 = subprocess.call('bash '+cmdstr, shell=True)

    def execute(self, cmdstr):
        """As the individual exectute methods are all identical, we can replace
        them with a single execute method
        """
        p = subprocess.call('bash '+cmdstr, shell=True)
        return


class Slurm(Local):
    # !!!!! NOTE !!!!!
    # for the time being I have made Slurm a child of Local
    # so that it inherits all of the methods from Local
    # Sabine - when you implement this, change the above line to:
    #     class Slurm(Executor):
    # !!!!! end note !!!!!

    def __init__(self,
                 **kwargs):
        super().__init__(**kwargs)
        # example of how to access the executor settings:
        # print(self.executor_settings['example_slurm_setting_1'])
        # print(self.executor_settings['example_slurm_setting_2'])




# end
