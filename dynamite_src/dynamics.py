import os
import glob
import numpy as np
import subprocess


class OrbitLibrary(object):

    def __init__(self,
                 system=None,
                 settings=None):
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
        
        #set the current directory to the directory in which the models are computed
        cur_dir=os.getcwd()
        os.chdir(self.mod_dir)

        #generate the initial conditions for the orbit library 
        self.generate_ics()
        self.integrate_orbits()

        #set the current directory to the dynamite directory
        os.chdir(cur_dir)


    def generate_ics(self):
        
        cmdstr=self.write_executable_for_ics()
        self.execute_ics(cmdstr)

    def write_executable_for_ics(self):
        
        cmdstr = 'cmda' + str(self.system.name) + '_' + str(int(np.random.uniform(0, 1) * 100000.0))
        print(cmdstr)
        
        #create the fortran executable
        txt_file = open(cmdstr, "w")
        
        txt_file.write('#!/bin/bash' + '\n')
        txt_file.write(
            'grep finished datfil/orbstart.dat || ' + self.legacy_directory +'/orbitstart < infil/orbstart.in >> datfil/orbstart.log' + '\n')
        
        txt_file.close()
        
        #returns the name of the executable
        return cmdstr


    def execute_ics(self, cmdstr):

        p1 = subprocess.call('bash ' + cmdstr, shell=True)

    def read_ics(self):
        # ...
        pass

    def integrate_orbits(self):
        
        cmdstr_tube,cmdstr_box =self.write_executable_for_integrate_orbits()
        self.execute_integrate_orbits(cmdstr_tube)
        self.execute_integrate_orbits(cmdstr_box)


    def write_executable_for_integrate_orbits(self):
        
        #tubeorbits
        cmdstr_tube = 'cmdb' + str(self.system.name) + '_' + str(int(np.random.uniform(0, 1) * 100000.0))
        
        txt_file = open(cmdstr_tube, "w")
        
        txt_file.write('#!/bin/bash' + '\n')
        txt_file.write('grep Writing datfil/orblib.dat.tmp && rm -f datfil/orblib.dat.tmp datfil/orblib.dat' + '\n')
        txt_file.write( self.legacy_directory +'/orblib < infil/orblib.in >> datfil/orblib.log' + '\n')
        txt_file.write('rm datfil/mass_qgrid.dat datfil/mass_radmass.dat datfil/mass_aper.dat infil/mass_aper.dat' + '\n')
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
        txt_file.write(
            'grep Writing datfil/orblibbox.dat.tmp && rm -f datfil/orblibbox.dat.tmp datfil/orblibbox.dat' + '\n')
        txt_file.write(self.legacy_directory + '/orblib < infil/orblibbox.in >> datfil/orblibbox.log' + '\n')
        txt_file.write('# if the gzipped orbit library does not exist zip it' + '\n')
        txt_file.write('test -e datfil/orblibbox.dat.bz2 || bzip2 -k datfil/orblibbox.dat' + '\n')
        txt_file.write('rm datfil/orblibbox.dat' + '\n')  
        txt_file.close()


        #returns the name of the executable
        return cmdstr_tube, cmdstr_box 
        

    def execute_integrate_orbits(self, cmdstr):
                
        p2 = subprocess.call('bash ' + cmdstr, shell=True)
        

    def read_integrate_orbits(self):
        # ...
        pass




# end
