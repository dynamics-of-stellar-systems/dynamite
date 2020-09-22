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
                 legacy_directory=None,
                 executor=None):
        self.system = system
        self.mod_dir = mod_dir
        self.settings = settings
        self.legacy_directory = legacy_directory
        self.executor = executor
        #set the current directory to the directory in which the models are computed
        cur_dir=os.getcwd()
        os.chdir(self.mod_dir)
        print("Calculating the orbit library for the proposed potential.")
        #generate the initial conditions for the orbit library
        self.generate_ics()
        self.integrate_orbits()
        #set the current directory to the dynamite directory
        os.chdir(cur_dir)
        print("Orbit integration is finished.")

    def generate_ics(self):
        cmdstr = self.executor.write_executable_for_ics()
        self.executor.execute(cmdstr)

    def read_ics(self):
        # ...
        pass

    def integrate_orbits(self):
        cmdstrs = self.executor.write_executable_for_integrate_orbits()
        cmdstr_tube, cmdstr_box = cmdstrs
        self.executor.execute(cmdstr_tube)
        self.executor.execute(cmdstr_box)

    def read_orbits(self):
        # ...
        pass




# end
