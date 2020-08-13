class OrbitLibrary(object):

    def __init__(self,
                 system=None,
                 settings=None):
        self.system = system
        self.settings = settings
        self.generate_ics()
        self.integrate_loop(timesteps)

    def generate_ics(self):
        pass

    def integrate_orbits(self):
        pass


class LegacyOrbitLibrary(OrbitLibrary):

    def __init__(self,
                 mod_dir=None,
                 settings=None):
        self.mod_dir = mod_dir
        self.settings = settings
        self.generate_ics()
        self.integrate_orbits()

    def generate_ics(self):
        # self.write_executable_for_ics()
        # self.execute_ics()
        pass

    def write_executable_for_ics(self):
        # ...
        pass

    def execute_ics(self):
        # ...
        pass

    def read_ics(self):
        # ...
        pass

    def integrate_orbits(self):
        # self.write_executable_for_integrate_orbits()
        # self.execute_integrate_orbits()
        pass

    def write_executable_for_integrate_orbits(self):
        # ...
        pass

    def execute_integrate_orbits(self):
        # ...
        pass

    def read_integrate_orbits(self):
        # ...
        pass




# end
