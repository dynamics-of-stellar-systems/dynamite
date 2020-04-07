class Potential(object):

    def __init__(self,
                 system=None,
                 parset=None):
        self.potential_mge = []
        # for component in system.cmp_list:
        #     if component.contributes_to_potential:
        #         self.potential_mge += [component.mge]

    def interpolate_accelaration_field(self):
        # ...
        pass

    def integrate_orbit(self, ics, timesteps):
        # code to integrate orbit given self.ics in potential
        # self.ts = timesteps
        pass


class OrbitLibrary(object):

    def __init__(self,
                 potential=Potential([]),
                 nE=10,
                 nI2=5,
                 nI3=5,
                 n_dither=3,
                 timesteps=[]
                 ):
        self.potential = potential
        self.nE = nE
        self.nI2 = nI2
        self.nI3 = nI3
        self.n_dither = n_dither
        self.n_orbs = nE * nI2 * nI3 * n_dither**3      # ... I think?
        self.timesteps = timesteps
        self.generate_ics()
        self.integrate_loop(timesteps)

    def generate_ics(self):
        # code to generate ICs in self.potential
        # self.ics = ics
        pass

    def integrate_loop(self,
                       timesteps):
        # orbits = np.zeros(size=(self.n_orbs, len(timesteps), 6))
        # for idx, ic0 in enumerate(self.ics):
        #     orbit0 = self.potential.integrate(ic0, timesteps)
        #     orbits[idx, :, :] = orbit0
        # self.orbits = orbits
        pass


# end
