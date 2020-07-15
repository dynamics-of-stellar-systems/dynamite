# classes to hold the physical components of the system
# e.g. the stellar light, dark matter, black hole, globular clusters

import numpy as np

class System(object):

    def __init__(self, *args):
        self.n_cmp = 0
        self.cmp_list = []
        self.n_pot = 0
        self.n_kin = 0
        self.n_pop = 0
        self.n_par = 0
        for idx, component in enumerate(args):
            self.n_cmp += 1
            self.cmp_list += [component]
            self.n_pot += component.contributes_to_potential
            self.n_kin += len(component.kinematic_data)
            self.n_pop += len(component.population_data)
            self.n_par += len(component.parameters)

    def add_component(self, cmp):
        self.cmp_list += [cmp]
        self.n_cmp += 1
        self.npot += cmp.contributes_to_potential
        self.n_kin += len(cmp.kinematic_data)
        self.n_pop += len(cmp.population_data)
        self.n_par += len(cmp.parameters)
    

class Component(object):

    def __init__(self,
                 visible=None,                    # Boolean
                 contributes_to_potential=None,   # Boolean
                 symmetry=None,                   # spherical, axisymm, or triax
                 kinematic_data=[],             # a list of Kinematic objects
                 population_data=[],            # a list of Population objects
                 parameters=[]):                # a list of Parameter objects
        self.visible = visible
        self.contributes_to_potential = contributes_to_potential
        self.symmetry = symmetry
        self.kinematic_data = kinematic_data
        self.population_data = population_data
        self.parameters = parameters


class VisibleComponent(Component):

    def __init__(self,
                 mge=None,
                 **kwds):
         # visible components are MGE surface density
        self.mge = mge
        super(VisibleComponent, self).__init__(visible=True,
                                               **kwds)


class DarkComponent(Component):

    def __init__(self,
                 density=None,
                 **kwds):
        # these have no observed properties (MGE/kinematics/populations)
        # instead they are initialised with an input density function
        self.density = density
        self.mge = 'self.fit_mge()'
        super(DarkComponent, self).__init__(visible=False,
#                                            contributes_to_potential=True,
                                            kinematic_data=[],
                                            population_data=[],
                                            **kwds)

    def fit_mge(self,
                density,
                parameters,
                xyz_grid=[]):
        # fit an MGE for a given set of parameters
        # will be used in potential calculation
        rho = self.density.evaluate(xyz_grid, parameters)
#        self.mge = MGES.intrinsic_MGE_from_xyz_grid(xyz_grid, rho)


class Plummer(DarkComponent):

    def __init__(self,
                 parameter_M,
                 parameter_a,
                 **kwds):
        super(Plummer, self).__init__(symmetry='spherical',
                                      parameters=[parameter_M,
                                                  parameter_a],
                                      **kwds)

    def density(x, y, z, pars):
        M, a = pars
        r = (x**2 + y**2 + z**2)**0.5
        rho = 3*M/4/np.pi/a**3 * (1. + (r/a)**2)**-2.5
        return rho


class NFW(DarkComponent):

    def __init__(self,
                 parameter_M,
                 parameter_c,
                 **kwds):
        super(NFW, self).__init__(symmetry='spherical',
                                  parameters=[parameter_M,
                                              parameter_c],
                                  **kwds)

    def density(x, y, z, pars):
        M, c = pars
        # rho = ...
        # return rho


# end
