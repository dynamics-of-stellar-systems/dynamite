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
#        self.n_par = 0
        self.parameters = None
        self.distMPc = None
        self.name = None
        self.position_angle = None
#        for idx, component in enumerate(args):
        for component in args:
            self.add_component(component)
            # self.n_cmp += 1
            # self.cmp_list += [component]
            # self.n_pot += component.contributes_to_potential
            # self.n_kin += len(component.kinematic_data)
            # self.n_pop += len(component.population_data)
            # self.n_par += len(component.parameters)

    def add_component(self, cmp):
        self.cmp_list += [cmp]
        self.n_cmp += 1
        self.n_pot += cmp.contributes_to_potential
        self.n_kin += len(cmp.kinematic_data)
        self.n_pop += len(cmp.population_data)
#        self.n_par += len(cmp.parameters)

    def validate(self):
        if not(self.ml and self.distMPc and self.galname and self.position_angle):
            raise ValueError('System needs ml, distMPc, galname, and position_angle attributes')

    def __repr__(self):
        return f'{self.__class__.__name__} with {self.__dict__}'


class Component(object):

    def __init__(self,
                 name = None,                     # string
                 visible=None,                    # Boolean
                 contributes_to_potential=None,   # Boolean
                 symmetry=None,                   # OPTIONAL, spherical, axisymm, or triax
                 kinematic_data=[],               # a list of Kinematic objects
                 population_data=[],              # a list of Population objects
                 parameters=[]):                  # a list of Parameter objects
        if name == None:
            self.name = self.__class__.__name__
        else:
            self.name = name
#        self.symmetries = ['spherical', 'axisymm', 'triax']
        self.visible = visible
        self.contributes_to_potential = contributes_to_potential
        self.symmetry = symmetry
        self.kinematic_data = kinematic_data
        self.population_data = population_data
        self.parameters = parameters
        self.validate()

    def validate(self):
        # if self.symmetry not in self.symmetries:
        #     raise ValueError('Illegal symmetry ' + str(self.symmetry) + '. Allowed: ' + str(self.symmetries))
        if self.contributes_to_potential is None:
            raise ValueError(f'Component {self.__class__.__name__} needs contributes_to_potential entry')

    def __repr__(self):
        return (f'\n{self.__class__.__name__}({self.__dict__}\n)')


class VisibleComponent(Component):

    def __init__(self,
                 mge=None,
                 **kwds):
         # visible components have MGE surface density
        self.mge = mge
        super().__init__(visible=True, **kwds)

    # def validate(self):
    #     super(VisibleComponent, self).validate()
    #     check for valid MGE data


class AxisymmetricVisibleComponent(VisibleComponent):

    def __init__(self, **kwds):
        super().__init__(symmetry='axisymm', **kwds)

    # def validate(self):
    #     super(VisibleComponent, self).validate()
    #     check for valid MGE data


class TriaxialVisibleComponent(VisibleComponent):

    def __init__(self, **kwds):
        super().__init__(symmetry='triax', **kwds)

    # def validate(self):
    #     super(VisibleComponent, self).validate()
    #     check for valid MGE data
    
    def triax_pqu2tpp(self,p,q,qobs,u):
        """
        transfer (p, q, u) to the three viewing angles (theta, psi, phi)
        with known flatting qobs.
        Taken from schw_basics
        We should possibly revisit the expressions later
        
        """
        
        p2 = np.double(p) ** 2
        q2 = np.double(q) ** 2
        u2 = np.double(u) ** 2
        o2 = np.double(qobs) ** 2
        
        w1 = (u2 - q2) * (o2 * u2 - q2) / ((1.0 - q2) * (p2 - q2))
        w2 = (u2 - p2) * (p2 - o2 * u2) * (1.0 - q2) / ((1.0 - u2) * (1.0 - o2 * u2) * (p2 - q2))
        w3 = (1.0 - o2 * u2) * (p2 - o2 * u2) * (u2 - q2) / ((1.0 - u2) * (u2 - p2) * (o2 * u2 - q2))
        
        if w1 >=0.0 :
            theta = np.arccos(np.sqrt(w1)) * 180 /np.pi
        else:
            theta=np.nan
                
        if w2 >=0.0 :
            phi = np.arctan(np.sqrt(w2)) * 180 /np.pi
        else:
            phi=np.nan

        if w3 >=0.0 :
            psi = 180 - np.arctan(np.sqrt(w3)) * 180 /np.pi
        else:
            psi=np.nan

        # print("******************************")
        # print('theta, phi, psi')
        # print(theta, phi, psi)
        # print("******************************")

        return theta,psi,phi


class DarkComponent(Component):

    def __init__(self,
                 density=None,
                 **kwds):
        # these have no observed properties (MGE/kinematics/populations)
        # instead they are initialised with an input density function
        self.density = density
        # self.mge = 'self.fit_mge()'
        super().__init__(visible=False,
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
        # self.mge = MGES.intrinsic_MGE_from_xyz_grid(xyz_grid, rho)


class Plummer(DarkComponent):

    def __init__(self, **kwds):
        super().__init__(symmetry='spherical', **kwds)

    def density(x, y, z, pars):
        M, a = pars
        r = (x**2 + y**2 + z**2)**0.5
        rho = 3*M/4/np.pi/a**3 * (1. + (r/a)**2)**-2.5
        return rho


class NFW(DarkComponent):

    def __init__(self, **kwds):
        super().__init__(symmetry='spherical', **kwds)

    def density(x, y, z, pars):
        M, c = pars
        # rho = ...
        # return rho


# end
