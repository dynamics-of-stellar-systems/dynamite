# classes to hold the physical components of the system
# e.g. the stellar light, dark matter, black hole, globular clusters

import numpy as np

# some tricks to add the current path to sys.path (so the imports below work)

import os.path
import sys

this_dir = os.path.dirname(__file__)
if not this_dir in sys.path:
    sys.path.append(this_dir)

import mges as mge


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
        """
        Ensures the System has the required attributes, at least one component,
        and the ml parameter.
        Additionally, the sformat string for the ml parameter is set.

        Raises
        ------
        ValueError : if required attributes or components are missing, or if
                     there is no ml parameter

        Returns
        -------
        None.

        """
        if not(self.distMPc and self.name and self.position_angle):
            raise ValueError('System needs distMPc, name, '
                             'and position_angle attributes')
        if not self.cmp_list:
            raise ValueError('System has no components')
        if len(self.parameters) != 1 and self.parameters[0].name != 'ml':
            raise ValueError('System needs ml as its sole parameter')
        self.parameters[0].update(sformat = '6.2f')

    def __repr__(self):
        return f'{self.__class__.__name__} with {self.__dict__}'

    def get_component_from_name(self, cmp_name):
        cmp_list_list = np.array([cmp0.name for cmp0 in self.cmp_list])
        idx = np.where(cmp_list_list == cmp_name)
        error_msg = f"There should be 1 and only 1 component named {cmp_name}"
        assert len(idx[0]) == 1, error_msg
        component = self.cmp_list[idx[0][0]]
        return component

class Component(object):

    def __init__(self,
                 name = None,                    # string
                 visible=None,                   # Boolean
                 contributes_to_potential=None,  # Boolean
                 symmetry=None,       # OPTIONAL, spherical, axisymm, or triax
                 kinematic_data=[],              # a list of Kinematic objects
                 population_data=[],             # a list of Population objects
                 parameters=[]):                 # a list of Parameter objects
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
        # self.validate()

    def validate(self, par_format=None):
        """
        Ensures the Component has the required attributes and parameters.
        Additionally, the sformat strings for the parameters are set.

        Parameters
        ----------
        par_format : a dict with parameter_name:sformat pairs. Mandatory.

        Raises
        ------
        ValueError : if a required attribute is missing or the required
                     parameters do not exist

        Returns
        -------
        None.

        """
        # if self.symmetry not in self.symmetries:
        #     raise ValueError('Illegal symmetry ' + str(self.symmetry) + \
        #                      '. Allowed: ' + str(self.symmetries))
        par = par_format.keys()
        errstr = f'Component {self.__class__.__name__} needs attribute '
        if self.visible is None:
            raise ValueError(errstr + 'visible')
        if self.contributes_to_potential is None:
            raise ValueError(errstr + 'contributes_to_potential')
        if not self.parameters:
            raise ValueError(errstr + 'parameters')

        # if len(self.parameters) != len(par):
        #     raise ValueError(f'{self.__class__.__name__} needs exactly '
        #         f'{len(par)} paramater(s), not {len(self.parameters)}')
        if set([p.name for p in self.parameters]) != set(par):
            raise ValueError(f'{self.__class__.__name__} needs parameters '
                             f'{par}, not {[p.name for p in self.parameters]}')
        self.set_format(par_format)

    def set_format(self, par_format=None):
        if par_format is None:
            raise ValueError(f'{self.__class__.__name__}: no format string')
        for p in self.parameters:
            p.update(sformat = par_format[p.name])

    def __repr__(self):
        return (f'\n{self.__class__.__name__}({self.__dict__}\n)')


class VisibleComponent(Component):

    def __init__(self,
                 mge=None,
                 **kwds):
         # visible components have MGE surface density
        self.mge = mge
        super().__init__(visible=True, **kwds)

    def validate(self, **kwds):
        super().validate(**kwds)
        if not isinstance(self.mge, mge.MGE):
            raise ValueError(f'{self.__class__.__name__}.mge must be '
                             'mges.MGE object')


class AxisymmetricVisibleComponent(VisibleComponent):

    def __init__(self, **kwds):
        super().__init__(symmetry='axisymm', **kwds)

    def validate(self):
        par_format = {'par1':'6.2f', 'par2':'6.3f'}
        super().validate(par_format)


class TriaxialVisibleComponent(VisibleComponent):

    def __init__(self, **kwds):
        super().__init__(symmetry='triax', **kwds)

    def validate(self):
        par_format = {'q':'6.2f', 'p':'6.3f', 'u':'7.4f'}
        super().validate(par_format=par_format)

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

    def validate(self):
        par_format = {'mass':'6.2f', 'a':'7.2f'}
        super().validate(par_format)
        # if len(self.parameters) != 2:
        #     raise ValueError(f'{self.__class__.__name__} needs exactly 2 '
        #                      f'paramaters, not {len(self.parameters)}')


class NFW(DarkComponent):

    def __init__(self, **kwds):
        self.legacy_code = 1
        super().__init__(symmetry='spherical', **kwds)

    def validate(self):
        par_format = {'dc':'6.2f', 'f':'6.2f'}
        super().validate(par_format)
        # if len(self.parameters) != 2:
        #     raise ValueError(f'{self.__class__.__name__} needs exactly 2 '
        #                      f'paramaters, not {len(self.parameters)}')


class Hernquist(DarkComponent):

    def __init__(self, **kwds):
        self.legacy_dm_code = 2
        super().__init__(symmetry='spherical', **kwds)

    def validate(self):
        par_format = {'rhoc':'6.2f', 'rc':'6.2f'}
        super().validate(par_format)
        # if len(self.parameters) != 2:
        #     raise ValueError(f'{self.__class__.__name__} needs exactly 2 '
        #                      f'paramaters, not {len(self.parameters)}')


class TriaxialCoredLogPotential(DarkComponent):

    def __init__(self, **kwds):
        self.legacy_dm_code = 3
        super().__init__(symmetry='triaxial', **kwds)

    def validate(self):
        par_format = {'Vc':'6.2f', 'rho':'6.2f', 'p':'6.2f', 'q':'6.2f'}
        super().validate(par_format)
        # if len(self.parameters) != 4:
        #     raise ValueError(f'{self.__class__.__name__} needs exactly 4 '
        #                      f'paramaters, not {len(self.parameters)}')


class GeneralisedNFW(DarkComponent):

    def __init__(self, **kwds):
        self.legacy_dm_code = 5
        super().__init__(symmetry='triaxial', **kwds)

    def validate(self):
        par_format = {'concentration':'6.2f', 'Mvir':'6.2f',
                      'inner_log_slope':'6.2f'}
        super().validate(par_format)
        # if len(self.parameters) != 3:
        #     raise ValueError(f'{self.__class__.__name__} needs exactly 3 '
        #                      f'paramaters, not {len(self.parameters)}')





# end
