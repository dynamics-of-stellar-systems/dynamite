# classes to hold the physical components of the system
# e.g. the stellar light, dark matter, black hole, globular clusters

import numpy as np
from scipy import special

# some tricks to add the current path to sys.path (so the imports below work)

import os.path
import sys
import logging

this_dir = os.path.dirname(__file__)
if not this_dir in sys.path:
    sys.path.append(this_dir)

import mges as mge


class System(object):

    def __init__(self, *args):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
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
        no duplicate component names, and the ml parameter.
        Additionally, the sformat string for the ml parameter is set.

        Raises
        ------
        ValueError : if required attributes or components are missing, or if
                     there is no ml parameter

        Returns
        -------
        None.

        """
        if len(self.cmp_list) != len(set(self.cmp_list)):
            raise ValueError('No duplicate component names allowed')
        # if self.parameters is not None: # Restriction should not be needed...
        #     for p in self.parameters:
        #         if any([p.name.endswith(c.name) for c in self.cmp_list]):
        #             raise ValueError('System parameter cannot end with '
        #                              f'"component": {p.name}')
        if not(self.distMPc and self.name and self.position_angle):
            text = 'System needs distMPc, name, and position_angle attributes'
            self.logger.error(text)
            raise ValueError(text)
        if not self.cmp_list:
            text = 'System has no components'
            self.logger.error(text)
            raise ValueError(text)
        if any(['_' in c.name for c in self.cmp_list]):
            self.logger.warning('System components should not contain '
                'underscores - model directory names may get confusing')
        if len(self.parameters) != 1 and self.parameters[0].name != 'ml':
            text = 'System needs ml as its sole parameter'
            self.logger.error(text)
            raise ValueError(text)
        self.parameters[0].update(sformat = '6.2f')

    def validate_parset(self, par):
        """
        Validates the system's parameter values. Kept separate from the
        validate method to facilitate easy calling from the parameter
        generator class.

        Parameters
        ----------
        par : dict
            { "p":val, ... } where "p" are the system's parameters and
            val are their respective values

        Returns
        -------
        isvalid : bool
            True if the parameter set is valid, False otherwise

        """
        isvalid = np.all(np.sign(tuple(par.values())) >= 0)
        return bool(isvalid)

    def __repr__(self):
        return f'{self.__class__.__name__} with {self.__dict__}'

    def get_component_from_name(self, cmp_name):
        cmp_list_list = np.array([cmp0.name for cmp0 in self.cmp_list])
        idx = np.where(cmp_list_list == cmp_name)
        self.logger.debug(f'Checking for 1 and only 1 component {cmp_name}...')
        error_msg = f"There should be 1 and only 1 component named {cmp_name}"
        assert len(idx[0]) == 1, error_msg
        self.logger.debug('...check ok.')
        component = self.cmp_list[idx[0][0]]
        return component

    def get_component_from_class(self, cmp_class):
        self.logger.debug('Checking for 1 and only 1 component of class '
                          f'{cmp_class}...')
        components = filter(lambda c: isinstance(c,cmp_class), self.cmp_list)
        component = next(components, False)
        if component is False or next(components, False) is not False:
            error_msg = 'Actually... there should be 1 and only 1 ' \
                        f'component of class {cmp_class}'
            self.logger.error(error_msg)
            raise ValueError(error_msg)
        self.logger.debug('...check ok.')
        return component

    def get_all_kinematic_data(self):
        all_kinematics = []
        for component in self.cmp_list:
            all_kinematics += component.kinematic_data
        return all_kinematics

class Component(object):

    def __init__(self,
                 name = None,                    # string
                 visible=None,                   # Boolean
                 contributes_to_potential=None,  # Boolean
                 symmetry=None,       # OPTIONAL, spherical, axisymm, or triax
                 kinematic_data=[],              # a list of Kinematic objects
                 population_data=[],             # a list of Population objects
                 parameters=[]):                 # a list of Parameter objects
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
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
            text = errstr + 'visible'
            self.logger.error(text)
            raise ValueError(text)
        if self.contributes_to_potential is None:
            text = errstr + 'contributes_to_potential'
            self.logger.error(text)
            raise ValueError(text)
        if not self.parameters:
            text = errstr + 'parameters'
            self.logger.error(text)
            raise ValueError(text)

        pars = [self.get_parname(p.name) for p in self.parameters]
        if set(pars) != set(par):
            text = f'{self.__class__.__name__} needs parameters ' + \
                   f'{list(par)}, not ' + \
                   f'{[self.get_parname(p.name) for p in self.parameters]}'
            self.logger.error(text)
            raise ValueError(text)

        self.set_format(par_format)

    def set_format(self, par_format=None):
        if par_format is None:
            text = f'{self.__class__.__name__}: no format string'
            self.logger.error(text)
            raise ValueError(text)
        for p in self.parameters:
            p.update(sformat=par_format[self.get_parname(p.name)])

    def validate_parset(self, par):
        """
        Validates the component's parameter values. Kept separate from the
        validate method to facilitate easy calling from the parameter
        generator class. This is a `placeholder` method which returns
        `True` if all parameters are non-negative. Specific validation
        should be implemented for each component subclass.

        Parameters
        ----------
        par : dict
            { "p":val, ... } where "p" are the component's parameters and
            val are their respective values

        Returns
        -------
        isvalid : bool
            True if the parameter set is valid, False otherwise

        """
        isvalid = np.all(np.sign(tuple(par.values())) >= 0)
        if not isvalid:
            self.logger.debug(f'Not a non-negative parset: {par}')
        return isvalid

    def get_parname(self, par):
        """
        Strips the component name suffix from the parameter name.

        Parameters
        ----------
        par : str
            The full parameter name "parameter-component".

        Returns
        -------
        pure_parname : str
            The parameter name without the component name suffix.

        """
        try:
            pure_parname = par[:par.rindex(f'-{self.name}')]
        except:
            self.logger.error(f'Component name {self.name} not found in '
                              f'parameter string {par}')
            raise
        return pure_parname

    def __repr__(self):
        return (f'\n{self.__class__.__name__}({self.__dict__}\n)')


class VisibleComponent(Component):

    def __init__(self,
                 mge_pot=None,
                 mge_lum=None,
                 **kwds):
         # visible components have MGE surface density
        self.mge_pot = mge_pot
        self.mge_lum = mge_lum
        super().__init__(visible=True, **kwds)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')

    def validate(self, **kwds):
        super().validate(**kwds)
        if not (isinstance(self.mge_pot, mge.MGE) and \
                isinstance(self.mge_lum, mge.MGE)):
            text = f'{self.__class__.__name__}.mge_pot and ' \
                   f'{self.__class__.__name__}.mge_lum ' \
                    'must be mges.MGE objects'
            self.logger.error(text)
            raise ValueError(text)
        if len(self.mge_pot.data) != len(self.mge_lum.data):
            text = f'{self.__class__.__name__}.mge_pot and ' \
                   f'{self.__class__.__name__}.mge_lum ' \
                    'must be of equal length'
            self.logger.error(text)
            raise ValueError(text)


class AxisymmetricVisibleComponent(VisibleComponent):

    def __init__(self, **kwds):
        super().__init__(symmetry='axisymm', **kwds)

    def validate(self):
        par_format = {'par1':'6.3g', 'par2':'6.4g'}
        super().validate(par_format)


class TriaxialVisibleComponent(VisibleComponent):

    def __init__(self, **kwds):
        super().__init__(symmetry='triax', **kwds)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.qobs = np.nan

    def validate(self):
        """
        In addition to validating parameter names and setting their sformat
        strings, also set self.qobs (minimal flattening from mge data)

        Returns
        -------
        None.

        """
        par_format = {'q':'6.3g', 'p':'6.4g', 'u':'7.5g'}
        super().validate(par_format=par_format)
        self.qobs = np.amin(self.mge_pot.data['q'])
        if self.qobs is np.nan:
            raise ValueError(f'{self.__class__.__name__}.qobs is np.nan')

    def validate_parset(self, par):
        """
        Validates the triaxial component's p, q, u parameter set. Requires
        self.qobs to be set. A parameter set is valid if the resulting
        (theta, psi, phi) are not np.nan.

        Parameters
        ----------
        par : dict
            { "p":val, ... } where "p" are the component's parameters and
            val are their respective values

        Returns
        -------
        bool
            True if the parameter set is valid, False otherwise

        """
        tpp = self.triax_pqu2tpp(par['p'], par['q'], par['u'])
        return bool(not np.any(np.isnan(tpp)))

    def triax_pqu2tpp(self,p,q,u):
        """
        transfer (p, q, u) to the three viewing angles (theta, psi, phi)
        with known flatting self.qobs.
        Taken from schw_basics, same as in vdB et al. 2008, MNRAS 385,2,647
        We should possibly revisit the expressions later

        """

        # avoid legacy_fortran's u=1 (rather, phi=psi=90deg) problem
        if u == 1:
            u *= (1-np.finfo(float).epsneg)  # same value as for np.double

        p2 = np.double(p) ** 2
        q2 = np.double(q) ** 2
        u2 = np.double(u) ** 2
        o2 = np.double(self.qobs) ** 2

        # Check for possible triaxial deprojection (v. d. Bosch 2004,
        # triaxpotent.f90 and v. d. Bosch et al. 2008, MNRAS 385, 2, 647)
        str = f'{q} <= {p} <= {1}, ' \
              f'{max((q/self.qobs,p))} <= {u} <= {min((p/self.qobs),1)}, ' \
              f'q\'={self.qobs}'
        # 0<=t<=1, t = (1-p2)/(1-q2) and p,q>0 is the same as 0<q<=p<=1 and q<1
        t = (1-p2)/(1-q2)
        if not (0 <= t <= 1) or \
           not (max((q/self.qobs,p)) <= u <= min((p/self.qobs),1)) :
            theta = phi = psi = np.nan
            self.logger.debug(f'DEPROJ FAIL: {str}')
        else:
            self.logger.debug(f'DEPROJ PASS: {str}')
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

        self.logger.debug(f'theta={theta}, phi={phi}, psi={psi}')
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

    def mass_enclosed(R, pars):
        M, a = pars
        Menc = M*R**3/a**3*(1 + R**2/a**2)**(-1.5)
        return Menc

    def validate(self):
        par_format = {'m':'6.3g', 'a':'7.3g'}
        super().validate(par_format)


class NFW(DarkComponent):

    def __init__(self, **kwds):
        self.legacy_code = 1
        super().__init__(symmetry='spherical', **kwds)

    def validate(self):
        par_format = {'c':'6.3g', 'f':'6.3g'}
        super().validate(par_format)


    # # c is concentration, f is dark mass fraction
    # def rhoc(c,f):
        # return 200/3 * rhocrit * c**3 / (log(1 + c) - c/(1+c))
    # def rc(c,f):
        # return (3*M200(c,f)/(800*pi*rhocrit*c**3))**(1/3)
    # def M200(c,f):
        # return 800*pi/3*rhocrit*(rc*c)**3
# 
    # def potential(x, y, z, pars):
        # c, f = pars
        # d2 = x**2 + y**2 + z**2
        # prefactor = 4*pi*G*rhoc(c,f)*(rc(c,f)**3)/sqrt(d2)
        # if sqrt(d2)/rc >= 1:
            # return prefactor * log(1 + sqrt(d2)/rc)
        # else:
            # return prefactor * 2 * atanh(sqrt(d2)/(2*rc(c,f) + sqrt(d2)))
# 
    # def density(x, y, z, pars):
        # c, f = pars
        # r = np.sqrt(x**2 + y**2 + z**2)
        # rho = rc(c,f)**3*rhoc(c,f)/(r*(r+rc(c,f))**2)
        # return rho
# 
    # def mass_enclosed(x, y, z, pars):
        # c, f = pars
        # r = np.sqrt(x**2 + y**2 + z**2)
        # Menc = 4*np.pi*rc(c,f)**3*rhoc(c,f)*(np.log(1 + r/rc(c,f)) - (r/rc(c,f))/(1 + r/rc(c,f)))
        # return Menc

class Hernquist(DarkComponent):

    def __init__(self, **kwds):
        self.legacy_dm_code = 2
        super().__init__(symmetry='spherical', **kwds)

    def validate(self):
        par_format = {'rhoc':'6.3g', 'rc':'6.3g'}
        super().validate(par_format)

    def potential(x, y, z, pars):
        rhoc, rc = pars
        r = np.sqrt(x**2 + y**2 + z**2)
        psi = 2*np.pi*G*rhoc*rc**2/(1 + r/rc)
        return psi

    def density(x, y, z, pars):
        rhoc, rc = pars
        r = np.sqrt(x**2 + y**2 + z**2)
        rho = rc**4*rhoc/(r*(r+rc)**3)
        return rho

    def mass_enclosed(x, y, z, pars):
        rhoc, rc = pars
        Menc = 2*np.pi*r**2*rc**3*rhoc/(r + rc)**2
        return Menc

class TriaxialCoredLogPotential(DarkComponent):

    def __init__(self, **kwds):
        self.legacy_dm_code = 3
        super().__init__(symmetry='triaxial', **kwds)

    def validate(self):
        par_format = {'Vc':'6.3g', 'rho':'6.3g', 'p':'6.3g', 'q':'6.3g'}
        super().validate(par_format)

    def potential(x, y, z, pars):
        rc, vc, p, q = pars
        m = x**2 + y**2/p**2 + z**2/q**2
        psi = -0.5*vc**2*np.log(rc**2 + m)
        return psi

    def density(x, y, z, pars):
        rc, vc, p, q = pars
        m = x**2 + y**2/p**2 + z**2/q**2
        rho = vc**2/(4*np.pi*G*(m+rc**2)**2)*( (m+rc**2)*(1 + 1/p**2 + 1/q**2) - 2*(x**2 + y**2/p**4 + z**2/q**4))
        return rho

    # this implementation assumes 1 > p > q
    def mass_enclosed(x, y, z, pars):
        rc, vc, p, q = pars
        r = np.sqrt(x**2 + y**2 + z**2)
        xx = r**2/(r**2/q**2 + rc**2)
        yy = r**2/(r**2/p**2 + rc**2)
        zz = r**2/(r**2 + rc**2)
        phi = np.arccos(np.sqrt(xx/zz))
        m = (zz-yy)/(zz-xx)
        Menc = r*vc**2/G * (1 - rc**2*r/np.sqrt((r**2+rc**2)*(r**2/p**2+rc**2)*(r**2/q**2+rc**2))*(zz-xx)**(-0.5)*special.ellipkinc(phi,m))
        

class GeneralisedNFW(DarkComponent):

    def __init__(self, **kwds):
        self.legacy_dm_code = 5
        super().__init__(symmetry='triaxial', **kwds)

    def validate(self):
        par_format = {'concentration':'6.3g', 'Mvir':'6.3g',
                      'inner_log_slope':'6.3g'}
        super().validate(par_format)

    def potential(x, y, z, pars):
        rhoc, rc, gamma = pars
        xi = r/(r + rc)
        psi = 4*np.pi*G*rhoc*rc**2*(rc/r*xi**(3-gamma)/(3-gamma)*special.hyp2f1(3-gamma,1,4-gamma,xi) + (1 - xi**(3-gamma))/(2-gamma))
        return psi

    def density(x, y, z, pars):
        rhoc, rc, gamma = pars
        rho = rhoc*rc**3*r**(-gamma)*(r + rc)**(gamma-3)
        return rho

    def mass_enclosed(x, y, z, pars):
        rhoc, rc, gamma = pars
        xi = r/(r + rc)
        Menc = 4*np.pi*rc**3*rhoc*xi**(3-gamma)/(3-gamma)*special.hyp2f1(3-gamma,1,4-gamma,xi)
        


# end
