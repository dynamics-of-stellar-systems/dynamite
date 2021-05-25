# classes to hold the physical components of the system
# e.g. the stellar light, dark matter, black hole, globular clusters

import numpy as np
import logging

from dynamite import mges as mge

class System(object):
    """The physical system being modelled

    e.g. system is a galaxy. A system is composed of ``Components`` e.g. the
    galaxy is composed of stars, black hole, dark matter halo. This object is
    automatically created when the configuration file is read.
    """
    def __init__(self, *args):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.n_cmp = 0
        self.cmp_list = []
        self.n_pot = 0
        self.n_kin = 0
        self.n_pop = 0
        self.parameters = None
        self.distMPc = None
        self.name = None
        self.position_angle = None
        for component in args:
            self.add_component(component)

    def add_component(self, cmp):
        """add a component to the system

        Parameters
        ----------
        cmp : a ``dyn.physical_system.Component`` object

        Returns
        -------
        None
            updated the system componenent attributes

        """
        self.cmp_list += [cmp]
        self.n_cmp += 1
        self.n_pot += cmp.contributes_to_potential
        self.n_kin += len(cmp.kinematic_data)
        self.n_pop += len(cmp.population_data)

    def validate(self):
        """
        Validate the system

        Ensures the System has the required attributes: at least one component,
        no duplicate component names, and the ml parameter, and that the
        sformat string for the ml parameter is set.

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
        if (self.distMPc is None) or (self.name is None) \
           or (self.position_angle is None):
            text = 'System needs distMPc, name, and position_angle attributes'
            self.logger.error(text)
            raise ValueError(text)
        if not self.cmp_list:
            text = 'System has no components'
            self.logger.error(text)
            raise ValueError(text)
        if len(self.parameters) != 1 and self.parameters[0].name != 'ml':
            text = 'System needs ml as its sole parameter'
            self.logger.error(text)
            raise ValueError(text)
        self.parameters[0].update(sformat = '01.2f') # sformat of ml parameter

    def validate_parset(self, par):
        """
        Validates the system's parameter values

        Kept separate from the validate method to facilitate easy calling from
        the ``ParameterGenerator`` class.

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
        """get_component_from_name

        Parameters
        ----------
        cmp_name : string
            component name (as specified in the congi file)

        Returns
        -------
        a ``dyn.physical_system.Component`` object

        """
        cmp_list_list = np.array([cmp0.name for cmp0 in self.cmp_list])
        idx = np.where(cmp_list_list == cmp_name)
        self.logger.debug(f'Checking for 1 and only 1 component {cmp_name}...')
        error_msg = f"There should be 1 and only 1 component named {cmp_name}"
        assert len(idx[0]) == 1, error_msg
        self.logger.debug('...check ok.')
        component = self.cmp_list[idx[0][0]]
        return component

    def get_component_from_class(self, cmp_class):
        """get_component_from_class

        Parameters
        ----------
        cmp_class : string
            name of the component type/class

        Raises
        -------
        ValueError : if there are more than one component of the same class.
            # TODO: remove this limit, e.g. if we had two MGE-based components
            one for stars, one for gas

        Returns
        -------
        a ``dyn.physical_system.Component`` object

        """
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
        """get_all_kinematic_data

        Loop over all components, extract their kinemtics into a list.

        Returns
        -------
        list
            all_kinematics in a list

        """
        all_kinematics = []
        for component in self.cmp_list:
            all_kinematics += component.kinematic_data
        return all_kinematics

class Component(object):
    """A component of the physical system

    e.g. the stellar component, black hole, or dark halo of a galaxy

    Parameters
    ----------
    name : string
        a short but descriptive name of the component
    visible : Bool
        whether this is visible <--> whether it has an associated MGE
    contributes_to_potential : Bool
        whether this contributes_to_potential **not currently used**
    symmetry : string
        one of 'spherical', 'axisymm', or 'triax' **not currently used**
    kinematic_data : list
        a list of ``dyn.kinemtics.Kinematic`` data for this component
    parameters  : list
        a list of ``dyn.parameter_space.Parameter`` objects for this component
    population_data : list
        a list of ``dyn.populations.Population`` data for this component **not
        currently used**

    """
    def __init__(self,
                 name = None,
                 visible=None,
                 contributes_to_potential=None,
                 symmetry=None,
                 kinematic_data=[],
                 population_data=[],
                 parameters=[]):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if name == None:
            self.name = self.__class__.__name__
        else:
            self.name = name
        self.visible = visible
        self.contributes_to_potential = contributes_to_potential
        self.symmetry = symmetry
        self.kinematic_data = kinematic_data
        self.population_data = population_data
        self.parameters = parameters

    def validate(self, par=None):
        """
        Validate the component

        Ensure it has the required attributes and parameters.
        Additionally, the sformat strings for the parameters are set.

        Parameters
        ----------
        par : a list with parameter names. Mandatory.

        Raises
        ------
        ValueError : if a required attribute is missing or the required
                     parameters do not exist

        Returns
        -------
        None.

        """
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
                   f'{par}, not {pars}.'
            self.logger.error(text)
            raise ValueError(text)

    def validate_parset(self, par):
        """
        Validates the component's parameter values.

        Kept separate from the
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
        Strip the component name suffix from the parameter name.

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
    """Any visible component of the sytem, with an MGE

    Parameters
    ----------
    mge_pot : a ``dyn.mges.MGE`` object
        describing the (projected) surface-mass density
    mge_lum : a ``dyn.mges.MGE`` object
        describing the (projected) surface-luminosity density

    """
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
        par = ['par1', 'par2']
        super().validate(par=par)


class TriaxialVisibleComponent(VisibleComponent):
    """Triaxial component with a MGE projected density

    Has parameters (p,q,u) = (b/a, c/a, sigma_obs/sigma_intrinsic) used for
    deprojecting the MGE. A given (p,q,u) correspond to a fixed set of
    `viewing angles` for the triaxial ellipsoid.

    """
    def __init__(self, **kwds):
        super().__init__(symmetry='triax', **kwds)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.qobs = np.nan

    def validate(self):
        """
        Validate the TriaxialVisibleComponent

        In addition to validating parameter names and setting their sformat
        strings, also set self.qobs (minimal flattening from mge data)

        Returns
        -------
        None.

        """
        par = ['q', 'p', 'u']
        super().validate(par=par)
        self.qobs = np.amin(self.mge_pot.data['q'])
        if self.qobs is np.nan:
            raise ValueError(f'{self.__class__.__name__}.qobs is np.nan')

    def validate_parset(self, par):
        """
        Validate the p, q, u parset

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
        transform axis ratios to viewing angles

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
    """Any dark component of the sytem, with no observed MGE or kinemtics

    This is an abstract layer and none of the attributes/methods are currently
    used.

    """
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
    """A Plummer sphere

    Defined with parameters: M [mass, Msol] and a [scale length, arcsec]

    """
    def __init__(self, **kwds):
        super().__init__(symmetry='spherical', **kwds)

    def density(x, y, z, pars):
        M, a = pars
        r = (x**2 + y**2 + z**2)**0.5
        rho = 3*M/4/np.pi/a**3 * (1. + (r/a)**2)**-2.5
        return rho

    def validate(self):
        par = ['m', 'a']
        super().validate(par=par)


class NFW(DarkComponent):
    """An NFW halo

    Defined with parameters: c [concentration, R200/scale] and f
    [dm-fraction, M200/total-stellar-mass]

    """
    def __init__(self, **kwds):
        self.legacy_code = 1
        super().__init__(symmetry='spherical', **kwds)

    def validate(self):
        par = ['c', 'f']
        super().validate(par=par)


class Hernquist(DarkComponent):
    """A Hernquist sphere

    Defined with parameters: rhoc [central density, Msun/km^3] and rc [scale
    length, km]

    """
    def __init__(self, **kwds):
        self.legacy_dm_code = 2
        super().__init__(symmetry='spherical', **kwds)

    def validate(self):
        par = ['rhoc', 'rc']
        super().validate(par=par)


class TriaxialCoredLogPotential(DarkComponent):
    """A TriaxialCoredLogPotential

    see e.g. Binney & Tremaine second edition p.171
    Defined with parameters: p [B/A], q [C/A], Rc [core radius, kpc], Vc
    [asympt. circular velovity, km/s]

    """
    def __init__(self, **kwds):
        self.legacy_dm_code = 3
        super().__init__(symmetry='triaxial', **kwds)

    def validate(self):
        par = ['Vc', 'Rc', 'p', 'q']
        super().validate(par=par)


class GeneralisedNFW(DarkComponent):
    """A GeneralisedNFW halo

    from Zhao (1996)
    Defined with parameters: concentration [R200/NFW scale length], Mvir [Msol],
    inner_log_slope []

    """
    def __init__(self, **kwds):
        self.legacy_dm_code = 5
        super().__init__(symmetry='triaxial', **kwds)

    def validate(self):
        par = ['concentration', 'Mvir', 'inner_log_slope']
        super().validate(par=par)


# end
