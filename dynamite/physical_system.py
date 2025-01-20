# classes to hold the physical components of the system
# e.g. the stellar light, dark matter, black hole, globular clusters

import numpy as np
from scipy import special

# some tricks to add the current path to sys.path (so the imports below work)

import os.path
import sys
import logging

from dynamite import mges as mge
from dynamite import constants as const
from dynamite import orblib as orb

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
        sformat string for the ml parameter is set. Also validates
        the parameter generator settings' minstep value for ml if it is a
        non-fixed parameter.

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
        if (self.distMPc is None) or (self.name is None):
            text = 'System needs distMPc and name attributes'
            self.logger.error(text)
            raise ValueError(text)
        if not self.cmp_list:
            text = 'System has no components'
            self.logger.error(text)
            raise ValueError(text)
        #if len(self.parameters) != 1 and self.parameters[0].name != 'ml':
        if self.parameters[0].name != 'ml':
            text = 'System needs ml as its first parameter'
            self.logger.error(text)
            raise ValueError(text)
        ml = self.parameters[0]
        ml.update(sformat = '05.2f') # sformat of ml parameter
        if not ml.fixed and 'minstep' in ml.par_generator_settings:
            generator_settings = ml.par_generator_settings
            if generator_settings['minstep'] > generator_settings['step']:
                text = f"{self.__class__.__name__} parameter {ml.name}'s " \
                       "parameter generator settings have minstep > step, " \
                       f"setting minstep=step={generator_settings['step']}."
                self.logger.warning(text)
                generator_settings['minstep'] = generator_settings['step']
        if len(self.parameters) > 1:
            omega = self.parameters[1]
            if self.parameters[1].name != 'omega':
                text = 'System needs omega as its second parameter'
                self.logger.error(text)
                raise ValueError(text)

    def validate_parset(self, par):
        """
        Validates the system's parameter values

        Kept separate from the validate method to facilitate easy calling from
        the ``ParameterGenerator`` class. Returns `True` if all parameters are
        non-negative, except for logarithmic parameters which are not checked.

        Parameters
        ----------
        par : dict
            { "p":val, ... } where "p" are the system's parameters and
            val are their respective raw values

        Returns
        -------
        isvalid : bool
            True if the parameter set is valid, False otherwise

        """
        p_raw_values = [par[p.name]
                        for p in self.parameters if not p.logarithmic]
        isvalid = np.all(np.sign(p_raw_values) >= 0)
        if not isvalid:
            self.logger.debug(f'Invalid system parameters {par}: at least '
                              'one negative non-log parameter.')
        return bool(isvalid)

    def get_par_by_name(self, n):
        """
        Get a parameter using its name.

        Parameters
        ----------
        n : str
            The parameter name.

        Returns
        -------
        p : a ``dyn.parameter_space.Parameter`` object
            The parameter in question.

        """
        ps = self.parameters
        return ps[[p.name for p in ps].index(n)]

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
        ------
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

    def get_all_mge_components(self):
        """Get all components which contain MGEs

        Returns
        -------
        list
            a list of Component objects

        """
        mge_cmp = [c for c in self.cmp_list if isinstance(c, TriaxialVisibleComponent) or isinstance(c, BarDiskComponent)]
        return mge_cmp

    def get_unique_triaxial_visible_component(self):
        """Return the unique non-bar MGE component (raises an error if there are zero or multiple such components)

        Returns
        -------
            a ``dyn.physical_system.TriaxialVisibleComponent`` object

        """
        mges = self.get_all_mge_components()
        if len(mges) > 1:
            self.logger.error('Found more than one triaxial visible component')
            raise ValueError('Found more than one triaxial visible component')
        if len(mges) == 0:
            self.logger.error('Found zero triaxial visible components')
            raise ValueError('Found zero triaxial visible components')
        return mges[0]

    def get_all_bar_components(self):
        """Get all components which are rotating MGEs (i.e. bars)

        Returns
        -------
        list
            a list of Component objects, keeping only the rotating MGE components

        """
        bar_cmp = [c for c in self.cmp_list if isinstance(c, BarDiskComponent)]
        return bar_cmp

    def get_unique_bar_component(self):
        """Return the unique rotating bar component (raises an error if there are zero or multiple such components)

        Returns
        -------
            a ``dyn.physical_system.BarDiskComponent`` object

        """
        bars = self.get_all_bar_components()
        if len(bars) > 1:
            self.logger.error('Found more than one bar component')
            raise ValueError('Found more than one bar component')
        if len(bars) == 0:
            self.logger.error('Found zero bar components')
            raise ValueError('Found zero bar components')
        return bars[0]

    def get_all_dark_components(self):
        """Get all components which are Dark

        Returns
        -------
        list
            a list of Component objects, keeping only the dark components

        """
        dark_cmp = [c for c in self.cmp_list if isinstance(c, DarkComponent)]
        return dark_cmp

    def get_all_dark_non_plummer_components(self):
        """Get all Dark components which are not plummer

        Useful in legacy orbit libraries for finding the dark halo component.
        For legacy models, the black hole is always a plummer, so any Dark but
        non plummer components must represent the dark halo.

        Returns
        -------
        list
            a list of Component objects, keeping only the dark components

        """
        dark_cmp = self.get_all_dark_components()
        dark_non_plum_cmp = [c for c in dark_cmp if not isinstance(c, Plummer)]
        return dark_non_plum_cmp

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

    def is_bar_disk_system(self):
        """is_bar_disk_system

        Check if the system contains at least one bar component and at least one disk component.

        Returns
        -------
        isbardisk : Bool
            System contains a bar and a disk.
        """
        isbardisk = len(self.get_all_bar_components()) > 0
        return isbardisk

    def is_bar_disk_system_with_angles(self):
        """is_bar_disk_system_with_angles

        Check if the system is a bar-disk with phi, psi, theta specified directly in the configuration file.

        Returns
        -------
        hasangles : Bool
            System is specified by angles.
        """
        hasangles = self.is_bar_disk_system() and (type(self.get_unique_bar_component()) is BarDiskComponentAngles)
        return hasangles

    def number_of_visible_components(self):
        return sum(1 for i in self.cmp_list if isinstance(i, VisibleComponent))

    def number_of_bar_components(self):
        return sum(1 for i in self.cmp_list if isinstance(i, BarDiskComponent))

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
        self.logger.info(f'{self.name}: DYNAMITE will currently ignore the '
                         'mandatory attribute contributes_to_potential.')
        self.symmetry = symmetry
        self.kinematic_data = kinematic_data
        self.population_data = population_data
        self.parameters = parameters

    def validate(self, par=None):
        """
        Validate the component

        Ensure it has the required attributes and parameters. Also validates
        the parameter generator settings' minstep value for non-fixed
        parameters.

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
            text = f'{self.__class__.__name__} needs parameter(s) ' + \
                   f'{par}, not {pars}.'
            self.logger.error(text)
            raise ValueError(text)

        for p in [p for p in self.parameters
                  if not p.fixed and 'minstep' in p.par_generator_settings]:
            generator_settings = p.par_generator_settings
            if generator_settings['minstep'] > generator_settings['step']:
                text = f"{self.__class__.__name__} parameter {p.name}'s " \
                       "parameter generator settings have minstep > step, " \
                       f"setting minstep=step={generator_settings['step']}."
                self.logger.warning(text)
                generator_settings['minstep'] = generator_settings['step']

    def validate_parset(self, par):
        """
        Validates the component's parameter values.

        Kept separate from the
        validate method to facilitate easy calling from the parameter
        generator class. This is a `placeholder` method which returns
        `True` if all parameters are non-negative, except for logarithmic
        parameters which are not checked. Specific validation
        should be implemented for each component subclass.

        Parameters
        ----------
        par : dict
            { "p":val, ... } where "p" are the component's parameters and
            val are their respective raw values

        Returns
        -------
        isvalid : bool
            True if the parameter set is valid, False otherwise

        """
        p_raw_values = [par[self.get_parname(p.name)]
                    for p in self.parameters if not p.logarithmic]
        isvalid = np.all(np.sign(p_raw_values) >= 0)
        if not isvalid:
            self.logger.debug(f'Invalid parameters {par}: at least one '
                              'negative non-log parameter.')
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

    def get_par_by_name(self, n):
        """
        Get a parameter using its (unsuffixed) name.

        Parameters
        ----------
        n : str
            The parameter name (without the component name suffix)

        Returns
        -------
        p : a ``dyn.parameter_space.Parameter`` object
            The parameter in question.

        """
        ps = self.parameters
        return ps[[p.name for p in ps].index(n + '-' + self.name)]

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

    def get_M_stars_tot(self, distance, parset):
        """
        Calculates and returns the total stellar mass via the mge.

        Parameters
        ----------
        distance : float
            Distance of the system in MPc
        parset : astropy table row
            must contain mass-to-light ratio ml

        Returns
        -------
        float
            Total stellar mass

        """
        mgepar = self.mge_pot.data
        mgeI = mgepar['I']
        mgesigma = mgepar['sigma']
        mgeq = mgepar['q']

        arctpc = distance*np.pi/0.648
        sigobs_pc = mgesigma*arctpc

        return 2 * np.pi * np.sum(mgeI * mgeq * sigobs_pc ** 2) * parset['ml']

    def intrin_spher(self, distance, parset):
        totalmass = self.get_M_stars_tot(distance, parset)
        mge = self.mge_pot
        ngauss_mge = len(mge.data)

        quad_nr = 10 # size of (r, th, ph) grid, hardcoded in Fortran
        quad_nth = 6
        quad_nph = 6

        quad_grid = np.zeros([quad_nph, quad_nth, quad_nr])
        quad_lr = np.zeros(quad_nr + 1)
        quad_lth = np.zeros(quad_nth + 1)
        quad_lph = np.zeros(quad_nph + 1)

    def intrin_radii(self, distance, parset, orblib_settings):
        totalmass = self.get_M_stars_tot(distance, parset)
        mge = self.mge_pot
        ngauss_mge = len(mge.data)
        qobs = mge.data['q']
        rlogmin = orblib_settings['logrmin']
        rlogmax = orblib_settings['logrmax']
        orbit_dithering = orblib_settings['dithering']
        nener = orblib_settings['nE']
        pintr, qintr, sigintr_km, dens = self.triax_deproj()

    def validate(self, **kwds):
        super().validate(**kwds)
        if not (isinstance(self.mge_pot, mge.MGE) and \
                isinstance(self.mge_lum, mge.MGE)):
            text = f'{self.__class__.__name__}.mge_pot and ' \
                   f'{self.__class__.__name__}.mge_lum ' \
                    'must be mges.MGE objects'
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
        self.par = ['q', 'p', 'u']

    def validate(self):
        """
        Validate the TriaxialVisibleComponent

        Validates parameter names and sets self.qobs
        (minimal flattening from mge data).

        Returns
        -------
        None.

        """
        super().validate(par=self.par)
        self.qobs = np.amin(self.mge_pot.data['q'])
        if np.isnan(self.qobs):
            raise ValueError(f'{self.__class__.__name__}.qobs is np.nan')

    def validate_parset(self, par):
        """
        Validate the p, q, u parameters

        Validates the triaxial component's p, q, u parameters. Requires
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
        # Check for p=0
        if np.isclose(p, 0.):
            theta = phi = psi = np.nan
            self.logger.debug(f'DEPROJ FAIL: p=0')
        # Check for q=0
        if np.isclose(q, 0.):
            theta = phi = psi = np.nan
            self.logger.debug(f'DEPROJ FAIL: q=0')
        # Check for u=0
        if np.isclose(u, 0.):
            theta = phi = psi = np.nan
            self.logger.debug(f'DEPROJ FAIL: u=0')
        # Check for u>1
        if u>1:
            theta = phi = psi = np.nan
            self.logger.debug(f'DEPROJ FAIL: u>1')
        if np.isclose(u,p):
            u=p
        # Check for possible triaxial deprojection (v. d. Bosch 2004,
        # triaxpotent.f90 and v. d. Bosch et al. 2008, MNRAS 385, 2, 647)
        str = f'{q} <= {p} <= {1}, ' \
              f'{max((q/self.qobs,p))} < {u} <= {min((p/self.qobs),1)}, ' \
              f'q\'={self.qobs}'
        # 0<=t<=1, t = (1-p2)/(1-q2) and p,q>0 is the same as 0<q<=p<=1 and q<1
        t = (1-p2)/(1-q2)
        if not (0 <= t <= 1) or \
           not (max((q/self.qobs,p)) < u <= min((p/self.qobs),1)) :
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

    def find_grid_of_valid_pqu(self, n_grid=200):
        """Find valid values of the parameters (p,q,u)

        Creates a grid of all values of 0<(p,q,u)<=1, and finds those which have
        a valid deprojection subject to the fulfilment of all three criteria:
        1. 0 < q <= p <=1
        2. max(q/qobs, p) < u
        3. u< min(p/qobs, 1)
        where qobs is the smallest value of q for the MGE.

        Parameters
        ----------
        n_grid : int
            grid size used for p,q,u

        Returns
        -------
        (p,q,u), valid
            3d grids of p,q,u values, and boolen array `valid` which is `True`
            for valid values

        """
        # make grid of possible p,q,u values
        p = np.linspace(0, 1, n_grid)[1:]
        q = np.linspace(0, 1, n_grid)[1:]
        u = np.linspace(0, 1, n_grid)[1:]
        p, q, u = np.meshgrid(p, q, u, indexing='ij')
        # check three conditions for whether (p,q,u) give a valid deprojection
        invalid_a = q>p
        invalid_b = np.maximum(q/self.qobs, p) >= u
        invalid_c = u >= np.minimum(p/self.qobs, 1.)
        # combine the conditions
        invalid_ab = np.logical_or(invalid_a, invalid_b)
        invalid_abc = np.logical_or(invalid_ab, invalid_c)
        valid = np.logical_not(invalid_abc)
        return (p,q,u), valid

    def suggest_parameter_values(self, target_u=0.9):
        """Suggest valid values of the parameters (p,q,u)

        Find valid values using the mehtod `find_grid_of_valid_pqu`. Then for
        each of (p,q,u), we suggest values:
        - lo/hi : the min/max of all valid values
        - value : u=target_u, and p/q = mean of all valid p/q values where u is
        close to target value
        - step/minstep : a fifth/twentieth of the range of valid values

        Parameters
        ----------
        target_u : float
            Desired value of the parameter u

        Returns
        -------
        string
            text to print out suggesting quantities for (p,q,u)
        """
        (p, q, u), valid = self.find_grid_of_valid_pqu()
        text = "No deprojection possible for the specificed values of (p,q,u)."
        text += "Here are some suggestions: \n"
        # take avg of valid p's and q's where u is close to targer value
        target_u = 0.9
        idx = np.where(np.abs(u[valid]-target_u)<0.005)
        if idx[0].shape==(0,):
            text = f"Cannot suggest valid (p,q,u) for a target u={target_u}"
            self.logger.error(text)
            raise ValueError(text)
        suggest_p = np.mean(p[valid][idx])
        suggest_q = np.mean(q[valid][idx])
        suggested_values = [suggest_p, suggest_q, target_u]
        for (symbol, array, val) in zip(['p', 'q', 'u'],
                                        [p[valid], q[valid], u[valid]],
                                        suggested_values):
            lo, hi = np.min(array), np.max(array)
            step = (hi-lo)/5.
            minstep = (hi-lo)/20.
            text += f'\t{symbol}:\n'
            text += f'\t\t lo : {lo:.2f}\n'
            text += f'\t\t hi : {hi:.2f}\n'
            text += f'\t\t step : {step:.2f}\n'
            text += f'\t\t minstep : {minstep:.2f}\n'
            text += f'\t\t value : {val:.2f}\n'
        return text

class BarDiskComponent(TriaxialVisibleComponent):
    """Rotating triaxial component with a MGE projected density (i.e. a bar)

    Note: all bar components are constrained to have the same omega.

    """
    def __init__(self,
                 mge_pot=None,
                 mge_lum=None,
                 disk_pot=None,
                 disk_lum=None,
                 **kwds):
        super().__init__(**kwds)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.qobs = np.nan
        self.par = ['q', 'p', 'u', 'qdisk']

class BarDiskComponentAngles(BarDiskComponent):
    """Rotating triaxial component with a MGE projected density (i.e. a bar),
    with viewing angles specified.

    Note: all bar components are constrained to have the same omega.

    """
    def __init__(self,
                 mge_pot=None,
                 mge_lum=None,
                 disk_pot=None,
                 disk_lum=None,
                 **kwds):
        super().__init__(**kwds)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.qobs = np.nan
        self.par = ['theta', 'psi', 'phi']

    def validate_parset(self, par):
        # Skip validation as we already know the angles
        return True

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

    def get_dh_legacy_strings(self, parset):
        """
        Generates and returns two strings needed for the legacy Fortran files.

        This method only applies to dark halo components.

        Parameters
        ----------
        parset : astropy table row
            Holds the parameter set.

        Returns
        -------
        specs : str
            A string with the legacy code and the number of parameters, space
            separated.
        par_vals : str
            The parameter values in the sequence legacy Fortran expects them,
            space separated.

        """
        try:
            legacy_code = self.legacy_code
            specs = f'{legacy_code} {len(self.parameters)}'
            par_vals = ''
            for par in self.par_names:
                p = f'{par}-{self.name}'
                par_vals += f'{parset[p]} '
            par_vals = par_vals[:-1]
            self.logger.debug(f'DH {self.__class__.__name__} legacy strings: '
                              f'{specs} / {par_vals}.')
            return specs, par_vals
        except AttributeError: # Only dh has a legacy code, Plummer: do nothing
            pass


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

    def mass_enclosed(R, pars):
        M, a = pars
        Menc = M*R**3/a**3*(1 + R**2/a**2)**(-1.5)
        return Menc

    def validate(self):
        par = ['m', 'a']
        super().validate(par=par)


class NFW(DarkComponent):
    """An NFW halo

    Defined with parameters: c [concentration, R200/scale] and f
    [dm-fraction, M200/total-stellar-mass]

    """
    par_names = ['c', 'f'] # parameter names in legacy sequence

    def __init__(self, **kwds):
        self.legacy_code = 1
        super().__init__(symmetry='spherical', **kwds)

    def validate(self):
        super().validate(par=self.par_names)


class NFW_m200_c(DarkComponent):
    """An NFW halo with m200-c relation from Dutton & Maccio 14

    The relation: log10(c200) = 0.905 - 0.101 * log10(M200/(1e12/h)).
    Component defined with parameter f [dm-fraction, M200/total-stellar-mass]

    """
    par_names = ['f'] # parameter names in legacy sequence

    def __init__(self, **kwds):
        self.legacy_code = 1
        super().__init__(symmetry='spherical', **kwds)

    def validate(self):
        super().validate(par=self.par_names)

    def get_c200(self, system, parset):
        """
        Calculates and returns c200 (see Dutton & Maccio 2014).

        Parameters
        ----------
        system : a ``dyn.physical_system.System`` object
        parset : astropy table row
            Must contain dark matter fraction f-{self.name} and ml

        Returns
        -------
        float
            c200

        """
        stars = system.get_component_from_class(TriaxialVisibleComponent)
        M_stars_tot = stars.get_M_stars_tot(system.distMPc, parset)
        f = parset[f'f-{self.name}']
        h=0.671 #add paper
        #total mass of dark matter
        MvDM = f * M_stars_tot
        #dutton&maccio2014 (https://arxiv.org/pdf/1402.7073.pdf) Eq. (8)
        lc200 = 0.905 - 0.101*np.log10( MvDM/(1e12/h))

        return 10.**lc200

    def get_dh_legacy_strings(self, parset, system):
        """
        Generates and returns two strings needed for the legacy Fortran files.

        This method overrides the parent class' method because for legacy
        Fortran purposes, NFW_m200_c has two parameters. Note that NFW_m200_c
        needs an addiional parameter ``system``.

        Parameters
        ----------
        parset : astropy table row
            Holds the parameter set.

        Returns
        -------
        specs : str
            A string with the legacy code and the number of parameters, space
            separated.
        par_vals : str
            The parameter values in the sequence legacy Fortran expects them,
            space separated.

        """
        specs, par_vals = super().get_dh_legacy_strings(parset)
        c200 = self.get_c200(system, parset)
        specs = f'{self.legacy_code} 2'
        par_vals = f'{c200} {par_vals}'
        self.logger.debug(f'DH {self.__class__.__name__} legacy strings '
                          f'amended to {specs} / {par_vals}.')
        return specs, par_vals


    # c is concentration, f is dark mass fraction
    ## fixme: should derive rhocrit from (c,f) (?)
    rhocrit = 1
    def rhoc(c,f):
        return 200/3 * rhocrit * c**3 / (log(1 + c) - c/(1+c))
    def rc(c,f):
        return (3*M200(c,f)/(800*pi*rhocrit*c**3))**(1/3)
    def M200(c,f):
        return 800*pi/3*rhocrit*(rc*c)**3

    def potential(x, y, z, pars):
        c, f = pars
        d2 = x**2 + y**2 + z**2
        prefactor = 4*pi*G*rhoc(c,f)*(rc(c,f)**3)/sqrt(d2)
        if sqrt(d2)/rc >= 1:
            return prefactor * log(1 + sqrt(d2)/rc)
        else:
            return prefactor * 2 * atanh(sqrt(d2)/(2*rc(c,f) + sqrt(d2)))

    def density(x, y, z, pars):
        c, f = pars
        r = np.sqrt(x**2 + y**2 + z**2)
        rho = rc(c,f)**3*rhoc(c,f)/(r*(r+rc(c,f))**2)
        return rho

    def mass_enclosed(x, y, z, pars):
        c, f = pars
        r = np.sqrt(x**2 + y**2 + z**2)
        Menc = 4*np.pi*rc(c,f)**3*rhoc(c,f)*(np.log(1 + r/rc(c,f)) - (r/rc(c,f))/(1 + r/rc(c,f)))
        return Menc

class Hernquist(DarkComponent):
    """A Hernquist sphere

    Defined with parameters: rhoc [central density, Msun/km^3] and rc [scale
    length, km]

    """
    par_names = ['rhoc', 'rc'] # parameter names in legacy sequence

    def __init__(self, **kwds):
        self.legacy_code = 2
        super().__init__(symmetry='spherical', **kwds)

    def validate(self):
        super().validate(par=self.par_names)

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
    """A TriaxialCoredLogPotential

    see e.g. Binney & Tremaine second edition p.171
    Defined with parameters: p [B/A], q [C/A], Rc [core radius, kpc], Vc
    [asympt. circular velovity, km/s]

    """
    par_names = ['Vc', 'Rc', 'p', 'q'] # parameter names in legacy sequence

    def __init__(self, **kwds):
        self.legacy_code = 3
        super().__init__(symmetry='triaxial', **kwds)

    def validate(self):
        super().validate(par=self.par_names)

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
    """A GeneralisedNFW halo

    from Zhao (1996)
    Defined with parameters: concentration [R200/NFW scale length], Mvir [Msol],
    inner_log_slope []

    """
    par_names = ['c', 'Mvir', 'gam'] # parameter names in legacy sequence

    def __init__(self, **kwds):
        self.legacy_code = 5
        super().__init__(symmetry='triaxial', **kwds)

    def validate(self):
        super().validate(par=self.par_names)

    def validate_parset(self, par):
        """
        Validates the GeneralisedNFW's parameters.

        Requires c and Mvir >0, and gam leq 1

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
        if (par['c']<0.) or (par['Mvir']<0.) or (par['gam']>1):
            is_valid = False
        else:
            is_valid = True
        return is_valid

    ## fixme: should actually derive (rhoc,rc) from (c,Mvir)
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
