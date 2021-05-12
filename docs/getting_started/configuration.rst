.. _configuration:

******************
Configuration
******************

All settings for running DYNAMITE can be controlled from a single configuration file in `yaml <https://en.wikipedia.org/wiki/YAML>`_ format. The basic structure of a yaml file are pairs of keys and values::

  key : value

which can be organised into hierarchical levels separated by tabs::

  main_key:
      sub_key1 : value1
      sub_key2 : value2

Comments are allowed, and begin with a ``#``. Values can be any type of variable e.g. integers, floats, strings, booleans etc. DYNAMITE configuration files must have the following hierarchy::

  # Example DYNAMITE configuration file e.g. called `config_file.yaml`
  ---
  ###############################################################################
  # SECTION 1 : SYSTEM
  # Define the physical system i.e. the galaxy being modelled
  ###############################################################################

  system_attributes:
      distMPc: ...        # distance in MPc
      name:  ...          # name for your galaxy
      position_angle:     # in degrees

  system_components:      # components of the system e.g. black hole, dark halo
      component_1: ...
      component_2: ...

  system_parameters:      # extra parameters of the system, unrelated to components
      ...

  ###############################################################################
  # SECTION 2: SETTINGS
  # Define other settings e.g. for the orbit library and weight solver
  ###############################################################################

  orblib_settings:
      ...

  weight_solver_settings:
      ...

  parameter_space_settings:
      ...

  io_settings:              # paths can be given with or without trailing slash
      input_directory: "input_files/"
      output_directory: "output/"
      all_models_file: "all_models.ecsv"

  multiprocessing_settings:
      ncpus: 4  # int or 'all_available'

  legacy_settings:          # where to find the `legacy` fortran programs
      directory: "default"  # "default" or an absolute path to alternative location

You can read the configuration file into a configuration object ``c`` as follows

.. code-block:: python

  import dynamite as dyn
  c = dyn.config_reader.Configuration('config_file.yaml') # read the configuration fie

The following sections go through each section of the configuration file and enumerate all the options available for that section,

1. `system_attributes`_
2. `system_components`_
3. `system_parameters`_
4. `orblib_settings`_
5. `weight_solver_settings`_
6. `parameter_space_settings`_
7. `io_settings`_
8. `multiprocessing_settings`_
9. `legacy_settings`_

Examples of completed configuration files for different scenarios can be found in the `tutorials <https://github.com/dynamics-of-stellar-systems/dynamite_release/tree/master/docs/tutorial_notebooks>`_. You may like to use these as templates for your own models.

``system_attributes``
=====================

This section requires three entries::

  system_attributes:
      distMPc: ...        # distance in MPc
      name:  ...          # name for your galaxy
      position_angle:     # in degrees

which can be accessed in the configuration object as ``c.system.distMPc`` etc.

``system_components``
=====================

The system is specified as a number of components. Each component needs entries for

- ``component name``: descriptive, but preferably short as this will be used in output directory names
    - ``type``: a string corresponding to one of the types in `component types`_
    - ``contributes_to_potential``: Boolean
    - ``include``: Boolean, whether to include this component or not. If False, equivalent to omitting this component entirely
    - ``parameters``, each requiring entries for
        - ``fixed``: Boolean, whether the parameter is to be kept fixed
        - ``value``: an initial value for the parameter
        - ``par_generator_settings``: settings controlling parameter search (can be omitted if ``fixed=True``)
            - ``lo``: minimum value
            - ``hi``: maximum value
            - ``step``: initial step size for parameter search
            - ``minstep``: minimum allowed stepsize for this parameter
        - ``logarithmic``: Boolean, whether logarithmic steps should be used for parameter search. If true, then (``value``, ``lo``, ``hi``) must all have log units
        - ``LaTeX``: LaTeX string for this parameter to be used for plots

``component types``
^^^^^^^^^^^^^^^^^^^^

The following component types are available, listed with their parameters:

- ``TriaxialVisibleComponent``
    - ``p``: intrinsic axis ratio B/A
    - ``q``: intrinsic axis ratio C/A
    - ``u``: sigma_observed / sigma_intrinsic
    - additionally, you must specify `observed data`_ for this component
- ``Plummer``
    - ``a``: scale length [CHECK UNITS]
    - ``m``: mass [solar masses]
- ``NFW``
    - ``c``: concentration parameter
    - ``f``: dark matter fraction ``M_200/M_{star}``
- ``Hernquist``
    - ``rhoc``: central density [CHECK UNITS]
    - ``rc``: scale length [CHECK UNITS]
- ``TriaxialCoredLogPotential``
    - ``p``: intrinsic axis ratio B/A
    - ``q``: intrinsic axis ratio C/A
    - ``rho``: ...? [CHECK UNITS]
    - ``Vc``: circular velocity at ...? [km/s CHECK UNITS]
- ``GeneralisedNFW``
    - ``concentration``: concentration parameter
    - ``Mvir``: virial mass [solar masses]
    - ``inner_log_slope``: central density slope [log [Msol/kpc^3] / log [kpc] CHECK...?]

**Note**: currently (v1.0) there is a strict requirement on the types of components that must be present in order to be compatible with the legacy - i.e. Fortran - implementation of the orbit integrator. The system must contain exactly:

- 1 ``Plummer`` component
    - representing the black hole (therefore with scale length ``a`` fixed to some arbitarily small value)
- 1 ``TriaxialVisibleComponent`` component
    - representing the stars
- 1 component representing the dark halo
    - either ``NFW``, ``Hernquist``, ``TriaxialCoredLogPotential``, or ``GeneralisedNFW``

``observed data``
^^^^^^^^^^^^^^^^^^^^

The ``TriaxialVisibleComponent`` represents the galaxy's stars, and therefore has associated observations. You must specify the following entries with filenames for observed data:

- ``TriaxialVisibleComponent``
    - ``mge_lum``: string, filename for the MGE of the projected luminosity density
    - ``mge_pot``: string, filename for the MGE of the projected mass density. This can be the same as ``mge_lum``, in which case there is an assumption that stellar-mass follows stellar-light
    - ``kinematics``
        - ``name of the kinematic set``
            - ``type``: type of kinematics - either ``GaussHermite`` or ``BayesLOSVD``
            - ``weight``: float, weighting applied to this kinematic set in chi2 calculation
            - ``datafile``: string, filename for the kinematics ECSV data file
            - ``aperturefile``: string, filename of the aperture file for this kinematic set
            - ``binfile``: string, filename of the bin file for this kinematic set
            - ``hist_width``: *optional*, float or 'default', the width (i.e. min. to max. value) of the velocity histogram for storing orbits. The default option is slightly wider than range of observed kinematics.
            - ``hist_center``: *optional*, float or 'default', the center of the velocity histogram for storing orbits. The default option = 0. This should not be changed.
            - ``hist_bins``: *optional*, int or 'default', the number of bins in the velocity histogram for storing orbits. The default option gives about 10 times better velocity sampling than the data.


``accessing the components``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The components are stored in the configuration object ``c`` as a list ``c.system.cmp_list``, with one entry for each of the components.

``system_parameters``
=====================

These are parameters of the system unrelated to any particular component. Currently there is only one such parameter, ``ml``. This acronym stands for mass-to-light, and indeed this parameter is the scaling between the surface luminosity of the galaxy (in $L_\odot \mathrm{pc}^{-3}$ CHECK UNITS) and the surface density (in $M_\odot \mathrm{pc}^{-3}$ CHECK UNITS). The ``ml`` parameter, however scales **all** mass components of the system, not just the stars. That is, say the system has a ``GeneralisedNFW`` component with ``Mvir=100`` but the system's ``ml`` parameter is equal to 2. The ``GeneralisedNFW`` would therefore physically represent a halo with mass ``Mvir=200``. This trick is employed in order to be able to re-use orbit-libraries to investigate a range of different potentials.

Specifying ``ml`` in the configuration file follows the same pattern as other parameters,

- ``system_parameters``
    - ``ml``
        - ``fixed``: Boolean, whether ``ml`` is to be kept fixed
        - ``value``: an initial value for ``ml``
        - ``par_generator_settings``: settings controlling parameter search (can be omitted if ``fixed=True``)
            - ``lo``: minimum value
            - ``hi``: maximum value
            - ``step``: initial step size for parameter search
            - ``minstep``: minimum allowed stepsize for this parameter
        - ``logarithmic``: Boolean, whether logarithmic steps should be used for parameter search. If true, then (``value``, ``lo``, ``hi``) must all have log units
        - ``LaTeX``: LaTeX format string for this parameter to be used for plots, e.g. in axis labels


``orblib_settings``
=====================

Hello!

``weight_solver_settings``
==========================

Hello!

``parameter_space_settings``
============================

Hello!

``io_settings``
=====================

Hello!

``multiprocessing_settings``
============================

Hello!


``legacy_settings``
=====================

Hello!
