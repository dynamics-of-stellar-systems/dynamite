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

Comments are allowed, and begin with a ``#``. Values can be any type of variable e.g. integers, floats, strings, booleans etc. DYNAMITE configuration files must have the following structure::

  # Example DYNAMITE configuration file e.g. called `config_file.yaml`
  ---
  ###############################################################################
  # SECTION 1 : SYSTEM
  # Define the physical system i.e. the galaxy being modelled
  ###############################################################################

  system_attributes:      # e.g. name, distance, ...
      ...

  system_components:      # components of the system e.g. black hole, dark halo
      component_1: ...
      component_2: ...
      ...

  system_parameters:      # extra parameters of the system
      ...

  ###############################################################################
  # SECTION 2: SETTINGS
  # Define other settings e.g. for the orbit library and weight solver
  ###############################################################################

  orblib_settings:            # settings for the orbit library calculation
      ...

  weight_solver_settings:     # settings for solving for orbital weights
      ...

  parameter_space_settings:   # settings for parameter search
      ...

  multiprocessing_settings:   # settings for multiprocessing
      ...

  io_settings:                # settings for input/output locations
      ...

  legacy_settings:            # location of Fortran programs
      ...

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

Examples of completed configuration files can be found in the tutorials, which you may like to use as templates for your own models.

Information about how the DYNAMITE configuration settings can be accessed, which may be useful for developers, can be found in the API documentation.

``system_attributes``
=====================

This section lists the following attributes of the system::

  system_attributes:
      distMPc: ...        # distance in MPc
      name:  ...          # name for your galaxy

``system_components``
=====================

The system consists of a number of physical components - e.g. the stars, black hole, dark halo. For each component the following values must be specified

- ``component name``: a descriptive name, but preferably short as this will be used to refer to the component in the code (e.g. ``bh`` for black hole)
    - ``type``: a string corresponding to one of the options in in `component types`_
    - ``contributes_to_potential``: Boolean (not currently used)
    - ``include``: Boolean, whether to include this component or not. If False, equivalent to omitting this component entirely
    - ``parameters``. The required parameters for each component are listed in `component types`_. Each  parameter must have values specified for
        - ``fixed``: Boolean, whether the parameter is to be kept fixed
        - ``value``: an initial value for the parameter
        - ``par_generator_settings``: settings controlling parameter search (can be omitted if ``fixed=True``). Note that if these settings are given, then ``value`` must be consistent with ``lo`` and ``hi``.
            - ``lo``: minimum value
            - ``hi``: maximum value
            - ``step``: initial step size for parameter search
            - ``minstep``: minimum allowed stepsize for this parameter
        - ``logarithmic``: Boolean, whether logarithmic steps should be used for parameter search. If true, then (``value``, ``lo``, ``hi``) must all have log units.
        - ``LaTeX``: LaTeX string for this parameter to be used for plots.

There is a special component type ``Chi2Ext`` that allows external modules to calculate a :math:`\chi^2` value for DYNAMITE models. These external :math:`\chi^2` values are calculated after DYNAMITE's internal :math:`\chi^2` calculations then appended to DYNAMITE's all_models, then can be used to inform the parameter search for later model iterations. The ``Chi2Ext`` component has different settings from the ones built in DYNAMITE:

- ``component name``: a descriptive name, but preferably short as this will be used to refer to the component in the code (e.g. ``gas``)
    - ``type``: ``Chi2Ext``
    - ``contributes_to_potential``: optional, will always be set to False (not currently used)
    - ``include``: Boolean, whether to include this component or not. If False, equivalent to omitting this component entirely
    - ``parameters``: Choose any parameter names as you see fit (e.g., par1, par2,...). Each parameter must have values specified for
        - ``fixed``: Boolean, whether the parameter is to be kept fixed
        - ``value``: an initial value for the parameter
        - ``par_generator_settings``: settings controlling parameter search (can be omitted if ``fixed=True``). Note that if these settings are given, then ``value`` must be consistent with ``lo`` and ``hi``.
            - ``lo``: minimum value
            - ``hi``: maximum value
            - ``step``: initial step size for parameter search
            - ``minstep``: minimum allowed stepsize for this parameter
        - ``logarithmic``: Boolean, whether logarithmic steps should be used for parameter search. If true, then (``value``, ``lo``, ``hi``) must all have log units.
        - ``LaTeX``: LaTeX string for this parameter to be used for plots.
    - ``ext_module``: a string which is the name of the module implementing the external :math:`\chi^2` calculation. The associated .py file should be in the Python path.
    - ``ext_class``: a string denoting the class name in the external module implementing the external :math:`\chi^2` calculation. It will be instantiated once, at the time the config file is read.
    - ``ext_class_args``: a dict holding the class parameters, can be empty (``{}``). Example: if ``ext_class_args: {arg1:47, arg2:"val2"}`` and ``ext_class: Chi2``, it will be instantiated via ``Chi2(arg1=47, arg2="val2")``.
    - ``ext_chi2``: a string which is the name of the ``ext_class`` method returning :math:`\chi^2` as a single ``float``. In DYNAMITE, it will be called after weight solving, passing the entire current parameter set as a dict, e.g. ``{a-bh:<val1>, m-bh:<val2>, ... ,p-stars:<valx>, q-stars:<valy>, ..., par1-gas:<valg1>, par2-gas:<valg2>}``.

``component types``
^^^^^^^^^^^^^^^^^^^^

The following types of component are available, listed with their parameters:

- ``TriaxialVisibleComponent``, a triaxial ellipsoid with surface density specified as an MGE,
    - ``p``: intrinsic axis ratio B/A (i.e. intermediate-to-major), where :math:`0<p<1`
    - ``q``: intrinsic axis ratio C/A (i.e. minor-to-major), where :math:`0<q<p`
    - ``u``: ratio between 2D observed and 3D intrinsic Gaussian widths of the MGE, i.e. :math:`\sigma_{2D}/\sigma_{3D}`
    - additionally, you must specify `observed data`_ for this component
- ``Plummer``
    - ``a``: scale length [arcsec]
    - ``m``: mass [:math:`M_\odot`]
- ``NFW``
    - ``c``: concentration parameter [:math:`R_{200}` / NFW-scale-length]
    - ``f``: dark matter fraction [:math:`M_{200}` / total-stellar-mass]
- ``NFW_m200_c``, an NFW halo with mass-concentration from `Dutton & Maccio (2014) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.441.3359D/abstract>`_
    - ``f``: dark matter fraction [:math:`M_{200}` / total-stellar-mass]
- ``Hernquist``
    - ``rhoc``: central density [:math:`M_\odot/\mathrm{km}^3`]
    - ``rc``: scale length [km]
- ``TriaxialCoredLogPotential``, see e.g. Binney \& Tremaine second edition p.171
    - ``p``: intrinsic intermediate-to-major axis ratio, where :math:`0<p<1`
    - ``q``: intrinsic minor-to-major axis ratio, where :math:`0<q<p`
    - ``Rc``: core radius [kpc]
    - ``Vc``: circular velocity for :math:`r>>R_c` [km/s]
- ``GeneralisedNFW`` from `Zhao (1996) <https://ui.adsabs.harvard.edu/abs/1996MNRAS.278..488Z/abstract>`_
    - ``c``: concentration parameter [:math:`R_{200}` / NFW-scale-length]
    - ``Mvir``: virial mass :math:`M_{200}` [:math:`M_\odot`]
    - ``gam``: AKA gamma, the inner logarithmic density slope, must be :math:`\leq 1`
- ``Chi2Ext``
    - Choose any parameter names as you see fit (e.g., ``par1``, ``par2``,...)

.. note::
  currently, there are limited combinations of component types that are valid. This is to ensure compatibility with the Fortran implementation of the orbit integrator. Later implementations may offer more flexibility. The current valid combinations of components are:

  - one ``Plummer`` component
      - representing the black hole
      - the scale length ``a`` should be fixed to some arbitrarily small value
  - one ``TriaxialVisibleComponent`` component
      - representing the stars
  - either no dark halo or exactly one out of [``NFW``, ``NFW_m200_c``, ``Hernquist``, ``TriaxialCoredLogPotential``, ``GeneralisedNFW``]
      - representing the dark halo
  - either no or one ``Chi2Ext`` component
      - calculating the external additive :math:`\chi^2`

.. _observed_data:

``observed data``
^^^^^^^^^^^^^^^^^^^^

The ``TriaxialVisibleComponent`` represents the galaxy's stars, and therefore has associated observations. You must specify the following entries with filenames for observed data:

- ``TriaxialVisibleComponent``
    - ``mge_lum``: string, filename for the MGE of the projected luminosity density, with intensity units of :math:`L_\odot \mathrm{pc}^{-2}`.
    - ``mge_pot``: string, filename for the MGE of the projected mass density, with intensity units of :math:`M_\odot \mathrm{pc}^{-2}`. If you assume that stellar-mass follows stellar-light, then the files ``mge_lum`` and ``mge_pot`` will be identical.
    - ``kinematics``
        - ``name_of_the_kinematic_set``: a descriptive name, best without spaces as it will be part of the kinematic plot file name.
            - ``type``: type of kinematics - either ``GaussHermite`` or ``BayesLOSVD``
            - ``weight``: float, weighting applied to this kinematic set in chi2 calculation; weights don't need to add up to 1.0.
            - ``datafile``: string, filename for the kinematics ECSV data file
            - ``aperturefile``: string, filename of the aperture file for this kinematic set
            - ``binfile``: string, filename of the bin file for this kinematic set
            - ``hist_width``: *optional*, float or 'default', the width (i.e. min. to max. value) of the velocity histogram for storing orbits. The default option is a width slightly wider than that of the observed kinematics.
            - ``hist_center``: *optional*, float or 'default', the center of the velocity histogram for storing orbits. The default option is 0.
            - ``hist_bins``: *optional*, int or 'default', the number of bins in the velocity histogram for storing orbits. The default option gives about 10 times better velocity sampling than the data.
        - ``name_of_next_kinematic_set`` (if any...)
            - ...

For more information on the input file formats, please refer to the :ref:`input_files` section of the Overview page.

``system_parameters``
=====================

This section is used for *global* parameters of the system i.e. those which are unrelated to any particular component.

Currently there is only one such parameter, ``ml``, which is a scale factor for the **total mass** of the system. Note that this scales the mass of **every** component of the system i.e. not just the stellar component (despite the acronym ``ml`` resembling *mass-to-light*). This is a time-saving trick: by scaling the total mass of the system, we are able to cheaply re-use orbit-libraries by re-scaling their velocity axes.

Care must be taken when interpreting mass parameters for models with different ``ml``. For example, say the system has a ``GeneralisedNFW`` component with ``Mvir=100`` but the system's ``ml`` parameter is equal to 2. The ``GeneralisedNFW`` would therefore *actually* represent a halo with mass ``Mvir=200``. Further note that the ``NFW`` component is parameterised with a mass *fraction* ``f`` rather than an absolute mass, and this fraction does **not** need to be re-scaled by ``ml``.

Specifying the ``ml`` parameter in the configuration file follows the same pattern as other parameters,

- ``system_parameters``
    - ``ml``
        - ``fixed``: Boolean, whether ``ml`` is to be kept fixed
        - ``value``: an initial value for ``ml``
        - ``par_generator_settings``: settings controlling parameter search (can be omitted if ``fixed=True``). Note that if these settings are given, then ``value`` must be consistent with ``lo`` and ``hi``.
            - ``lo``: minimum value
            - ``hi``: maximum value
            - ``step``: initial step size for parameter search
            - ``minstep``: minimum allowed stepsize for this parameter
        - ``logarithmic``: Boolean, whether logarithmic steps should be used for parameter search. If true, then (``value``, ``lo``, ``hi``) must all have log units
        - ``LaTeX``: LaTeX format string for this parameter to be used for plots, e.g. in axis labels.


``orblib_settings``
=====================

This section is used for settings relevant for the calculation of orbit libraries.

.. note::
  The size of the orbit library is controlled by 4 parameters: :math:`(n_E, n_{I2}, n_{I3})` and ``dithering``. The parameters :math:`(n_E, n_{I2}, n_{I3})` are the grid-dimensions in the three *integrals-of-motion* used for generating orbit initial conditions. Each initial-condition is used three times: once to seed a *box-orbit*, and twice to seed *tube-orbits* with opposing senses of rotation. The parameter ``dithering`` then seeds a *mini-grid* of orbits around each set of initial conditions, of size ``dithering``:math:`^3`. The total number of orbits in the library is thus

  .. math::

    \text{total number of orbits} = 3 \; n_E \; n_{I2} \; n_{I3} \; \mathrm{(dithering)}^3

  .. table:: Recommendations for orbit library sizes
     :width: 95%
     :widths: auto

     =====  =====  =====  ===========  ================================================
       nE    nI2    nI3    dithering    Example use
     =====  =====  =====  ===========  ================================================
       5      4      3         1        Test-orbit library.
                                        A good orbit library for fast checking whether DYNAMITE runs, but it is too small to be used for scientific analyses
      21     10      7         5        Good orbit library for CALIFA and ATLAS3D-like data quality.
                                        Examples: CALIFA (`Zhu et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.473.3000Z/abstract>`_),
                                        SAMI (`Santucci et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...930..153S/abstract>`_),
                                        ATLAS3D (`Thater et al. 2023b <https://arxiv.org/abs/2305.09344>`_), MANGA
      still bigger orbit libraries      MUSE-like data: e.g., `Poci et al. 2021 <https://arxiv.org/pdf/2102.02449.pdf>`_,
                                        `Ding et al. 2023 <https://arxiv.org/abs/2301.05532>`_,
                                        `Thater et al. 2023 <https://arxiv.org/abs/2304.13310>`_
     --------------------------------  ------------------------------------------------
     =====  =====  =====  ===========  ================================================

- ``orblib_settings``
    - ``nE``: integer, size of grid in integral-of-motion :math:`E`
    - ``nI2``: integer, size of grid in second integral-of-motion :math:`I_2` (similar to :math:`L_z`). Must be at least 4.
    - ``nI3``: integer, size of grid in third integral-of-motion :math:`I_3`
    - ``dithering``: integer, size of mini-grid of orbits around each set of initial conditions
    - ``logrmin``: log10 of minimum orbit radius in arcsecs
    - ``logrmax``: log10 of maximum orbit radius in arcsecs
    - ``random_seed``: integer, used for stochastically blurring orbit library by the PSF. Any value :math:`\leq 0` gives a stochastic seed.
    - ``quad_nr``: integer, sampling of grid recording the intrinsic moments in :math:`r`, default if missing: 10
    - ``quad_nth``: integer, sampling of grid recording the intrinsic moments in :math:`\theta`, default if missing: 6
    - ``quad_nph``: integer, sampling of grid recording the intrinsic moments in :math:`\phi`, default if missing:  6

The following settings must also be set in the configuration files but have *typical* values which should generally be sufficient and should not be changed,

- ``orblib_settings``
    - ``orbital_periods``: integer, typical 200, the number of orbital periods to integrate orbits
    - ``sampling``: integer, typical 50000, number of points to sample for each orbit in the meridional plane
    - ``starting_orbit``: integer, typically 1, the index of which  orbit to start integrating orbits
    - ``number_orbits``: integer, the number of orbits to integrate, if -1 then integrate all orbits
    - ``accuracy``: typical ``1.0d-5``, the accuracy of the orbit integrator

``weight_solver_settings``
==========================

Settings relevant for solving for orbital weights.

.. note::
  If any kinematic set has type ``BayesLOSVD``, then the ``weight_solver_settings`` must have type ``NNLS``

- ``weight_solver_settings``
    - ``type``: string, one of ``LegacyWeightSolver`` to use Fortran implementations of Lawson and Hanson non-negative least-squares algorithm, or ``NNLS`` to use Python implementations
    - ``nnls_solver``: options depend on the ``type`` selected. If
        - ``type = LegacyWeightSolver`` then set ``nnls_solver : 1``
        - ``type = NNLS`` then ``nnls_solver`` can be one of the strings,
            - ``scipy`` to use the `scipy NNLS function <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html>`_
            - ``cvxopt`` to use an implementation using the `CVXOPT <https://cvxopt.org/>`_ package
    - ``lum_intr_rel_err``: float, typical 0.01, the systematic error (fraction) applied to the intrinsic luminosity constraint
    - ``sb_proj_rel_err``: float, typical 0.01, the systematic error (fraction) applied to the projected surface brightness constraint
    - ``CRcut``: Boolean, default False, whether to use the ``CRcut`` solution for the counter-rotating orbit problem. See `Zhu et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.473.3000Z/abstract>`_ for more details.

If any kinematics have of type ``GaussHermite``, the following additional settings are needed.

- ``weight_solver_settings``
    - ``number_GH``: integer, the highest order kinematics to be used when solving for orbital weights. Note that this can be different from the order of the input data you provide. If ``number_GH`` is lower than in the data, then higher order kinematics are ignored while weight solving. Alternatively, if ``number_GH`` is higher than in the data, then we (fictitiously) assume that the higher-order kinematics were observed to be zero, with a *nonzero* systematic error that must be specified in the ``GH_sys_err`` setting. The latter option can be considered as a form of regularisation, penalising solutions where higher-order kinematics (although unobserved) reach unrealistically high values.
    - ``GH_sys_err``: a string of floats, must contain at least ``number_GH`` entries. These are systematic errors applied to ``V``, ``sigma``, ``h3``, ..., ``hN``. During weight solving, these systematic errors are added in quadrature to the random errors which you provide in the data file. If ``number_GH`` is larger than the kinematic order of the observed data, then the corresponding systematic errors must be > 0 and can be interpreted as a typical value for higher order kinematics; models with higher-order kinematics which exceed this typical value will be penalised.

If any kinematic set has type ``BayesLOSVD``, then the ``weight_solver_settings`` must have type ``NNLS``, and no additional settings are required.

If DYNAMITE shall recover from an unsuccessful weight solving attempt, the following option can be used:

- ``weight_solver_settings``
    - ``reattempt_failures``: if True, DYNAMITE will use a model's existing orblibs from an earlier run to reattempt weight solving.

``parameter_space_settings``
============================

Settings relevant for parameter search.

- ``parameter_space_settings``
    - ``generator_type``: string, specifying which algorithm to use for parameter search. Note that all generator types will exclude invalid or already-executed parameter combinations by default. The different options are:
        - ``GridWalk``: Start at the initial point. Start the iteration: (i) find the model with the minimum :math:`\chi^2`, (ii) for each free parameter, seed new models by independently take a step :math:`\pm 1` of size ``step`` (cartesian grid in one step size, so if 2 parameters are free, 8 new models will be created). Repeat until :math:`\chi^2` is improved by less than min_delta_chi2. This may result in a large number of models.
        - ``LegacyGridSearch``: Start at the initial point. Start the iteration: (i) find all models with :math:`|\chi^2 - \chi_\mathrm{min}^2|` within the threshold (specified with ``threshold_del_chi2_XXX``), (ii) for each model within the threshold, seed new models by independently take a step :math:`\pm 1` of size ``step`` (i.e. as done for ``GridWalk``).  If no new models are seeded at the end of an iteration, then divide all parameter stepsizes by two till their specified ``minstep`` are reached.
        - ``FullGrid``: Create a *full* grid, i.e. a Cartesian grid in all free parameters, with bounds ``lo/hi`` and stepsize ``step``. **Warning**: If several (>3) parameters are free, this will result in a large number of models.
    - ``which_chi2``: string, specifies which :math:`\chi^2` value to consider when generating new parameters, must be one of the following:
        - ``kinchi2``: this includes contributions from only the kinematics. If ``GaussHermite`` kinematics are used then this is includes terms from all Hermite coefficients :math:`h_1, h2, h3, ..., h_N`. If ``BayesLOSVD`` kinematics are used, then this includes contributions from all LOSVD bins.
        - ``chi2``: this includes contributions from the observed surface density, de-projected 3D density, and kinematics (as specified above).
        - ``kinmapchi2``: the :math:`\chi^2` directly calculated from the ``GaussHermite`` kinematic maps (not available for ``BayesLOSVD`` kinematics).
    - ``generator_settings``: if ``generator_type = LegacyGridSearch``, then one of the following two settings must be set. These are the :math:`|\chi^2|` thresholds used for in ``LegacyGridSearch``,
        - ``threshold_del_chi2_abs``: an absolute :math:`|\chi^2|` threshold
        - ``threshold_del_chi2_as_frac_of_sqrt2nobs``: a threshold given as a fraction of :math:`\sqrt{2N_\mathrm{obs}}` where :math:`N_\mathrm{obs}` is the total number of kinematic observations, which is equal to the number of spatial apertures multiplied by (i) ``number_GH`` if ``GaussHermite`` kinematics are used, or (ii) the number of LOSVD bins if ``BayesLOSVD`` kinematics are used.
    - ``stopping_criteria``: all of the following must be specified. If any of the criteria are met, then the parameter generation will stop:
        - One of ``min_delta_chi2_abs`` or ``min_delta_chi2_rel`` must be set: float, absolute or relative tolerance for ending the parameter search. If an iteration does not improve the minimum chi2 by this threshold, no new iteration will be performed.
        - ``n_max_mods``: int, maximum number of models desired
        - ``n_max_iter``: int, maximum number of iterations to be run. The iteration a model was created in is listed under the ``which_iter`` column of the ``all_models`` table, and these are indexed from ``0,... n_max_iter-1``. The ``n_max_iter`` setting controls the total *cumulative* number of iterations to run i.e. if you specify ``n_max_iter=10`` and there are existing models which ``which_iter=9``, then no new iterations will be run. Note that the first two iterations are always run together i.e. whether you specify ``n_max_iter=1`` or ``n_max_iter=2``, iterations 0 and 1 will both be run.

``io_settings``
=====================

Settings specifying the location of input and output directory names. Paths are relative to the current working directory, and can be given with or without trailing slash::

    io_settings
        input_directory: "input_files/"     # directory holding input data
        output_directory: "output/"         # directory (will be created) for output
        all_models_file: "all_models.ecsv"  # filename for the summary file of models run so far

``multiprocessing_settings``
============================

Settings for multiprocessing. Models can be evaluated in parallel, with the number of parallel processes specified by the ``ncpus*`` settings::

  multiprocessing_settings:
      ncpus: 4                              # integer or string 'all_available' (default: 'all_available')
      ncpus_weights: 4                      # int or 'all_available', optional (default: ncpus), not used by all iterators
      orblibs_in_parallel: True             # calculate tube and box orbits in parallel (default: False)
      modeliterator: 'SplitModelIterator'   # optional (default: 'ModelInnerIterator')

Due to very different CPU and memory consumption of orbit integration and weight solving, there are two different settings: while orbit integration will use ``ncpus``, weight solving will use ``ncpus_weights`` parallel processes, with ``ncpus`` ≥ ``ncpus_weights`` in general. Note that ``ncpus_weights`` will default to ``ncpus`` if not specified. Currently, only the ``SplitModelIterator`` model iterator and recovering from an unsuccessful weight solving attempt (``reattempt_failures=True``) use the ``ncpus_weights`` setting.

If ``orblibs_in_parallel`` is set to ``False``, DYNAMITE will first integrate the tube orbits and then the box orbits. If it is set to ``True``, the tube and box orbits will be integrated in parallel, which will use 2 parallel processes per model.

If ``ncpus : 'all_available'`` or ``ncpus_weights : 'all_available'`` is set, then DYNAMITE automatically detects the number of available cpus :math:`N_\mathrm{CPU}` for parallelisation and will set ``ncpus`` = ``ncpus_weights`` = :math:`N_\mathrm{CPU}`.

Important performance hint:

- Most ``numpy`` and ``scipy`` implementations are compiled for shared-memory parallelism (e.g., involving blas/openblas). This can be verified by inspecting the ``MAX_THREADS`` values in the output of ``numpy.__config__.show()`` and ``scipy.__config__.show()``, respectively. The number of threads to be used by ``numpy`` and ``scipy`` can be limited by setting the environment variable ``OMP_NUM_THREADS`` to the desired value before executing DYNAMITE.
- Recommendation: ``OMP_NUM_THREADS=n`` with ``ncpus * n`` ≤ :math:`N_\mathrm{CPU}` if ``orblibs_in_parallel`` is set to ``False`` and ``ncpus * n`` ≤ :math:`\frac{1}{2}\,N_\mathrm{CPU}` if ``orblibs_in_parallel`` is set to ``True``.

``legacy_settings``
=====================

Location of the *legacy* Fortran programs::

  legacy_settings:
      directory: "default"  # or an alternative directory

If ``default``, then the Fortran programs created during installation are used. Can be set to an alternative directory if required.
