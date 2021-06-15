.. _changelog:

****************
Change Log
****************

Version:
================

- New feature: the new method ``model.AllModels.get_n_best_models(...)`` returns the ``n`` best models based on their ``chi2``/``kinchi2`` values
- New feature: the new method ``model.AllModels.get_mods_within_chi2_thresh(...)`` returns all models within a given ``chi2``/``kinchi2`` threshold
- Updated: All tests in ``dev_tests/`` now use ``kinchi2`` rather than ``chi2``

Version: 1.0
================

- New feature: Added Galahad compilation script that auto-magically downloads and installs the latest galahad + it's dependencies
- New feature: Added a script for the preparation of the kinematic data and a tutorial
- New features: Added Bayes LOSVD solver and a tutorial
- New feature: Added gridSearch that searches in a regular grid for the bestfit parameters
- New feature: In addition to the NFW profile, DYNAMITE can fit now a generalised NFW, Hernquist and a Triaxial cored log potential dark matter profile. The type is chosen in the “dh” component of the config file.
- New feature: All plotting routines from schwpy are implemented in DYNAMITE now
- New feature: Added multiprocessing such that DYNAMITE can run multiple models simultaneously. The keyword “multiprocessing_settings: ncpus:” is added in the config file
- Improvement: New (python-based) NNLS solvers are added. The type can be chosen in “weight solver”
- Improvement: Multiple kinematics data set can be fitted simultaneously
- Improvement: Changed paramsb and parameter file to “parameters_lum” and “parameters_pot” to avoid confusion. The mass mge and the lum mge can be different now and are added separately in the config file
- Improvement: Changed the model directory names to avoid directory naming inconsistencies in the future
- Improvement: Logging added
- Improvement: The DYNAMITE scrips no longer change the system path
- Improvement: Added “validate_parset" to the system and its components to avoid incorrect use of DYNAMITE
- Improvement: Relative/absolute stopping criteria in LegacyGridSearch and GridWalk
- Improvement: Option for threshold_del_chi2 to be given as fraction of sqrt(2*n_obs)
- Improvement: “get_orbit_ics” and “get_orbit_library” are split now in LegacyOrbitLibrary
- Updated: Installation guide and documentation were updated
- Updated: Replaced Plotbin4dyn with the latest version from plotbin (https://pypi.org/project/plotbin/)
- Updated: We added a randomNumberGenerator to get reproducible orbit libraries. This number called “random_seed” is included in the config file. Users should put this number to -1
- Bugfix: A galaxy with position angle of 0 does not cause error when reading in the config file anymore
- Bugfix: Fix the stars component bug: stellar component was called system.cmp_list[2] before and relied on the stars being the third component. Now this component is called “stars”
- Bugfix: The plotting did not work correctly in VSC where X11 does not work. We therefore put the matplotlib backend to “Agg”
- Bugfix: Removed unused import statements and code clean-up
