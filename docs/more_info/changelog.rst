.. _changelog:

****************
Change Log
****************

- Improvement: more robust directory naming when failed models are present in the all_models table
- Bugfix: don't crash with pops data when using ModelInnerIterator.
- Bugfix: don't crash when continuing a run with just one valid model.
- Improvement: added counter-rotating components to decomposition.
- Bugfix: fix a bug that prevents Bayes LOSVD kinemtic maps from being created more than once
- New feature: a dark halo is no longer mandatory (models can consist of either zero or one dark halo component)
- Improvement: Reduced disk space requirements and performance by splitting orbit library files while maintaining backward compatibility.
- Bugfix: avoid DYNAMITE crashing because it gets stuck in a model directory due to a Fortran error
- Improvement: better error messages to improve debugging with multiprocessing
- Improvement: minor improvements of tutorials 2 and 3

Version: 4.3
================

- Improvement: Updated requirements to Python 3.9 or later (Python 3.8 end of life was Oct 7, 2024) and numpy<1.27.0, scipy>=1.11,<1.12 to avoid Python NNLS freeze
- Improvement: made DYNAMITE compatible with numpy 2.0.0 and matplotlib 3.9.0 (removed use of deprecated features)
- Bugfix: fix a bug that prevents Bayes LOSVD kinematic maps from being created more than once
- Improvement: misc. improvements in documentation and tutorials.
- Improvement: Sampling of grid recording the intrinsic moments (in r, theta, phi) can optionally be defined in the config file's orblib settings.
- Improvement: The installation procedure has been changed to make DYNAMITE compatible with Python 3.12. Installing and uninstalling is now done using pip.
- Improvement: added a new tutorial notebook ``7_orbital_distributions.ipynb`` which takes a closer look at orbit distributions.
- Improvement: the ``weight`` attribute of kinematics is now officially DEPRECATED as it has always been ignored by DYNAMITE.
- Improvement: DYNAMITE now checks for nan values in the kinematics and mges when first reading the data
- Improvement: prevent DYNAMITE from crashing if NNLS weight solving fails.
- Improvement: the Gauss Hermite kinematic maps new parameter value `cbar_lims='user'` allows user-defined velocity and velocity dispersion limits (see `Plotter.plot_kinematic_maps()`).
- Improvement: reduced the main memory requirements of the Python NNLS solvers.
- Bugfix: fix a crash when creating the BayesLOSVD kinematics file in rare cases where the completed bins were determined incorrectly.
- New feature: ``data_prep/generate_kin_input.py`` implements reading NIFS kinematics with an arbitrary number of GH moments.
- Improvement: improved checks and error messages for velocity and spatial bin input data inconsistencies.
- Improvement: save disk space by cleaning up decompressed files after a crash and removing unused legacy file nn_orbmat.out after solving.
- Improvement: stability fix in MGE: if q>0.9999 it will be set to 0.9999 (before, it was 0.99999).
- Bugfix: the chi2 plot now shows correct axis ticks for log quantities.
- Bugfix: fixed colorbar overlap with x-axis in the chi2 plot if only two parameters are varied and added label to chi2 plot colorbar.
- Improvement: ``LegacyWeightSolver`` is now DEPRECATED and will be removed along with GALAHAD in a future version of DYNAMITE. Use weight solver ``type: "NNLS"`` instead if you can.

Version: 4.2
================

- Improvement: if ``number_GH`` in the config file is larger than the kinematic order of the observed data, then DYNAMITE ensures that the corresponding systematic errors are > 0.
- Bugfix: fixed a bug in the kinematics errors (affects NNLS solves).
- New feature: Gauss-Hermite kinematic maps can now be plotted for any number of Gauss-Hermite coefficients.
- Improvement: removed broken link from tutorial 2 and added some data preparation comments to tutorials 1 and 2
- Bugfix: fixed crash when different kinematics had different numbers of PSF components.
- Bugfix: fixed a bug in retrofitting kinmapchi2 in old all_models tables.
- Improvement: removed deprecated silent option from config reader.
- Improvement: the Plotter's new optional argument ``dpi`` (default: 100) allows to change the resolution of all saved figures except the kinematic maps (always 300 dpi).
- Improvement: the beta plots now work for all implemented weight solvers.

Version: 4.1
================

- Improvement: calculation of kinmapchi2 now aligns with number_GH in config file's weight_solver_settings
- Bugfix: fixed crash when the number of GH coefficients a kinematics file does not match number_GH in config file's weight_solver_settings
- Improvement: The bash test script ``dev_tests/test_notebooks.sh`` executes all tutorial notebooks for testing a valid DYNAMITE installation
- Improvement: updated tutorial notebooks
- Bugfix: Re-enable support for directly instantiating a Model object (bypassing ModelIterator) if the all_models table is empty. Only recommended for testing.

Version: 4.0
================

- New feature: kinmapchi2 (directly calculated from the kinematic maps) is now also available for the python NNLS solver
- New feature: added support for bar/disk decomposition
- New feature: added support for getting intrinsic model moments for both Gauss Hermite and a BayesLOSVD models
- Improvement: Eliminated unused position_angle system attribute from the configuration file (the angle is read from aperture.dat)
- Improvement: DYNAMITE can now be built without GALAHAD (LegacyWeightSolver will not be available then)
- Improvement: plotting gh kinematic maps is more efficient and now works for all weight solvers
- New feature: New parameter generator SpecificModels generates and runs a predefined list of models or models resulting from a cartesian product of parameter values
- Improvement: the orbit plot (Plotter.orbit_plot) now works for all implemented weight solvers
- Bugfix: fixed a bug that under certain circumstances prevented the staging files from being deleted after a successful iteration
- Bugfix: removed logging from the list of requirements because it is in the Python standard library
- New feature: added a new DYNAMITE module analysis, its class Decompostion creates decomposition plots
- Improvement: Changed the ml directory name format to '05.2f' so all model directory names have the same length
- New feature: Added AllModels.remove_unused_orblibs() utility method to free up disk space
- New feature: Added a new method AllModels.make_best_models_table() that creates a table of the best models (best n models or models within a chi2-threshold of the best) and saves it to disk
- Bugfix: If reattempt_failures=False, in certain cases it could occur that orblibs of successful models were deleted
- Bugfix: Fixed a bug related to a nonexistent model directory if a crash occurs between the parameter generator adding a model and starting to solve it
- Improvement: Dynamite will no longer crash upon Legacy Fortran errors (except when executables are not found), but issue warnings and assign nan to the affected chi2 values
- Improvement: When executing a dummy run (do_dummy_run==True), model_iterator will set both kinchi2 and kinmapchi2 to nan (instead of zero)
- Improvement: DYNAMITE will retrofit existing all_models tables with the new column kinmapchi2 and calculate its values for existing models whenever possible
- New feature: chi2 can now be directly calculated from the kinematic maps when using the LegacyWeightSolver via which_chi2: "kinmapchi2"
- Improvement: when instantiating the Configuration object, the user can now specify the name of the logfile (several options), avoiding log conflicts with multiple DYNAMITE runs in the same directory
- Bugfix: Fixed a bug that may cause a crash in case a parameter does not have a minstep value
- Improvement: DYNAMITE will catch and correct the erroneous parameter generator setting minstep>step by setting minstep=step for non-fixed component parameters
- Bugfix: Fixed a bug that may occur in the parameter generators (ensures that DYNAMITE creates all possible models)
- Improvement: now the models of the first two iterations are computed together, better utilizing parallel computing
- Bugfix: included cmasher in the list of required packages
- Bugfix: reattempt_failures will no longer result in an error if multiple to-delete models share the same orblib or the orblib directory does not exist
- Improvement: made DYNAMITE compatible with more Linux distributions
- Improvement: update publication list
- Bugfix: fixed wrong version number and copyright year in documentation

Version: 3.0
================

- Improvement: DYNAMITE now works with newer versions of Astropy. The new requirement is astropy>=5.0.4
- New feature: Integrate tube and box orbits in parallel by setting the multiprocessing option orblibs_in_parallel
- New feature: Added support for new dark halo component type NFW_m200_c (fixed m200_c relation)
- New feature: The Configuration object parameter reset_existing_output will delete previously existing data and create a new output directory tree
- Improvement: The presence of datfil/orblib.dat.bz2 and datfil/orblibbox.dat.bz2 is now a more reliable indicator for existing orblibs. In the past, a crash may have resulted in corrupt bz2 files.
- New feature: add new data-preparation method `BayesLOSVD.center_v_systemic`
- New feature: Each model writes a file model_done_staging.ecsv upon completion. After a crash, DYNAMITE will update the all_models table with the completed models' data and delete any "all_done==False" models
- New feature: New weight_solver_settings option reattempt_failures for reattempting failed weight solving when an orbit library already exists
- Improvement: For better tracking, each model folder holds a copy of the config file now (instead of saving the config file in the output folder)
- New feature: The new model iterator SplitModelIterator calculates orbit libraries and weights consecutively with independently adjustable number of threads
- Improvement: Cleaned up the legacy_fortran folder and the makefile in it, unused orbgen.f90 and partgen.f90 moved to subfolder
- Bugfix: Fixed a bug that on rare occasion caused an error when updating the timestamp entry when continuing an aborted run
- Implement the correction to orbit mirroring introduced in `Quenneville et al 2021 <https://arxiv.org/abs/2111.06904>`_
- Implement kinematic maps for BayesLOSVD data

Version: 2.0
================

- New feature: Wherever appropriate, the configuration object is now passed to instantiated classes like Plotter, AllModels, Model, the weight solvers, and LegacyOrbitLibrary. This changes the DYNAMITE user interface! Please refer to the tutorials and ``dev_tests/`` scripts on how to use.
- Bugfix: fixed sorting of the chisquare values in chi2plot so that the best-fit value is plotted last and always visible
- New feature: upon reading mge data, q values too close to 1 are set to q=0.99999 for numerical stability
- Updated: Require astropy v4.2 due to ecsv file compatibility problems with later versions
- Improvement: Make sure DYNAMITE stops with an error if a legacy Fortran executable does not succeed, even if its return code is zero
- Bugfix: Fixed a bug preventing negative values of logarithmic parameters
- New feature: the number of configuration file backups can be better controlled by ``config_reader.Configuration.backup_config_file(...)`` options
- Improvement: The bash test script ``dev_tests/test_all.sh`` executes a grid of test scenarios (different base scripts with different parameter generators and weight solvers) either locally or via Slurm
- Improvement: Greatly improved performance of the chi2 plot
- Bugfix: Fixed a bug in the path in ``model.Model.get_model_directory()``
- New feature: The new method ``model.AllModels.get_n_best_models(...)`` returns the ``n`` best models based on their ``chi2``/``kinchi2`` values
- New feature: The new method ``model.AllModels.get_mods_within_chi2_thresh(...)`` returns all models within a given ``chi2``/``kinchi2`` threshold
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
