.. _overview:

******************
Overview
******************

This page gives a descriptive overview of everything you will need to run orbit-based models in DYNAMITE:

1. `Directory Structure`_
2. `Input Files`_
3. `Configuration File`_
4. `The Main Script`_
5. `Plotting`_
6. `Multiprocessing + Slurm Submission`_
7. `Managing output`_
8. `Logging output`_

The tutorials show an example of running DYNAMITE from start to finish - this would also be a great place to start getting acquainted with the code.

Directory Structure
===================

Here is an example of the directory structure needed to run DYNAMITE. After installing DYNAMITE on your system, you can run the code from any location. You should create a ``main_directory`` with the following structure::

      | main_directory
      | ├── input_files     # contains all input data-files
      | │   ├── input_file_1.txt
      | │   ├── input_file_2.txt
      | │   └── ...
      | ├── config_file.yaml
      | ├── main_script.py
      |

To run the script, you must be in ``main_directory`` to execute the command ::

    python main_script.py

After running the script, the following directories/files will be created::

  | main_directory
  | ├── input_files     # contains all input data-files
  | │   ├── input_file_1.txt
  | │   ├── input_file_2.txt
  | │   └── ...
  | ├── config_file.yaml
  | ├── main_script.py
  | ├── dynamite.log    # a log file
  | ├── output
  | │   ├── models/     # output model directory
  | │   ├── plots/      # output plot directory
  | │   ├── all_models.ecsv        # summarises all models run so far
  |

Subsequent runs of scripts from the main directory (e.g. after you have altered configuration settings) will not change this directory structure.
Instead, each run will add more output to the existing output directory.
To keep track of your configuration settings, a copy of the configuration file will be created in the output model directories e.g. ``output/models/orblib_XXX_YYY/mlZZZ/config_file.yaml``.

Note: DYNAMITE can also be run interactively, e.g. from a Jupyter notebook, but this must be launched from ``main_directory``.

.. _input_files:

Input Files
===================

The following input files are required::

  | main_directory
  | ├── input_files
  | │   ├── mge.ecsv          # the MGE of stellar surface density
  | │   ├── kinematics.ecsv   # file of kinematic data
  | │   ├── bins.dat          # info about binning of kinematics
  | │   ├── aperture.dat      # more info about binning of kinematics
  |

The exact filenames can be freely defined in the configuration file as described in the :ref:`observed_data` section of the Configuration page.

The Multi Gaussian Expansion (MGE) in ``mge.ecsv`` describes the galaxy's 2D surface-brightness distributon. To generate this, fit an MGE to a photometric image e.g. using `mge <http://www-astro.physics.ox.ac.uk/~mxc/software/#mge>`_. The data must be given as a table in `Astropy ECSV format <https://docs.astropy.org/en/stable/api/astropy.io.ascii.Ecsv.html>`_, with columns

- ``I``: peak surface brightness values of the MGE Gaussians describing the surface brightness of the tracer population for which the kinematics is derived. Units: ``L_sun/pc^2`` (solar luminosities per ``parsec^2``).
- ``sigma``: dispersion of the MGE Gaussians describing the distribution of the kinematic-tracer population. Units: ``arcsec`` (arcseconds).
- ``q``: observed axial ratio q of the MGE Gaussians describing the distribution of the kinematic-tracer population.
- ``PA_twist``: observed position angle (psi=PA_twist) of the MGE Gaussians describing the distribution of the kinematic-tracer population. Units: ``deg`` (degrees).

It is also possible to provide two separate MGE's for the surface-brightness and surface mass-density (see the :ref:`observed_data` section of the Configuration page for details.

Two types of kinematic are supported: tables of Gauss Hermite expansion coefficients, or histogrammed LOSVDs output by `BayesLOSVD <https://github.com/jfalconbarroso/BAYES-LOSVD>`_.
These must be in the form of `Astropy ECSV files <https://docs.astropy.org/en/stable/api/astropy.io.ascii.Ecsv.html>`_ e.g., ``kinematics.ecsv``. The files ``aperture.dat`` and ``bins.dat`` contain information about the spatial binning of your kinematic data. Convenience functions are provided for creating and converting some standard kinematic data files, and examples demonstrating these can be found in the tutorials.
Note that the kinematics need to be centered at the center of the MGE.

The file ``aperture.dat`` file contains the spatial extent in arcseconds, the angle (in degrees ) ``90 - position_angle``, and size of the grid in pixels::

  #counter_rotation_boxed_aperturefile_version_2
        min_x               min_y
        extent_x            extent_y
        90 - position_angle
        npix_x              n_pix_y

To clarify the quantities in this file: the rectangular aperture in arcseconds extends from ``min_x`` to ``max_x = min_x + extent_x`` and ``min_y`` to ``max_y = min_y + extent_y``, respectively. ``90 - position_angle`` is the angle in degrees measured counter clockwise from the galaxy major axis to the x-axis of the input data. Finally, ``npix_x`` and ``npix_y`` define the number of pixels along the x and y axes, respectively.

As ``aperture.dat`` is also read by legacy Fortran components of DYNAMITE, it is important that its first line (starting with `#`) is exactly as displayed above, otherwise DYNAMITE will crash.

The file ``bins.dat`` encodes the spatial (e.g. Voronoi) binning: specifically, one header line with the total number of pixels in the grid, followed by the bin ID of each pixel in the grid::

    #Counterrotation_binning_version_1
    no of pixels in grid
    ...

Note that also for this file the first line (starting with `#`) needs to be exactly like displayed above to avoid legacy Fortran errors. An exception to this is the typo ``#Counterrotaton_binning_version_1`` which will be accepted in ``bins.dat`` due to historical reasons.

Comments on kinematics
----------------------

Weight solver issues
^^^^^^^^^^^^^^^^^^^^

LegacyWeightSolver can't be used with BayesLOSVD - use weight-solver type NNLS.

In some cases, the weight solver ``type: "NNLS"`` ``nnls_solver: "scipy"`` may fail for some models if the SciPy version is >= 1.12. In such cases, it is recommended to use ``type: "NNLS"`` ``nnls_solver: "cvxopt"`` or to reinstall DYNAMITE with ``scipy<1.12`` in ``requirements.txt``.

Simultaneous fitting of multiple kinematic sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to simultaneously fit multiple sets of kinematics in DYNAMITE, mixing kinematics of different types is supported as is mixing kinematics with differing histogram settings. In that case, all input files should be placed in the input directory::

  | main_directory
  | ├── input_files
  | │   ├── mge.ecsv            # the MGE of stellar surface density
  | │   ├── kinematics_1.ecsv   # file of kinematic data 1
  | │   ├── bins_1.dat          # info about binning of kinematics 1
  | │   ├── aperture_1.dat      # more info about binning of kinematics 1
  | │   ├── kinematics_2.ecsv   # file of kinematic data 2
  | │   ├── bins_2.dat          # info about binning of kinematics 2
  | │   ├── aperture_2.dat      # more info about binning of kinematics 2
  |

The specific names of the files given here are just examples - you can specify the names you would like to use in the configuration file.

In case of Gauss Hermite kinematics, the individual kinematics' tables need to have the same number of expansion coefficients. In case your kinematics have different numbers of Gauss Hermite expansion coefficients, we recommend to augment the respecive tables with zero values for the additional coefficients and set the respective coefficients' errors to a large number (e.g., 0.3 or 0.5).

Keep in mind when mixing different kinematic data sets that their contribution to the orbit weights optimisation depends on the relative quantity and quality of the different data sets. For example, BayesLOSVD kinematics typically provides many more constraints than Gauss-Hermite kinematics, so that BayesLOSVD can dominate the weight solving and the Gauss-Hermite kinematics contribute little. An ad-hoc way to compensate for this is by adjusting the relative quality of the data sets, e.g., by reducing the errors on the Gauss-Hermite kinematics by a constant factor and so increase their relative contribution to the weight solving.

Configuration File
===================

All settings for running DYNAMITE can be controlled from a single configuration file. This specifies:

- the components of the gravitational potential
- the potential parameter values and ranges
- the type of kinematic data, e.g Gauss Hermite vs BayesLOSVD histograms
- settings for the orbit library, e.g. number of orbits
- the location of the input and output files
- the number of models you want to run

amongst others. More details can be found on the :ref:`configuration page <configuration>`.

The Main Script
======================

The main script should contain all of the DYNAMITE commands you wish to execute. This may change from run to run. This script must be executed from the ``main_directory``. Below are two common examples of what you may have in your main script.

To run a single Schwarzschild model ``main_script.py`` should be the following,

.. code-block:: python

   import dynamite as dyn

   c = dyn.config_reader.Configuration('config_file.yaml') # read the configuration fie
   parset = c.parspace.get_parset()                        # extract a parameter set from configuration
   model = dyn.model.Model(config=c, parset=parset)        # make a model object
   model.setup_directories() # make directory tree
   model.get_orblib()        # make an orbit library
   model.get_weights()       # find orbital weights

If you want to run a grid of models, ``main_script.py`` should be,

.. code-block:: python

  import dynamite as dyn

  c = dyn.config_reader.Configuration('config_file.yaml') # read the configuration fie
  smi = dyn.model_iterator.ModelIterator(config=c)        # create and run an iterative grid of models

You may have additional commands in the main script related to e.g. (i) plotting, (ii) multiprocessing, (iii) managing output, and (iv) logging. DYNAMITE provides functions for these four activities, described below.

Plotting
========

To make plots, you can use the Plotter object:

.. code-block:: python

  p = dyn.plotting.Plotter(config=c) # make the plotter object

Here we propose a few examples of the plots that can be done with this object. First, you can generate maps of the surface brightness, mean line-of-sight velocity, velocity dispersion, and higher order Gauss-Hermite moments. The figure produced will show the maps relative to the data in the first row, those relative to the best-fit model in the second row and residuals in the third row; it can be obtained by using:

.. code-block:: python

  p.plot_kinematic_maps(kin_set=0, cbar_lims='data') # the limits of color bars are based on the data values, and only the first kinematic set is plotted

To explore how the :math:`\chi^2` changes as a function of the parameters or of the model ID, you can use the following two functions, respectively:

.. code-block:: python

  p.make_chi2_plot(which_chi2='kinchi2', n_excl=50, figtype='.pdf') # saves a .pdf figure of the 'kinchi2' chisquare, excluding the first 50 models (burn-in)
  p.make_chi2_vs_model_id_plot(which_chi2='kinchi2') # saves a .png figure (default) of the 'kinchi2' chisquare as a function of the model ID

You can also plot the cumulative mass and the (intrinsic and projected) anisotropy profiles, out to a radius of 30 arcsec:

.. code-block:: python

  p.mass_plot(Rmax_arcs=30) # cumulative mass plot, saved as a .png file
  p.beta_plot(Rmax_arcs=30) # anisotropy plots, saved as .png files

These plots are made by considering only models close to the :math:`\chi^2` minimum, within a certain confidence level. You can decide which :math:`\chi^2` to use for this (``kinchi2`` is the recommended option), and what type of figure to produce, by specifying a file extension in the parameter ``figtype``.

To see how orbits are distributed in the best-fit model (or in a model of your choice, to be specified in the variable ``model`` when calling the function), you can use:

.. code-block:: python

  p.orbit_plot(Rmax_arcs=30) # orbit plot, saved as a .png file

In this case, ``Rmax_arcs`` indicates the upper radial limit for orbit selection, meaning that only orbits extending up to ``Rmax_arcs`` are plotted.

Finally, you can make a plot of the intrinsic flattening of your best-fit model:

.. code-block:: python

  p.qpu_plot(Rmax_arcs=30,figtype='.pdf') # triaxiality plot, saved as a .pdf file

In the examples above, the figures are created and saved automatically. If you want to make some changes into the appearance of the plots, you can use the fact that all the above functions return a ``matplotlib.pyplot.figure`` instance. For the figures to appear in the interactive mode, you first need to run the following line:

.. code-block:: python

  matplotlib.use('TkAgg')

and you can then proceed to make figures that you can modify as you prefer, for example:

.. code-block:: python

  fig = p.mass_plot(Rmax_arcs=30)

Please note that a copy of the figure as produced by DYNAMITE is always saved in the ``plots`` folder.


Multiprocessing + Slurm Submission
======================================

Different models can be run as separate processes. The number of processes which can be run simultaneously should be specified in the configuration file::

  multiprocessing_settings:
      ncpus: 4 # an integer or 'all_available'

If ``ncpus: 'all_available'`` is selected, the program will automatically detect the total number of disposable cpus.

If you use the Slurm job submission system on a cluster, then you must add a Python `shebang line <https://en.wikipedia.org/wiki/Shebang_(Unix)>`_ and any Slurm settings to the top of ``main_script.py`` e.g.

.. code-block:: python

  #!/bin/env python
  #SBATCH --job-name=my_dynamite_run
  #SBATCH --mem-per-cpu=50
  #SBATCH --qos={NAM OF YOUR QOS}
  #SBATCH -N {NUMBER OF NODES TO USE}
  #SBATCH --output="dyn_%j.out"
  #SBATCH --error="dyn_%j.err"

  import dynamite as dyn
  # etc ...

You can then submit this job as::

  sbatch main_script.py

So far we have not used job submission systems other than Slurm. If you need these, or have experience doing this yourself, please let us know and we will update the docs.

Note: multiprocessing is handled by the `pathos <https://pypi.org/project/pathos/>`_ module, specifically using ``pathos.multiprocessing.Pool``. This is very similar to the native Python ``multiprocessing.pool`` but can work with class methods as well as functions.

Managing output
===================

We provide utility functions to manage output, e.g. if you want to remove output from previous runs, change some configuration settings, before running again. These are methods of the configuration object, i.e.

.. code-block:: python

   import dynamite as dyn

   c = dyn.config_reader.Configuration('config_file.yaml') # read the configuration fie

where ``c`` has the following utility functions,

.. code-block:: python

  c.remove_existing_orblibs()
  c.remove_existing_orbital_weights()
  c.remove_existing_plots(...)
  c.remove_existing_all_models_file(...)
  c.remove_all_existing_output(...)
  c.backup_config_file(...)

which you can add to your main script, with caution! The different options may be useful if you want to delete some but not all previous output, e.g. to re-calculate weights but keep old orbit libraries. The API documentation has more information on the different options.

Logging output
===================

Logging is handled by the Python `logging <https://docs.python.org/3/library/logging.html>`_ module and by default uses your logging settings in the main script.

If you don't want to think about logging, you can activate the DYNAMITE standard logging settings by specifying ``reset_logging=True`` when reading the configuration file:

.. code-block:: python

  import dynamite as dyn
  c = dyn.config_reader.Configuration('config_file.yaml', reset_logging=True)

This will write logging messages of at least level ``INFO`` to the console and messages of at least level ``DEBUG`` to the log-file ``dynamite.log``. The levels, in increasing level of detail, are ``CRITICAL``, ``ERROR``, ``WARNING``, ``INFO``, ``DEBUG`` (currently, DYNAMITE does not use ``CRITICAL``).
If you (optionally) wish to control the verbosity of the logging output, do not use ``reset_logging=True`` but add the following lines near the top of the main script,

.. code-block:: python

  import logging
  dyn.config_reader.DynamiteLogging(
                        logfile='dynamite.log',
                        console_level=logging.INFO,
                        logfile_level=logging.DEBUG)

then you change the name of the log-file, and the level of logging output sent to the console and to the logfile. The values shown above are the defaults.

By default, the logging output is recorded in the file ``dynamite.log``, but you can also specify a different file name by using the parameter ``user_logfile`` (in ``Configuration``) or ``logfile`` (in ``DynamiteLogging``). These parameters can take, as values: a string indicating your desired name for the logfile (``.log`` will be appended to the string you provide), ``False`` (which will create a UTC-timestamped logfile ``dynamiteYYMMDD-HHMMSSuuuuuu.log``), or ``None`` (which will **not** create a logfile). This option can be useful especially if you are launching multiple DYNAMITE runs from the same directory, because otherwise the logging from all the runs will be written in the same ``dynamite.log`` file.
