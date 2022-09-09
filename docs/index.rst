.. DYNAMITE documentation master file, created by
   sphinx-quickstart on Thu Sep 24 12:12:24 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DYNAMITE's documentation!
====================================

**DYNAMITE** (DYnamics, Age and Metallicity Indicators Tracing Evolution) is a tool for orbit-based dynamical modelling of stellar systems. It will soon be enhanced with stellar-population modelling tools.

How to cite
================

If you use DYNAMITE, please cite our `ASCL entry <http://www.ascl.net/code/v/2684>`_ using the following `BibTex citation <https://ui.adsabs.harvard.edu/abs/2020ascl.soft11007J/exportcitation>`_ and `Thater et al. 2022 <https://arxiv.org/abs/2205.04165>`_, using the following `BibTex citation <https://ui.adsabs.harvard.edu/abs/2022arXiv220504165T/exportcitation>`_.

===================
Orbit mirroring bug
===================

`Quenneville et al. 2022 <https://iopscience.iop.org/article/10.3847/1538-4357/ac3e68>`_ reported a bug in the orbit calculation of the original `van den Bosch et al. 2008 <https://academic.oup.com/mnras/article/385/2/647/1068433>`_ code that is used in DYNAMITE. We corrected the old mirroring in DYNAMITE, starting from version 3.0.0. From this version on, the default DYNAMITE run uses the correct mirrroring. We caution the user to not use the old orbit calculation routine.

In `Thater et al. 2022 <https://arxiv.org/abs/2205.04165>`_, our team provides a thorough quantification of how this bug has affected the results of dynamical analyses performed with previous versions of the code. Focusing on the typical scientific applications of the Schwarzschild triaxial code, in all our tests we find that differences are negligible with respect to the statistical and systematic uncertainties.

The bug occurred in the orbit calculation (legacy-fortran/orblib_f.f90), where a few velocity components needed a different sign (as reported in Table 1 of Quenneville et al. 2022). The corrected orbit calculation is found in (legacy-fortran/orblib_f_new_mirror.f90).

Getting Started
================

To get started with DYNAMITE,

1. Get the latest stable version from our `GitHub release page <https://github.com/dynamics-of-stellar-systems/dynamite/releases>`_. If you want the current version in development, you can also download this from our `GitHub page <https://github.com/dynamics-of-stellar-systems/dynamite>`_.
2. Install. The `installation page <https://dynamics.univie.ac.at/dynamite_docs/getting_started/installation.html>`_ has the full instructions. An overview is:

  a. Install Galahad: in the directory ``legacy_fortran/galahad-2.3/`` run the command ``./install_Galahad``
  b. Compile the Fortran programs: in the directory ``legacy_fortran/``  run the command ``make all``
  c. Install DYNAMITE Python package: in the main directory run the command ``python setup.py install``

3. Run orbit-based models! The following program will run a single model,

.. code-block:: python

   import dynamite as dyn

   c = dyn.config_reader.Configuration('my_config.yaml') # read configuration fie
   parset = c.parspace.get_parset()                      # extract a parameter set from configuration
   model = dyn.model.Model(config=c,parset=parset)       # make a model object
   model.setup_directories() # make directory tree
   model.get_orblib()        # make an orbit library
   model.get_weights()       # find orbital weights

The following pages give all the information needed to get started,

.. toctree::
  :maxdepth: 1

  getting_started/installation.rst
  getting_started/code_overview.rst
  getting_started/configuration.rst
  getting_started/getting_help.rst

The tutorials also show an example of running DYNAMITE from start to finish - this could also be a great place to start getting acquainted with the code. Further sections show other API-documentation for specific classes and methods, and other miscellaneous information.

Tutorials
=========

The following tutorials give detailed walkthroughs for using DYNAMITE.
Each page is an ipython notebook which you can either view in the browser, or download and interact with yourself.

.. toctree::
   :maxdepth: 1
   :caption: Tutorials:

   tutorial_notebooks/1_data_prep_for_gauss_hermites.ipynb
   tutorial_notebooks/2_quickstart.ipynb
   tutorial_notebooks/3_model_iterations_and_plots.ipynb
   tutorial_notebooks/4_BayesLOSVD.ipynb
   tutorial_notebooks/5_parameter_space.ipynb
   tutorial_notebooks/6_orbits_and_weights.ipynb


API Documentation
=================

These pages contain DYNAMITE's API documentation, i.e. the classes and method definitions and descriptions.
This may be useful for code developers, or anyone who intends to modify the code for their own personal use.
The API overview may be a useful starting point.

.. toctree::
   :maxdepth: 1
   :caption: API Documentation:

   api_docs/overview
   api_docs/configuration
   api_docs/data
   api_docs/mges
   api_docs/kinematics
   api_docs/physical_system
   api_docs/parameter_space
   api_docs/model
   api_docs/orblib
   api_docs/weight_solver
   api_docs/model_iterator
   api_docs/plotting

More Information
================

.. toctree::
   :maxdepth: 1
   :caption: More information:

   more_info/team.rst
   more_info/publications.rst
   more_info/changelog.rst
   more_info/getting_involved.rst
   more_info/license.rst
   more_info/making_the_docs.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
