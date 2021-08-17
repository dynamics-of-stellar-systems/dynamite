.. DYNAMITE documentation master file, created by
   sphinx-quickstart on Thu Sep 24 12:12:24 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DYNAMITE's documentation!
====================================

**DYNAMITE** (DYnamics, Age and Metallicity Indicators Tracing Evolution) is a tool for orbit-based dynamical modelling of stellar systems. It will soon be enhanced with stellar-population modelling tools.

How to cite
================

If you use DYNAMITE, please cite our `ASCL entry <http://www.ascl.net/code/v/2684>`_ using the following `BibTex citation <https://ui.adsabs.harvard.edu/abs/2020ascl.soft11007J/exportcitation>`_

Getting Started
================

To get started with DYNAMITE,

1. Get the code from our `GitHub page <https://github.com/dynamics-of-stellar-systems/dynamite_release/releases>`_
2. Install. The `installation page <https://www.univie.ac.at/dynamics/dynamite_docs/installation.html>`_ has the full instructions. An overview is:

  a. Install Galahad: in the directory ``legacy_fortran/galahad-2.3/`` run the command ``./install_Galahad``
  b. Compile the Fortran programs: in the directory ``legacy_fortran/``  run the command ``make all``
  c. Install DYNAMITE Python package: in the main directory run the command ``python setup.py install``

3. Run orbit-based models! The following program will run a single model,

.. code-block:: python

   import dynamite as dyn

   c = dyn.config_reader.Configuration('my_config.yaml') # read configuration fie
   parset = c.parspace.get_parset() # extract a parameter set from configuration
   model = dyn.model.Model(
     system=c.system,
     settings=c.settings,
     parspace=c.parspace,
     parset=parset)          # make a model object
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

   tutorial_notebooks/data_prep_for_gauss_hermites.ipynb
   tutorial_notebooks/Quickstart.ipynb
   tutorial_notebooks/BayesLOSVD_and_DYNAMITE.ipynb

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
