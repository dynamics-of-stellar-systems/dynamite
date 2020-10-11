.. DYNAMITE documentation master file, created by
   sphinx-quickstart on Thu Sep 24 12:12:24 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DYNAMITE's documentation!
====================================

**DYNAMITE** (DYnamics, Age and Metallicity Indicators Tracing Evolution) is a tool for stellar population and dynamical modelling.

How to cite
================

This is how to cite our code!

Getting Started
================

To get started with DYNAMITE,

1. Get the code from our `GitHub page <https://github.com/dynamics-of-stellar-systems/dynamite>`_
2. Install. Detailed instructions can be found `here <installation.rst>`_. An overview of the steps involved is:

   i. install Galahad
   ii. Make Fortran executables in the legacy_fortran directory
   iii. install Python dependencies: astropy, pyyaml, numpy, matplotlib
   iv. add the directory containing dynamite to your PYTHONPATH variable

3. Here is an example of how you can run a Schwarzschild model in DYNAMITE

.. code-block:: python

   import dynamite as dyn

   # read the configuration
   c = dyn.config_reader.Configuration('config_file.yaml')
   # extract a parameter set
   parset = c.parspace.get_parset()
   # make and run a Schwarzschild model
   model = dyn.model.LegacySchwarzschildModel(
     system=c.system,
     settings=c.settings,
     parspace=c.parspace,
     executor=c.executor,
     parset=parset)
   model.setup_directories()
   model.get_orblib()
   model.get_weights()

.. toctree::
  :maxdepth: 1
  :caption: Getting Started:

  installation.rst
  getting_help.rst

Tutorials
=========

The following tutorials give detailed walkthroughs for using DYNAMITE.
Each page is an ipython notebook which you can either view in the browser, or - preferably! - download and interact with yourself.

.. toctree::
   :maxdepth: 1
   :caption: Tutorials:

   tutorial_notebooks/running_a_model.ipynb
   tutorial_notebooks/running_a_grid_of_models.ipynb
   tutorial_notebooks/exploring_model_output.ipynb
   tutorial_notebooks/parameter_space.ipynb

Documentation
=============

The following pages describe the classes used in DYNAMITE.
Start with the `overview <classes/overview>`_ page for a description of all the
classes and how they interact. The subsequent pages then discuss the individual
classes in more detail.

.. toctree::
   :maxdepth: 3
   :caption: Documentation:

   classes/overview
   classes/configuration
   classes/physical_system
   classes/data
   classes/model
   classes/model_iterator
   classes/parameter_space
   classes/executor
   classes/plotting

More Information
================

.. toctree::
   :maxdepth: 2
   :caption: More information:

   about.rst
   publications.rst
   team.rst
   getting_involved.rst
   license.rst
   making_the_docs.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
