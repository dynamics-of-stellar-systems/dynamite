.. _api_overview:

************
API overview
************

These pages contain DYNAMITE's API documentation, i.e. the classes and method definitions and descriptions. This may be useful for code developers, or anyone who intends to modify the code for their own personal use.

As a broad overview, a DYNAMITE run has the following steps. First, `configuration <configuration>`__ file is read. This will hold all of the user-specified settings. In particular, it will create an object to store information about the `physical_system <physical_system>`__ being modelled (e.g. the galaxy, and its constituent components), as well as any observational `data <data>`__, e.g. `mges <mges>`__ and `kinematics <kinematics>`__. A set of `parameters <parameter_space>`__ is required to define a `model <model>`__. For a given model, DYNAMITE calculates an `orbit library <orblib>`__, and then solves for the orbital `weights <weight_solver>`__ which best reproduce the observed kinematics.  This process is `iterated <model_iterator>`__, using algorithms to `vary the model parameters <parameter_space>`__ till an optimum parameter set is obtained. Throughout the process `plots <plotting>`__ are made to visualise the procedure, and additional plotting routines are available to use once you have found a satisfactory fit.

The following pages contain the API documentation for this procedure, and each page corresponds to a single file of DYNAMITE source code. The bottom of each page shows the inheritance diagram for that section of code.

.. toctree::
  :maxdepth: 1

  configuration
  data
  mges
  kinematics
  populations
  physical_system
  parameter_space
  model
  orblib
  weight_solver
  model_iterator
  plotting
  analysis
  coloring

Inheritance Diagram
===================

The inheritance diagram showing class relations between all DYNAMITE classes:

.. inheritance-diagram:: parameter_space model_iterator model populations weight_solvers orblib kinematics data physical_system config_reader plotter analysis
  :parts: 1
