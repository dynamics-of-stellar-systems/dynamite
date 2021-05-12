.. _api_configuration:

*************
Configuration
*************

Describe the configuration file and configuration object.

You can read the configuration file into a configuration object ``c`` as follows

.. code-block:: python

  import dynamite as dyn
  c = dyn.config_reader.Configuration('config_file.yaml') # read the configuration fie

``system_attributes``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

which can be accessed in the configuration object as ``c.system.distMPc`` etc.

``system_components``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The components are stored in the configuration object ``c`` as a list ``c.system.cmp_list``, with one entry for each of the components.

``system_parameters``
=====================

.. code-block:: python

  cmp_0 =  c.system.cmp_list[0]
  cmp_0.pars


``orblib_settings``
=====================

.. code-block:: python

  c.settings.orblib_settings

``weight_solver_settings``
==========================

.. code-block:: python

  c.settings.weight_solver_settings


``parameter_space_settings``
============================

.. code-block:: python

  c.settings.parameter_space_settings


``io_settings``
=====================

Hello!

``multiprocessing_settings``
============================

Hello!


``legacy_settings``
=====================

Hello!


Inheritance Diagram
===================

.. inheritance-diagram:: config_reader
