.. _api_configuration:

*************
Configuration
*************

All configuration settings are stored in the configuration object ``c``, which is structured very similarly to structure of the YAML configuration file,

.. code-block:: python

  import dynamite as dyn

  c = dyn.config_reader.Configuration('config_file.yaml') # read the configuration fie

  ######################################################
  # config. options relating to the physical system
  ######################################################

  # system_attributes accessed via:
  c.system.distMPc # etc

  # system_components stores in the list:
  c.system.cmp_list
  c.system.cmp_list[0]    # the 0'th component,
  c.system.cmp_list[1]    # the 1'st component, etc...

  # component parameters are stored in the list:
  c.system.cmp_list[0].parameters

  # system_parameters are stored in:
  c.system.parameters

  # all parameters (component and system) are stored together in the list:
  c.parspace

  ######################################################
  # config. options relating to other settings
  ######################################################

   # orblib_settings stored in the dictionary:
   c.settings.orblib_settings

   # weight_solver_settings stored in the dictionary:
   c.settings.weight_solver_settings

   # multiprocessing_settings stored in the dictionary:
   c.settings.multiprocessing_settings

   # io_settings stored in the dictionary:
   c.settings.io_settings

   # legacy_settings stored in the dictionary:
   c.settings.legacy_settings

API
===================

 .. automodule:: config_reader
  :members:

Inheritance Diagram
===================

.. inheritance-diagram:: config_reader
