# DYNAMITE

DYnamics, Age and Metallicity Indicators Tracing Evolution.

A package for Schwarzschild- and stellar-population modelling of stellar systems.

## Installation

See [the documentation](https://dynamics.univie.ac.at/dynamite_docs/getting_started/installation.html) for full instructions.

There are two steps.

1. Compile and install the Galahad library and the DYNAMITE Fortran executables.
   - If you are using a modern Linux distribution with gfortran, then follow the instructions in ``legacy_fortran/README.linux``.
   - If you are using Mac OS X (or Linux with a nonstandard compiler) then
     1. Install the Galahad optimization library
       - ``./install_galahad`` in the directory ``legacy_fortran/galahad-2.3/``
     2. Compile the Fortran programs
       - ``make all`` in the directory ``legacy_fortran/``
2. Install the DYNAMITE Python package
   - ``python setup.py install`` in the root directory

## Usage

```python
import dynamite as dyn

c = dyn.config_reader.Configuration('my_config.yaml') # read configuration
parset = c.parspace.get_parset() # extract a parameter set from configuration
model = dyn.model.Model(
  system=c.system,
  settings=c.settings,
  parspace=c.parspace,
  parset=parset)          # make a Schwarzschild model
model.setup_directories() # make directory tree
model.get_orblib()        # make an orbit library
model.get_weights()       # find orbital weights
```

## License

[MIT](https://choosealicense.com/licenses/mit/)
