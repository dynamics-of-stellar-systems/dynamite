# DYNAMITE

DYnamics, Age and Metallicity Indicators Tracing Evolution.

A package for Schwarzschild- and stellar-population modelling of stellar systems.

## Installation

See [the documentation](https://www.univie.ac.at/dynamics/dynamite_docs/installation.html) for full instructions. An overview is:
1. Install the Galahad optimization library
   - ``./install_galahad`` in the directory ``legacy_fortran/galahad-2.3/``
2. Compile the Fortran programs
   - ``make all`` in the directory ``legacy_fortran/``
3. Install DYNAMITE Python package
   - ``python setup.py install`` in the root directory

## Usage

```python
import dynamite as dyn

c = dyn.config_reader.Configuration('my_config.yaml') # read configuration
parset = c.parspace.get_parset() # extract a parameter set from configuration
model = dyn.model.LegacySchwarzschildModel(
  system=c.system,
  settings=c.settings,
  parspace=c.parspace,
  executor=c.executor,
  parset=parset)          # make a Schwarzschild model
model.setup_directories() # make directory tree
model.get_orblib()        # make an orbit library
model.get_weights()       # find orbital weights
```

## License

[MIT](https://choosealicense.com/licenses/mit/)
