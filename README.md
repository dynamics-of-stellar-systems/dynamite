# DYNAMITE

DYnamics, Age and Metallicity Indicators Tracing Evolution.

A package for Schwarzschild- and stellar-population modelling of stellar systems.

## Installation

See [the documentation](https://dynamics.univie.ac.at/dynamite_docs/getting_started/installation.html) for full instructions. An overview is:
1. Install the Galahad optimization library
   - ``./install_galahad`` in the directory ``legacy_fortran/galahad-2.3/``
2. Compile the Fortran programs
   - ``make all`` in the directory ``legacy_fortran/``
3. Install DYNAMITE Python package
   - ``python -m pip install .`` in the root directory (note the dot at the end of the line)

## Usage

```python
import dynamite as dyn

c = dyn.config_reader.Configuration('my_config.yaml') # read configuration
_ = dyn.model_iterator.ModelIterator(c) # the model iterator populates the all_models.ecsv table
```

## License

[MIT](https://choosealicense.com/licenses/mit/)
