# %load_ext autoreload
# %autoreload 2

import numpy as np
import dynamite_src as dyn
import matplotlib.pyplot as plt
import os
from astropy.io import ascii

# read configuration
fname = 'pj_model/config_legacy_example.yaml'
c = dyn.config_reader.ConfigurationReaderYaml(fname, silent=True)
parspace = dyn.parameter_space.ParameterSpace(c.system)

# read all existing models, or make empty all_models object if none
fname = c.settings.io_settings['output_directory']+'all_models.ecsv'
if os.path.isfile(fname):
    all_models = dyn.schwarzschild.AllModels(from_file=True,
                                             filename=fname,
                                             parspace=parspace,
                                             settings=c.settings)
else:
    all_models = dyn.schwarzschild.AllModels(from_file=False,
                                             parspace=parspace,
                                             settings=c.settings)

# run models
smi = dyn.shwarzschild_model_iterator.SchwarzschildModelIterator(
    system=c.system,
    all_models=all_models,
    settings=c.settings)

# save the all_models table
fname = c.settings.io_settings['output_directory']+'all_models.ecsv'
all_models.table.write(fname, format='ascii.ecsv')










# end
