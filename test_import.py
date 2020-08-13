%load_ext autoreload
%autoreload 2

import numpy as np
import dynamite_src as dyn

# fname = './datafiles/config_example.yaml'
fname = './datafiles/config_legacy_example.yaml'
c = dyn.config_reader.ConfigurationReaderYaml(fname)

# extract parameter space
parspace = dyn.parameter_space.ParameterSpace(c.system)

len(parspace)
parspace[0]
parspace[0]

# extract parameter space
all_models = dyn.schwarzschild.AllModels(from_file=False,
                                         parspace=parspace,
                                         config=c.config)
all_models.table

all_models.convert_legacy_chi2_file(
    legacy_filename='outputs/legacy/NGC6278/griddata/_chi2.cat',
)
all_models = dyn.schwarzschild.AllModels(config=c.config)

all_models.table

# take the first row of completed models for an example parameter set
parset0 = all_models.table[0]
parset0 = parset0[parspace.par_names]

# create a model object
mod = dyn.schwarzschild.LegacySchwarzschildModel(
    system=c.system,
    config=c.config,
    parspace=parspace,
    parset=parset0)
mod.chi2
mod.kinchi2















# end
