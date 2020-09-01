%load_ext autoreload
%autoreload 2

import numpy as np
import dynamite_src as dyn

# fname = './datafiles/config_example.yaml'
fname = './datafiles/config_legacy_example.yaml'
c = dyn.config_reader.ConfigurationReaderYaml(fname)

# extract parameter space
parspace = dyn.parameter_space.ParameterSpace(c.system)

# parspace is a list of parameter objects
print(len(parspace))

parspace[0]

# parspace is a list of parameter objects
parspace[0]

c.settings.orblib_settings


# extract parameter space
all_models = dyn.schwarzschild.AllModels(from_file=False,
                                         parspace=parspace,
                                         settings=c.settings)

print(all_models.table)

print('parameter_space start...')

# Instantiate GridSearch object
g = dyn.parameter_space.GridSearch(parspace)

# generate first model list: grid walk "iteration 0"
print('grid walk "iteration 0"')
g.generate(current_models = all_models)
print(all_models.table)
print(g.status)

# generate second model list: grid walk "iteration 1"
print('grid walk "iteration 1"')
g.generate(current_models=all_models)
print(all_models.table)
print(g.status)

# generate second model list: grid walk "iteration 2"
print('grid walk "iteration 2"')
g.generate(current_models=all_models)
print(all_models.table)
print(g.status)

all_models.convert_legacy_chi2_file(
    legacy_filename='outputs/legacy/NGC6278/griddata/_chi2.cat',
)
all_models = dyn.schwarzschild.AllModels(settings=c.settings)

print(all_models.table)

# generate model list based on current models: grid walk "iteration n"
# NOTE: the read-in data only has 7 parameters and does not match the parspace
#       read from the config file...
g.generate(current_models = all_models)
print(all_models.table)



# take the first row of completed models for an example parameter set
parset0 = all_models.table[0]

parset0

parspace.par_names

parset0 = parset0[parspace.par_names]
parset0
# create a model object
mod = dyn.schwarzschild.LegacySchwarzschildModel(
    system=c.system,
    settings=c.settings,
    parspace=parspace,
    parset=parset0)

c.system.cmp_list[0].parameters

mod.chi2
mod.kinchi2




# end
