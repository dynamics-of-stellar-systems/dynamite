#%load_ext autoreload
#%autoreload 2

import numpy as np
import dynamite_src as dyn
#from astropy.table import vstack

# fname = './datafiles/config_example.yaml'

fname = './model_example/NGC6278/input_data/config_legacy_example.yaml'
fname = './model_example/NGC6278/config_legacy_example.yaml'
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
parspace_settings = c.settings.parameter_space_settings
#g = dyn.parameter_space.GridWalk(parspace, parspace_settings=parspace_settings)
g = dyn.parameter_space.LegacyGridSearch(parspace, parspace_settings=parspace_settings)

# generate first model list: "iteration 0"
print('grid search / grid walk "iteration 0"')
g.generate(current_models = all_models)
print(all_models.table)
print(g.status)

# generate second model list: "iteration 1"
print('grid search / grid walk "iteration 1"')
g.generate(current_models=all_models)
print(all_models.table)
print(g.status)

# generate second model list: "iteration 2"
print('grid search / grid walk "iteration 2"')
g.generate(current_models=all_models)
print(all_models.table)
print(g.status)

print(all_models.table['f', 'ml', 'chi2', 'kinchi2', 'time_modified', 'which_iter'])

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

c.system.cmp_list[0].parameters


#run a first dynamite model
print('-----------------------------')
print('Run a dynamical model')
print('-----------------------------')


#directory of the models, need to double check where this is done
#c.config.output_settings['directory']='model_example/'

# create a model object
mod = dyn.schwarzschild.LegacySchwarzschildModel(
    system=c.system,
    settings=c.settings,
    parspace=parspace,
    parset=parset0)


#print(mod.directory)
#print(mod.directory[:-7])
##mod.chi2
##mod.kinchi2

#for sabines tests
#print('-----------------------------')
#print('-----------------------------')
#print(c.system.cmp_list[2].mge.data)
#print(c.system.cmp_list[2].mge.data[0][2])
#print(c.system.distMPc)
#print(parspace.par_names)
#print(parspace[6].value)
#print(c.config.orblib_settings['nE'])
#print(c.system.cmp_list[2].kinematic_data[0].PSF['weight'])

#.__dict__
#.--dict--.keys()


#os.chdir("/Users/sabine/Work/Dynamics/DYNAMITE/Version1/Fortran/model_example/NGC6278/bh6.00dc1.00f2.00q0.54p0.970u0.9999/infil")
#dir="/Users/sabine/Work/Dynamics/DYNAMITE/Version1/Fortran/model_example/NGC6278/bh6.00dc1.00f2.00q0.54p0.970u0.9999/infil/"
#for file in glob.glob(dir+"*.dat"):
#    print(file)


# end
