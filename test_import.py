#%load_ext autoreload
#%autoreload 2

import numpy as np
import dynamite_src as dyn
#from astropy.table import vstack

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
    config=c.config,
    parspace=parspace,
    parset=parset0)


#print(mod.directory)
#print(mod.directory[:-7])
##mod.chi2
##mod.kinchi2

print('-----------------------------')
print('-----------------------------')
print(c.system.cmp_list[2].mge.data)
#print(c.system.cmp_list[2].mge.data[0][2])
print(c.system.distMPc)
print(parspace.par_names)
print(parspace[6].value)
print(c.config.orblib_settings['nE'])
print(c.system.cmp_list[2].kinematic_data[0].PSF['weight'])

#.__dict__
#.--dict--.keys()


#os.chdir("/Users/sabine/Work/Dynamics/DYNAMITE/Version1/Fortran/model_example/NGC6278/bh6.00dc1.00f2.00q0.54p0.970u0.9999/infil")
#dir="/Users/sabine/Work/Dynamics/DYNAMITE/Version1/Fortran/model_example/NGC6278/bh6.00dc1.00f2.00q0.54p0.970u0.9999/infil/"
#for file in glob.glob(dir+"*.dat"):
#    print(file)


#psf=np.loadtxt(path+'datafiles/kinpsffile.dat',skiprows=2)

#path='/Users/sabine/Work/Dynamics/DYNAMITE/Version1/Fortran/model_example/NGC6278/bh6.00dc1.00f2.00q0.54p0.970u0.9999/infil/'
#len_mge=len(c.system.cmp_list[2].mge.data)
#np.savetxt(path+'parameters.in',c.system.cmp_list[2].mge.data,header=str(len_mge),comments='')


# end
