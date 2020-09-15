%load_ext autoreload
%autoreload 2

import numpy as np
import dynamite_src as dyn
import matplotlib.pyplot as plt
import os
from astropy.io import ascii

'''
To run this example, the current directory should contain a directory
model_example with the following contents:
    model_example/NGC6278/
                    config_legacy_example.yaml
                    input_data/
                        aperture.dat
                        bins.dat
                        gauss_hermite_kins.ecsv
                        mge.ecsv
The fowllowing entry in the config file should be changed to the location of
your local triaxialschwarzschild directory:
    legacy_settings:
        directory: xxx
'''

# read configuration
fname = 'model_example/NGC6278/config_legacy_example.yaml'
c = dyn.config_reader.ConfigurationReaderYaml(fname, silent=True)

parspace = dyn.parameter_space.ParameterSpace(c.system)
print('Free parameters are:')
for p in parspace:
    if p.fixed is False:
        print(f'... {p.name}')

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

# for speed of testing, I'll set the following kw arguments in order to:
# (i) not exectute the models
# (ii) use a fake chi2, as defined by this fake_chi2_function:
def fake_chi2_function(parset):
    fake_chi2 = np.log10(parset['mass']) + parset['ml'] + np.log10(parset['f'])
    return fake_chi2
model_kwargs = {'execute_run':False,
                'use_fake_chi2':True,
                'fake_chi2_function':fake_chi2_function}

model_kwargs = {'execute_run':True,
                'use_fake_chi2':False}

print('Instantiate executor object')
print('    executor type: ', c.settings.executable_settings['type'])
kw_executor = {'system':c.system,
               'legacy_directory':c.settings.legacy_settings['directory']}
executor = getattr(dyn.executor,
                   c.settings.executable_settings['type'])(**kw_executor)

# "run" the models
smi = dyn.shwarzschild_model_iterator.SchwarzschildModelIterator(
    system=c.system,
    all_models=all_models,
    settings=c.settings,
    model_kwargs=model_kwargs,
    executor=executor)

all_models.table['which_iter']

# plot the models: f vs ml at each iteration:
for iter in np.unique(all_models.table['which_iter']):
    table = all_models.table
    table = table[table['which_iter']==iter]
    plt.scatter(table['f'],
                table['ml'],
                c=table['chi2'],
                cmap=plt.cm.viridis_r,
                s=200)
    # plt.colorbar()
    plt.gca().set_xlim(1e-2, 1e3)
    plt.gca().set_ylim(0, 7)
    plt.gca().set_xscale('log')
    plt.show()

# plot the models: f vs ml altogether
plt.scatter(all_models.table['f'],
            all_models.table['ml'],
            c=all_models.table['chi2']
            cmap=plt.cm.viridis_r,
            s=200)
plt.gca().set_xscale('log')
plt.show()

# save the all_models table
fname = c.settings.io_settings['output_directory']+'all_models.ecsv'
all_models.table.write(fname, format='ascii.ecsv')












# end
