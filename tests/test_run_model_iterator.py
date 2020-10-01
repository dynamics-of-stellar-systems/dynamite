%load_ext autoreload
%autoreload 2

# import sys
# sys.path.append("..")

import numpy as np
import dynamite as dyn
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
fname = 'tests/NGC6278/config.yaml'
c = dyn.config_reader.ConfigurationReaderYaml(fname, silent=True)
parspace = dyn.parameter_space.ParameterSpace(c.system)
all_models = dyn.model.AllModels(parspace=parspace, settings=c.settings)

len(all_models.table)

# for speed of testing, I'll set the following kw arguments in order to:
# (i) not exectute the models
# (ii) use a dummy chi2, as defined by this dummy_chi2_function:
do_dummy_run = True
def dummy_chi2_function(parset):
    chi2 = parset['ml'] + np.log10(parset['f'])
    return chi2

# Or, to actually run models, do not pass the above keywords:
do_dummy_run = False
dummy_chi2_function = None

print('Hi from Prash - :D !!!!!')
print('Instantiate executor object')
print('    executor type: ', c.settings.executor_settings['type'])
kw_executor = {'system':c.system,
               'legacy_directory':c.settings.legacy_settings['directory'],
               'executor_settings':c.settings.executor_settings}
executor_type = c.settings.executor_settings['type']
executor = getattr(dyn.executor, executor_type)(**kw_executor)

# "run" the models
smi = dyn.model_iterator.ModelIterator(
    system=c.system,
    all_models=all_models,
    settings=c.settings,
    executor=executor,
    do_dummy_run=do_dummy_run,
    dummy_chi2_function=dummy_chi2_function)




len(all_models.table)

# plot the models: f vs ml at each iteration:
for iter in np.unique(all_models.table['which_iter']):
    table = all_models.table
    table = table[table['which_iter']==iter]
    plt.scatter(table['f'],
                table['ml'],
                c=table['chi2'],
                cmap=plt.cm.viridis_r,
                s=200)
    plt.colorbar()
    plt.gca().set_title(f'iteration {iter}')
    plt.gca().set_xlim(1e0, 1e2)
    plt.gca().set_ylim(3, 7)
    plt.gca().set_xscale('log')
    plt.show()

# plot the models: f vs ml altogether
plt.scatter(all_models.table['f'],
            all_models.table['ml'],
            c=all_models.table['chi2'],
            cmap=plt.cm.viridis_r,
            s=200)
plt.colorbar()
plt.gca().set_title(f'all iterations')
plt.gca().set_xlim(1e0, 1e2)
plt.gca().set_ylim(3, 7)
plt.gca().set_xscale('log')
plt.show()

# save the all_models table
all_models.save()














# end
