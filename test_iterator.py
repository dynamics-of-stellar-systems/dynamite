# %load_ext autoreload
# %autoreload 2

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
    all_models = dyn.model.AllModels(from_file=True,
                                     filename=fname,
                                     parspace=parspace,
                                     settings=c.settings)
else:
    all_models = dyn.model.AllModels(from_file=False,
                                     parspace=parspace,
                                     settings=c.settings)

# for speed of testing, I'll set the following kw arguments in order to:
# (i) not exectute the models
# (ii) use a dummy chi2, as defined by this dummy_chi2_function:
def dummy_chi2_function(parset):
    chi2 = parset['ml'] + np.log10(parset['f'])
    return chi2
model_kwargs = {'dummy_run':True,
                'dummy_chi2_function':dummy_chi2_function}
# Or, to actually run models, do not pass the above keywords:
#model_kwargs = {}

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
    model_kwargs=model_kwargs)

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
fname = c.settings.io_settings['output_directory']+'all_models.ecsv'
all_models.table.write(fname, format='ascii.ecsv', overwrite=True)
print(all_models.table['chi2', 'kinchi2'])











# end
