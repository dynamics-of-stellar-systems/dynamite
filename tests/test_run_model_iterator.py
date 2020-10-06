%load_ext autoreload
%autoreload 2

# import sys
# sys.path.append("..")

import numpy as np
import dynamite as dyn
import matplotlib.pyplot as plt

# read configuration
fname = 'tests/NGC6278/config.yaml'
c = dyn.config_reader.ConfigurationReaderYaml(fname, silent=True)

# for speed of testing, we can set the following kw arguments in order to:
# (i) not exectute the models
# (ii) use a dummy chi2, as defined by this dummy_chi2_function:
# do_dummy_run = True
# def dummy_chi2_function(parset):
#     chi2 = parset['ml'] + np.log10(parset['f'])
#     return chi2
# ... or, to actually run models, use these values (which are the defaults)
do_dummy_run = False
dummy_chi2_function = None

# "run" the models
smi = dyn.model_iterator.ModelIterator(
    system=c.system,
    all_models=c.all_models,
    settings=c.settings,
    executor=c.executor,
    do_dummy_run=do_dummy_run,
    dummy_chi2_function=dummy_chi2_function)

print(len(c.all_models.table))

# plot the models: f vs ml at each iteration:
for iter in np.unique(c.all_models.table['which_iter']):
    table = c.all_models.table
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
plt.scatter(c.all_models.table['f'],
            c.all_models.table['ml'],
            c=c.all_models.table['chi2'],
            cmap=plt.cm.viridis_r,
            s=200)
plt.colorbar()
plt.gca().set_title(f'all iterations')
plt.gca().set_xlim(1e0, 1e2)
plt.gca().set_ylim(3, 7)
plt.gca().set_xscale('log')
plt.show()

# save the all_models table
c.all_models.save()

# end
