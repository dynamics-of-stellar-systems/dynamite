# import sys
# sys.path.append("..")

# assuming this script is in dynamite/tests, change to folder dynamite
import os
os.chdir('..')

import numpy as np
import dynamite as dyn
import matplotlib.pyplot as plt
import time

# read configuration
fname = 'tests/NGC6278/user_test_config.yaml'
c = dyn.config_reader.ConfigurationReaderYaml(fname, silent=True)

# "run" the models
t = time.perf_counter()
smi = dyn.model_iterator.ModelIterator(
    system=c.system,
    all_models=c.all_models,
    settings=c.settings,
    executor=c.executor)
delt = time.perf_counter()-t
print(f'Computation time: {delt} seconds = {delt/60} minutes')

# print all model results
c.all_models.table.pprint_all()

# plot the models
plt.scatter(c.all_models.table['which_iter'],
            c.all_models.table['mass'],
            c=c.all_models.table['chi2'],
            cmap=plt.cm.viridis_r,
            s=200)
cb=plt.colorbar()
cb.set_label('chi2', y=1.1, labelpad=-40, rotation=0)
plt.gca().set_title(f'all iterations')
plt.gca().set_yscale('log')
plt.xlabel('iteration')
plt.ylabel('mass')
#plt.gca().set_xlim(1e0, 1e2)
#plt.gca().set_ylim(3, 7)
#plt.gca().set_xscale('log')
plt.show()


# # plot the models: f vs ml at each iteration:
# for iter in np.unique(c.all_models.table['which_iter']):
#     table = c.all_models.table
#     table = table[table['which_iter']==iter]
#     plt.scatter(table['f'],
#                 table['ml'],
#                 c=table['chi2'],
#                 cmap=plt.cm.viridis_r,
#                 s=200)
#     plt.colorbar()
#     plt.gca().set_title(f'iteration {iter}')
#     plt.gca().set_xlim(1e0, 1e2)
#     plt.gca().set_ylim(3, 7)
#     plt.gca().set_xscale('log')
#     plt.show()

# # plot the models: f vs ml altogether
# plt.scatter(c.all_models.table['f'],
#             c.all_models.table['ml'],
#             c=c.all_models.table['chi2'],
#             cmap=plt.cm.viridis_r,
#             s=200)
# plt.colorbar()
# plt.gca().set_title(f'all iterations')
# plt.gca().set_xlim(1e0, 1e2)
# plt.gca().set_ylim(3, 7)
# plt.gca().set_xscale('log')
# plt.show()

# save the all_models table
c.all_models.save()

# end
