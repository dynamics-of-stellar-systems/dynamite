#!/bin/env python3
#SBATH --qos=p71474_0096
#SBATCH --job-name=test_dynamite
#SBATCH -N 1
#SBATCH --mem-per-cpu=10
#SBATCH --output=dynamite_output.log
#SBATCH --error=dynamite_error.log

# new in config file - multiprocessing_settings: ncpus
# set this either to an integer or to 'all_available'

# to submit this locally, run ``python test_slurm.py``
# to submit this on cluster with slurm run ``sbatch test_slurm.py``

import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os
import shutil

import dynamite as dyn

print('Using DYNAMITE version:', dyn.__version__)
print('Located at:', dyn.__path__)

# read configuration
fname = 'test_validate_pqu_config.yaml'
c = dyn.config_reader.Configuration(fname)
io_settings = c.settings.io_settings
outdir = io_settings['output_directory']
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# stars = c.system.get_component_from_name('stars')
#
# p = np.linspace(0.01, 0.99, 30)
# q = np.linspace(0.01, 0.99, 31)
# u = np.linspace(0.01, 0.99, 32)
# is_valid = np.zeros((30, 31, 32), dtype=bool)
# for i, p0 in enumerate(p):
#     for j, q0 in enumerate(q):
#         for k, u0 in enumerate(u):
#             pqu = {'p':p0, 'q':q0, 'u':u0}
#             is_valid[i,j,k] = stars.validate_parset(pqu)
# pp, qq, uu = np.meshgrid(p, q, u, indexing='ij')
# p_val = pp[is_valid]
# q_val = qq[is_valid]
# u_val = uu[is_valid]
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter3D(p_val, q_val, u_val, s=1, c=u_val)
# ax.set_xlabel('p')
# ax.set_ylabel('q')
# ax.set_zlabel('u')
# ax.set_title('Valid (p,q,u) combinations')
# fig.tight_layout()
# plt.show()

# run the models
t = time.perf_counter()
smi = dyn.model_iterator.ModelIterator(
    system=c.system,
    all_models=c.all_models,
    settings=c.settings,
    ncpus=c.settings.multiprocessing_settings['ncpus'])
delt = time.perf_counter()-t
print(f'Computation time: {delt} seconds = {delt/60} minutes')
