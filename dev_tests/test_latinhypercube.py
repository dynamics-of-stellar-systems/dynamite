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

import dynamite as dyn

print('Using DYNAMITE version:', dyn.__version__)
print('Located at:', dyn.__path__)

# read configuration
fname = 'test_latinhypercube.yaml'
c = dyn.config_reader.Configuration(fname)
c.remove_all_existing_output(wipe_all=True, create_tree=True)

# run the models
t = time.perf_counter()
smi = dyn.model_iterator.ModelIterator(c)
delt = time.perf_counter()-t
print(f'Computation time: {delt} seconds = {delt/60} minutes')
