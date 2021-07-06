#!/bin/env python3
#SBATCH --qos=p71474_0096
#SBATCH --job-name=test_dynamite
#SBATCH -N 1
#SBATCH --mem-per-cpu=50
#SBATCH --output="dyn_%j.out"
#SBATCH --error="dyn_%j.err"

# new in config file - multiprocessing_settings: ncpus
# set this either to an integer or to 'all_available'

# to submit this locally, run ``python test_slurm.py``
# to submit this on cluster with slurm run ``sbatch test_slurm.py``

import time
import os
import logging

# Set matplotlib backend to 'Agg' (compatible when X11 is not running
# e.g., on a cluster). Note that the backend can only be set BEFORE
# matplotlib is used or even submodules are imported!
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import dynamite as dyn

def run_user_test():

    print('Using DYNAMITE version:', dyn.__version__)
    print('Located at:', dyn.__path__)

    # read configuration
    fname = 'test_slurm_config.yaml'
    c = dyn.config_reader.Configuration(fname, silent=True)

    # delete previous output if available
    c.remove_existing_orblibs()
    c.remove_existing_all_models_file()

    plotdir = c.settings.io_settings['plot_directory']
    plotfile = plotdir + 'slurm_model_timings.png'
    if os.path.isfile(plotfile):
        os.remove(plotfile)

    # "run" the models
    t = time.perf_counter()
    smi = dyn.model_iterator.ModelIterator(
        system=c.system,
        all_models=c.all_models,
        settings=c.settings,
        ncpus=c.settings.multiprocessing_settings['ncpus'])
    delt = time.perf_counter()-t
    print(f'Computation time: {delt} seconds = {delt/60} minutes')

    # print all model results
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(6,4))
    ax[0].plot(c.all_models.table['which_iter'], '.')
    ax[1].plot(c.all_models.table['time_modified'], '.')
    ax[1].set_xlabel('model ID')
    ax[0].set_ylabel('Iteration')
    ax[1].set_ylabel('time modified')
    ax[0].set_ylabel('Iteration')
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    fig.savefig(plotfile, dpi=300)
    plt.close()

    # print all model results
    print(f'Look at {plotfile}')
    c.all_models.table.pprint(max_lines=-1, max_width=-1)

    # select and print all models with (kin)chi2 values within a threshold
    # of delta
    best_models = c.all_models.get_mods_within_chi2_thresh(delta=300000)
    best_models.pprint(max_lines=-1, max_width=-1)

if __name__ == '__main__':

    # Example for directly setting up DynamiteLogging
    dyn.config_reader.DynamiteLogging(logfile='test_slurm.log',
                                      console_level=logging.ERROR,
                                      logfile_level=logging.DEBUG)

    run_user_test()

# end
