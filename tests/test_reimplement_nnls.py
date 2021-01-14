#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# first, make sure the paths are set up
# we assume that this script is located and run in the folder dynamite/tests

import os
#import sys
import shutil

# if os.getcwd().rpartition('/')[0] not in sys.path:
#     sys.path.append(os.getcwd().rpartition('/')[0])

import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import table
import dynamite as dyn

def run_user_test(stat_mode=False):

    print('Using DYNAMITE version:', dyn.__version__)
    print('Located at:', dyn.__path__)

    # read configuration
    fname = 'test_reimplement_nnls.yaml'
    c = dyn.config_reader.Configuration(fname, silent=True)

    io_settings = c.settings.io_settings
    outdir = io_settings['output_directory']
    # delete previous output if available
    models_folder = outdir + 'models/'
    models_file = outdir + io_settings['all_models_file']
    shutil.rmtree(models_folder, ignore_errors=True)
    if os.path.isfile(models_file):
        os.remove(models_file)
    plotdir = outdir + 'plots/'
    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)
    # plotfile_ml = plotdir + 'ml_vs_iter_chi2.png'
    # if os.path.isfile(plotfile_ml):
    #     os.remove(plotfile_ml)

    # re-read configuration now that old output has been deleted
    c = dyn.config_reader.Configuration(fname, silent=True)

    # "run" the models
    t = time.perf_counter()
    smi = dyn.model_iterator.ModelIterator(
        system=c.system,
        all_models=c.all_models,
        settings=c.settings,
        ncpus=c.settings.multiprocessing_settings['ncpus'])
    delt = time.perf_counter()-t
    print(f'Computation time: {delt} seconds = {delt/60} minutes')

    c.all_models.table.pprint(max_lines=-1, max_width=-1)

    return

if __name__ == '__main__':
    run_user_test()

# end
