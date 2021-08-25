#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import logging
import time

import numpy as np

# Set matplotlib backend to 'Agg' (compatible when X11 is not running
# e.g., on a cluster). Note that the backend can only be set BEFORE
# matplotlib is used or even submodules are imported!
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from astropy import table
import dynamite as dyn

def run_user_test(make_comp=False):

    logger = logging.getLogger()
    logger.info(f'Using DYNAMITE version: {dyn.__version__}')
    logger.info(f'Located at: {dyn.__path__}')
    # print to console anyway...
    print('Using DYNAMITE version:', dyn.__version__)
    print('Located at:', dyn.__path__)

    # read configuration
    if '__file__' in globals():
        file_dir = os.path.dirname(__file__)
        if file_dir:
            os.chdir(file_dir)
    fname = 'dif_dm_halos_config.yaml'
    c = dyn.config_reader.Configuration(fname, reset_logging=True)
    # delete previous output if available
    c.remove_all_existing_output(wipe_all=True, create_tree=True)

    # "run" the models
    t = time.perf_counter()
    smi = dyn.model_iterator.ModelIterator(c)

    delt = time.perf_counter()-t
    logger.info(f'Computation time: {delt} seconds = {delt/60} minutes')
    # print to console regardless of logging level
    print(f'Computation time: {delt} seconds = {delt/60} minutes')

    # print all model results
    c.all_models.table.pprint(max_lines=-1, max_width=-1)

    return

if __name__ == '__main__':
    run_user_test()

# end
