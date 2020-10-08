#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# first, make sure the paths are set up
# we assume that this script is located and run in the folder dynamite/tests

import os
import sys

os.chdir('..')
this_dir = os.getcwd()
#print(f'working directory: {this_dir}')
if not this_dir in sys.path:
    sys.path.append(this_dir)

import numpy as np
import dynamite as dyn
import matplotlib.pyplot as plt
import time

def run_user_test():
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
    plt.xticks(range(0,max(c.all_models.table['which_iter'])+1))
    plt.ylabel('mass')
    plt.show()
    
    # save the all_models table
    c.all_models.save()

    return c.all_models.table

if __name__ == '__main__':
    run_user_test()

# end
