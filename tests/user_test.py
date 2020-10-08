#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# first, make sure the paths are set up
# we assume that this script is located and run in the folder dynamite/tests

import os
import sys
import shutil

os.chdir('..')
this_dir = os.getcwd()
if not this_dir in sys.path:
    sys.path.append(this_dir)

import time
import dynamite as dyn
import matplotlib.pyplot as plt
from astropy import table

def run_user_test():

    # delete previous output if available
    models_folder = this_dir+'/'+'tests/NGC6278/models'
    models_file = this_dir+'/'+'tests/NGC6278/all_models.ecsv'
    shutil.rmtree(models_folder, ignore_errors=True)
    if os.path.isfile(models_file):
        os.remove(models_file)

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
    plt.gca().set_title('all iterations')
    plt.gca().set_yscale('log')
    plt.xlabel('iteration')
    plt.xticks(range(0,max(c.all_models.table['which_iter'])+1))
    plt.ylabel('mass')
    plt.show()

    # save the all_models table
    c.all_models.save()

    # compare to chi2 in chi2_compare.dat
    chi2_compare = table.Table.read('tests/chi2_compare.dat', format='ascii')
    plt.vlines(chi2_compare['model_id'],
               chi2_compare['chi2_min'],
               chi2_compare['chi2_max'])
    plt.plot([i for i in range(len(c.all_models.table))],
             c.all_models.table['chi2'],
             'rx')
    plt.gca().set_title('calculated chi2 vs should-be range')
    plt.xlabel('model_id')
    plt.xticks(range(len(c.all_models.table)))
    plt.ylabel('chi2')
    plt.show()
    print(f'chi2 statistics for comparison:\n{chi2_compare}')

    return c.all_models.table

if __name__ == '__main__':
    run_user_test()

# end
