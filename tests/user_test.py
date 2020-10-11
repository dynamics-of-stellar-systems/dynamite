#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# first, make sure the paths are set up
# we assume that this script is located and run in the folder dynamite/tests

import os
import sys
import shutil

if os.getcwd().rpartition('/')[0] not in sys.path:
    sys.path.append(os.getcwd().rpartition('/')[0])

import time
import matplotlib.pyplot as plt
from astropy import table
import dynamite as dyn

def run_user_test(stat_mode=False):

    if stat_mode==False:
        # set working directory to dynamite
        os.chdir(os.getcwd().rpartition('/')[0])

    # delete previous output if available
    models_folder = 'tests/NGC6278/models'
    models_file = 'tests/NGC6278/all_models.ecsv'
    shutil.rmtree(models_folder, ignore_errors=True)
    if os.path.isfile(models_file):
        os.remove(models_file)

    # read configuration
    fname = 'tests/NGC6278/user_test_config.yaml'
    c = dyn.config_reader.Configuration(fname, silent=True)
    stat_file = "tests/chi2_compare_" \
                f"{c.settings.orblib_settings['nE']}" \
                f"{c.settings.orblib_settings['nI2']}" \
                f"{c.settings.orblib_settings['nI3']}"
    if stat_mode==False:
        stat_file += "_10.dat"

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
#    c.all_models.table.pprint_all()
    c.all_models.table.pprint(max_lines=-1, max_width=-1)

    if stat_mode==False:
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
        chi2_compare = table.Table.read(stat_file, format='ascii')
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
        print('chi2 statistics for comparison:\n')
        print(chi2_compare.pprint(max_lines=-1, max_width=-1))

    return c.all_models.table, \
        stat_file, \
        c.settings.parameter_space_settings['stopping_criteria']['n_max_mods']

if __name__ == '__main__':
    run_user_test()

# end
