#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# first, make sure the paths are set up
# we assume that this script is located and run in the folder dynamite/tests

import os
import sys
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
    fname = 'user_test_config_ml.yaml'
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
    plotfile_ml = plotdir + 'ml_vs_iter_chi2.png'
    if os.path.isfile(plotfile_ml):
        os.remove(plotfile_ml)
    plotfile_chi2 = plotdir + 'chi2_vs_model_id.png'
    if os.path.isfile(plotfile_chi2):
        os.remove(plotfile_chi2)

    stat_file = outdir + "chi2_compare_ml_" \
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
        plt.figure()
        plt.scatter(c.all_models.table['which_iter'],
                    c.all_models.table['ml'],
                    c=c.all_models.table['chi2'],
                    cmap=plt.cm.viridis_r,
                    s=200)
        cb=plt.colorbar()
        cb.set_label('chi2', y=1.1, labelpad=-40, rotation=0)
        plt.gca().set_title('all iterations')
        # plt.gca().set_yscale('log')
        plt.xlabel('iteration')
        plt.xticks(range(0,max(c.all_models.table['which_iter'])+1))
        plt.ylabel('ml')
        plt.savefig(plotfile_ml)

        # compare to chi2 in chi2_compare.dat
        chi2_compare = table.Table.read(stat_file, format='ascii')
        radius = (np.max(chi2_compare['chi2_average']) - \
                 np.min(chi2_compare['chi2_average'])) / 10
        print(f'Radius={radius}')
        plt.figure()
        plt.scatter(chi2_compare['model_id'],
                    chi2_compare['chi2_average'],
                    s=400,
                    facecolors='none',
                    edgecolors='black')
        # plt.vlines(chi2_compare['model_id'],
        #            chi2_compare['chi2_min'],
        #            chi2_compare['chi2_max'])
        plt.plot([i for i in range(len(c.all_models.table))],
                  c.all_models.table['chi2'],
                  'rx')
        plt.gca().set_title('calculated chi2 (red) vs should-be range '
                            '(black circles)')
        plt.xlabel('model_id')
        plt.xticks(range(len(c.all_models.table)))
        plt.ylabel('chi2')
        plt.savefig(plotfile_chi2)

        print(f'Look at {plotfile_ml} and {plotfile_chi2}')
        print('chi2 statistics for comparison:\n')
        print(chi2_compare.pprint(max_lines=-1, max_width=-1))

    return c.all_models.table, \
        stat_file, \
        c.settings.parameter_space_settings['stopping_criteria']['n_max_mods']

def create_stats(n_chi2=10):
    chi2_all = []

    # run user_test n_chi2 times
    for i in range(n_chi2):
        model_results, output_file, n_max = run_user_test(stat_mode=True)
        chi2_all.append([])
        chi2_values = list(model_results['chi2'])
        for j in range(len(chi2_values)):
            chi2_all[i].append(chi2_values[j])
        # if j<n_max:
        #     chi2_all[i].append(float("NaN"))

    chi2_all = np.array(chi2_all).T

    t = table.Table()
    t['model_id'] = [i for i in range(len(chi2_all))]
    t['chi2_min'] = np.nanmin(chi2_all, axis=1)
    t['chi2_max'] = np.nanmax(chi2_all, axis=1)
    t['chi2_average'] = np.nanmean(chi2_all, axis=1)
    t['chi2_sdev'] = np.nanstd(chi2_all, axis=1)

    print(t)

    np.savetxt(output_file+f'_raw_{n_chi2}.dat', chi2_all)

    t.write(output_file+f'_{n_chi2}.dat', format='ascii')


if __name__ == '__main__':
    run_user_test()

# end