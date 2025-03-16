#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time

# Set matplotlib backend to 'Agg' (compatible when X11 is not running
# e.g., on a cluster). Note that the backend can only be set BEFORE
# matplotlib is used or even submodules are imported!
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from astropy import table
import dynamite as dyn
import dynamite.constants as const

def run_user_test(make_comp=False):

    print('Using DYNAMITE version:', dyn.__version__)
    print('Located at:', dyn.__path__)

    # read configuration
    if '__file__' in globals():
        file_dir = os.path.dirname(__file__)
        if file_dir:
            os.chdir(file_dir)
    else:
        file_dir = None
    fname = 'user_test_config_ml_gas.yaml'
    # fname = 'FCC047_2kin/FCC047_config.yaml'
    # fname = 'user_test_config_fixedvalues.yaml'
    c = dyn.config_reader.Configuration(fname,
                                        reset_logging=True,
                                        user_logfile='test_nnls',
                                        reset_existing_output=True)

    # delete previous output if available
    # c.remove_existing_orblibs()
    # c.remove_existing_all_models_file(wipe_other_files=False)
    # c.backup_config_file(keep=3, delete_other=True)
    # c.remove_existing_plots()

    print(f'{const.GRAV_CONST_KM=}, {const.PARSEC_KM=}, {const.RHO_CRIT=}.')

    which_chi2 = c.settings.parameter_space_settings['which_chi2']

    plotdir = c.settings.io_settings['plot_directory']
    plotfile_ml = plotdir + f'ml_vs_iter_{which_chi2}.png'
    if os.path.isfile(plotfile_ml):
        os.remove(plotfile_ml)
    plotfile_chi2 = plotdir + f'{which_chi2}_vs_model_id.png'
    if os.path.isfile(plotfile_chi2):
        os.remove(plotfile_chi2)

    compare_file = file_dir if file_dir else '.'
    compare_file += "/data/chi2_compare_ml_" \
                    f"{c.settings.orblib_settings['nE']}" \
                    f"{c.settings.orblib_settings['nI2']}" \
                    f"{c.settings.orblib_settings['nI3']}.dat"
    print(f'Comparing results to {compare_file}.')
    # "run" the models
    t = time.perf_counter()
    smi = dyn.model_iterator.ModelIterator(c)
    delt = time.perf_counter()-t
    print(f'Computation time: {delt} seconds = {delt/60} minutes')

    # print all model results
#    c.all_models.table.pprint_all() # This only works in astropy 3.2 or later
    c.all_models.table.pprint(max_lines=-1, max_width=-1)

    if not make_comp:
        # plot the models
        plt.figure()
        plt.scatter(c.all_models.table['which_iter'],
                    c.all_models.table['ml'],
                    c=c.all_models.table[which_chi2],
                    cmap=plt.cm.viridis_r,
                    s=200)
        cb=plt.colorbar()
        cb.set_label(which_chi2, y=1.1, labelpad=-40, rotation=0)
        plt.gca().set_title('all iterations')
        # plt.gca().set_yscale('log')
        plt.xlabel('iteration')
        plt.xticks(range(0,max(c.all_models.table['which_iter'])+1))
        plt.ylabel('ml')
        plt.savefig(plotfile_ml)

        # compare to chi2 in compare_file
        chi2_compare = table.Table.read(compare_file, format='ascii')
        plt.figure()
        plt.scatter(chi2_compare['model_id'],
                    chi2_compare[which_chi2],
                    s=2000,
                    facecolors='none',
                    edgecolors='black')
        plt.plot(range(len(c.all_models.table)),
                  c.all_models.table[which_chi2],
                  'rx')
        plt.gca().set_title(f'calculated {which_chi2} (red) vs '
                            'should-be range (black circles)')
        plt.xlabel('model_id')
        plt.xticks(range(len(c.all_models.table)))
        plt.ylabel(which_chi2)
        plt.savefig(plotfile_chi2)

        print(f'Look at {plotfile_ml} and {plotfile_chi2}')
        print(f'{which_chi2} comparison data:\n')
        chi2_compare.pprint(max_lines=-1, max_width=-1)
        print('The best 2 models:')
        c.all_models.get_best_n_models(n=2).pprint(max_lines=-1, max_width=-1)

    if make_comp:
        return c.all_models.table, \
          compare_file, \
          c.settings.parameter_space_settings['stopping_criteria']['n_max_mods']
    else:
        return c

def create_comparison_data():

    model_results, output_file, n_max = run_user_test(make_comp=True)
    t = table.Table()
    t['model_id'] = range(n_max)
    for which_chi2 in ('chi2', 'kinchi2', 'kinmapchi2'):
        values = list(model_results[which_chi2])
        for i in range(len(values),n_max): # just in case...
            values.append(float("NaN"))
        t[which_chi2] = values
    print(t)
    t.write(output_file, format='ascii', overwrite=True)

if __name__ == '__main__':
    # create_comparison_data()
    c = run_user_test()

# end
