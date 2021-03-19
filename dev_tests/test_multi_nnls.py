#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# first, make sure the paths are set up
# we assume that this script is located and run in the folder dynamite/tests

import os
import shutil
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
import physical_system as physys

def run_user_test(make_comp=False):

    logger = logging.getLogger()
    logger.info(f'Using DYNAMITE version: {dyn.__version__}')
    logger.info(f'Located at: {dyn.__path__}')
    # print to console anyway...
    print('Using DYNAMITE version:', dyn.__version__)
    print('Located at:', dyn.__path__)

    # read configuration
    if '__file__' in globals():
        os.chdir(os.path.dirname(__file__))
    fname = 'user_test_config_multi_ml.yaml'
    c = dyn.config_reader.Configuration(fname, silent=True, reset_logging=True)

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

    # re-read configuration now that old output has been deleted
    c = dyn.config_reader.Configuration(fname, silent=True)

    logger.info(f'{c.system.n_kin} kinematics data sets in system')
    print(f'{c.system.n_kin} kinematics data sets in system')

    compare_file = outdir + "chi2_compare_ml_" \
                            f"{c.settings.orblib_settings['nE']}" \
                            f"{c.settings.orblib_settings['nI2']}" \
                            f"{c.settings.orblib_settings['nI3']}.dat"

    # "run" the models
    t = time.perf_counter()
    smi = dyn.model_iterator.ModelIterator(
        system=c.system,
        all_models=c.all_models,
        settings=c.settings,
        ncpus=c.settings.multiprocessing_settings['ncpus'])
    delt = time.perf_counter()-t
    logger.info(f'Computation time: {delt} seconds = {delt/60} minutes')
    # print to console regardless of logging level
    print(f'Computation time: {delt} seconds = {delt/60} minutes')

    # print all model results
#    c.all_models.table.pprint_all() # This only works in astropy 3.2 or later
    c.all_models.table.pprint(max_lines=-1, max_width=-1)

    if make_comp==False:
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

        # compare to chi2 in compare_file
        chi2_compare = table.Table.read(compare_file, format='ascii')
        radius = (np.max(chi2_compare['chi2']) - \
                 np.min(chi2_compare['chi2'])) / 10
        logger.debug(f'Radius={radius}')
        # print(f'Radius={radius}')
        plt.figure()
        plt.scatter(chi2_compare['model_id'],
                    chi2_compare['chi2'],
                    s=400,
                    facecolors='none',
                    edgecolors='black')
        plt.plot([i for i in range(len(c.all_models.table))],
                  c.all_models.table['chi2'],
                  'rx')
        plt.gca().set_title('calculated chi2 (red) vs should-be range '
                            '(black circles)')
        plt.xlabel('model_id')
        plt.xticks(range(len(c.all_models.table)))
        plt.ylabel('chi2')
        plt.savefig(plotfile_chi2)

        logger.info(f'Look at {plotfile_ml} and {plotfile_chi2}')
        chi2stat = ''
        for s in chi2_compare.pformat(max_lines=-1, max_width=-1):
            chi2stat += '\n'+s
        logger.info(f'chi2 comparison data: {chi2stat}')
        # print to console anyway...
        print(f'Look at {plotfile_ml} and {plotfile_chi2}')
        print('chi2 comparison data:\n')
        chi2_compare.pprint(max_lines=-1, max_width=-1)

        plotter = dyn.plotter.Plotter(system=c.system,
                                      settings=c.settings,
                                      parspace=c.parspace,
                                      all_models=c.all_models)
        stars = c.system.get_component_from_class(physys.TriaxialVisibleComponent)
        for m in range(len(c.all_models.table)):
            model = c.all_models.get_model_from_row(m)
            for kinset in range(c.system.n_kin):
                plotfile = f'{plotdir}model{m}_kin_map_{stars.kinematic_data[kinset].name}.png'
                if os.path.isfile(plotfile):
                    os.remove(plotfile)
                figure = plotter.plot_kinematic_maps(model,kin_set=kinset)
                figure.savefig(plotfile)
                print(f'Look at {plotfile}')

    return c.all_models.table, \
        compare_file, \
        c.settings.parameter_space_settings['stopping_criteria']['n_max_mods']

def create_comparison_data():

    model_results, output_file, n_max = run_user_test(make_comp=True)
    chi2_values = list(model_results['chi2'])
    for i in range(len(chi2_values),n_max): # just in case...
        chi2_values.append(float("NaN"))
    t = table.Table()
    t['model_id'] = [i for i in range(len(chi2_values))]
    t['chi2']     = chi2_values
    print(t)
    t.write(output_file, format='ascii', overwrite=True)

if __name__ == '__main__':
    run_user_test()

# end