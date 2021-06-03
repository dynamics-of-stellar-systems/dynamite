#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import logging
import time

# Set matplotlib backend to 'Agg' (compatible when X11 is not running
# e.g., on a cluster). Note that the backend can only be set BEFORE
# matplotlib is used or even submodules are imported!
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import dynamite as dyn

def run_user_test():

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
    fname = 'user_test_config_multi_ml_FCC.yaml'
    c = dyn.config_reader.Configuration(fname, reset_logging=True)

    c.remove_existing_orblibs()
    c.remove_existing_all_models_file()
    plotdir = c.settings.io_settings['plot_directory']
    plotfile_ml = plotdir + 'ml_vs_iter_chi2.png'
    if os.path.isfile(plotfile_ml):
        os.remove(plotfile_ml)
    plotfile_chi2 = plotdir + 'chi2_vs_model_id.png'
    if os.path.isfile(plotfile_chi2):
        os.remove(plotfile_chi2)

    logger.info(f'{c.system.n_kin} kinematics data sets in system')
    print(f'{c.system.n_kin} kinematics data sets in system')

    compare_file = os.path.dirname(__file__) \
                    + "/data/chi2_compare_ml_" \
                      f"{c.settings.orblib_settings['nE']}" \
                      f"{c.settings.orblib_settings['nI2']}" \
                      f"{c.settings.orblib_settings['nI3']}.dat"

    # "run" the models
    t = time.perf_counter()
    smi = dyn.model_iterator.ModelIterator(c)
    delt = time.perf_counter()-t
    logger.info(f'Computation time: {delt} seconds = {delt/60} minutes')
    # print to console regardless of logging level
    print(f'Computation time: {delt} seconds = {delt/60} minutes')

    c.all_models.table.pprint(max_lines=-1, max_width=-1)

    # for one model, re-calculate solution with the new weight solver
    print('Recalculating orbit weights with scipy NNLS solver')

    fig, ax = plt.subplots(1, 3, sharey=True, figsize=(12,4))
    for i in [0,1,2]:
        mod0 = c.all_models.get_model_from_row(i)
        parset0 = c.all_models.get_parset_from_row(i)
        orblib0 = dyn.orblib.LegacyOrbitLibrary(config=c,
                                                mod_dir=mod0.directory_noml,
                                                parset=parset0)
        orblib0.read_losvd_histograms()
        weight_solver = mod0.get_weights()
        weights_old, chi2_tot_old, chi2_kin_old = weight_solver.solve(orblib0)
        weight_solver_new = dyn.weight_solvers.NNLS(
                config=c,
                directory_with_ml=mod0.directory,
                CRcut=True,
                nnls_solver='scipy')
        solution_new = weight_solver_new.solve(orblib0)
        weights_new = solution_new[0]
        ax[i].plot(weights_old, label='Legacy')
        ax[i].plot(weights_new, '--', label='NNLS scipy')
        ax[i].legend()
        ax[i].set_xlabel('orbit')
        ax[i].set_title(f'Model {i}')
    ax[0].set_ylabel('weight')
    fig.subplots_adjust(wspace=0)
    fig.tight_layout()
    plotfile = f'{plotdir}multikin_wsolver_compare.png'
    if os.path.isfile(plotfile):
        os.remove(plotfile)
    fig.savefig(plotfile)
    logger.info(f'Look at {plotfile}')
    print(f'Look at {plotfile}')

    return c.all_models.table, \
        compare_file, \
        c.settings.parameter_space_settings['stopping_criteria']['n_max_mods']

if __name__ == '__main__':
    run_user_test()

# end
