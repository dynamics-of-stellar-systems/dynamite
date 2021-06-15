#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import time
import numpy as np
import matplotlib.pyplot as plt
import dynamite as dyn
from dynamite import physical_system as physys

def remove_existing_output(config, remove_orblibs=False):
    # delete model directory if it exits
    if remove_orblibs:
        config.remove_existing_orblibs()
    else:
        # just remove all 'ml' directories
        config.remove_existing_orbital_weights()
    # delete the all_models file
    config.remove_existing_all_models_file()

def run_user_test():

    print('Using DYNAMITE version:', dyn.__version__)
    print('Located at:', dyn.__path__)

    # read configuration
    # this file uses the old 'LegcayWeight' weight solver
    fname = 'reimplement_nnls_config1.yaml'
    c1 = dyn.config_reader.Configuration(fname, reset_logging=True)
    remove_existing_output(c1, remove_orblibs=True)

    # run the models
    t = time.perf_counter()
    smi = dyn.model_iterator.ModelIterator(
        system=c1.system,
        all_models=c1.all_models,
        settings=c1.settings,
        ncpus=c1.settings.multiprocessing_settings['ncpus'])
    delt = time.perf_counter()-t
    print(f'Computation time: {delt} seconds = {delt/60} minutes')
    c1.all_models.table.pprint(max_lines=-1, max_width=-1)

    # read 2nd configuration file
    # this file uses the new 'NNLS' weight solver
    fname = 'reimplement_nnls_config2.yaml'
    c2 = dyn.config_reader.Configuration(fname)
    remove_existing_output(c2, remove_orblibs=False)

    # "run" the models
    t = time.perf_counter()
    smi = dyn.model_iterator.ModelIterator(
        system=c2.system,
        all_models=c2.all_models,
        settings=c2.settings,
        ncpus=c2.settings.multiprocessing_settings['ncpus'])
    delt = time.perf_counter()-t
    print(f'Computation time: {delt} seconds = {delt/60} minutes')
    c2.all_models.table.pprint(max_lines=-1, max_width=-1)

    # plot result
    outdir = c1.settings.io_settings['output_directory']
    plotdir = outdir + 'plots/'
    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)

    fig, ax = plt.subplots(1, 3, figsize=(9, 3))
    # plot weights from one model using LegacyWeightSolvery
    model_idx = 0
    mod1 = c1.all_models.get_model_from_row(model_idx)
    mod1.get_weights()
    ax[0].plot(mod1.weights, label='LegacyWeightSolver')
    # plot weights from model 0 using NNLS
    mod2 = c2.all_models.get_model_from_row(model_idx)
    orblib2 = mod2.get_orblib()
    mod2.get_weights(orblib2)
    ax[0].plot(mod2.weights, ':', label='NNLS')
    all_close = np.allclose(mod1.weights, mod2.weights, rtol=1e-10, atol=1e-6)
    if all_close:
        ax[0].set_title(f'weights of model {model_idx} ARE identical')
    else:
        ax[0].set_title(f'weights of model {model_idx} ARE NOT identical')
    ax[0].legend()
    ax[0].set_xlabel('orbit')
    ax[0].set_ylabel('weight')
    # plot chi2 grid using LegacyWeightSolver
    c1amt = c1.all_models.table
    dark_halo = c1.system.get_component_from_class(physys.NFW).name
    ax[1].scatter(c1amt[f'f-{dark_halo}'],
                  c1amt['ml'],
                  s=40,
                  c=c1amt['kinchi2'])
    ax[1].set_xscale('log')
    ax[1].set_title('LegacyWeightSolver $\chi^2$')
    ax[1].set_xlabel('$f_{DM}$')
    ax[1].set_ylabel('ML')
    # plot chi2 grid using new solver
    c2amt = c2.all_models.table
    ax[2].scatter(c2amt[f'f-{dark_halo}'],
                  c2amt['ml'],
                  s=40,
                  c=c2amt['kinchi2'])
    ax[2].set_xscale('log')
    all_close = np.allclose(c1amt['kinchi2'],
                            c2amt['kinchi2'],
                            rtol=1e-10,
                            atol=1e-6)
    if all_close:
        ax[2].set_title('new $\chi^2$ ARE IDENTICAL')
    else:
        ax[2].set_title('new $\chi^2$ ARE NOT IDENTICAL')
    ax[2].set_xlabel('$f_{DM}$')
    ax[2].set_ylabel('ML')
    fig.tight_layout()
    fig.savefig(f'{plotdir}compare_weight_solver.png')
    plt.close()

if __name__ == '__main__':
    run_user_test()

# end
