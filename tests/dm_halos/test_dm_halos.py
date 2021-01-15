#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil

import time
import numpy as np
import matplotlib.pyplot as plt
from astropy import table
import dynamite as dyn

def run_dm_test(config_file=None, remove_old_output=True):

    print('Using DYNAMITE version:', dyn.__version__)
    print('Located at:', dyn.__path__)
    # read config file
    c = dyn.config_reader.Configuration(config_file, silent=True)
    if remove_old_output:
        # remove all previous output
        io_settings = c.settings.io_settings
        outdir = io_settings['output_directory']
        if os.path.isdir(outdir):
            shutil.rmtree(outdir, ignore_errors=True)
        # re-read config file after previous output removed
        c = dyn.config_reader.Configuration(config_file, silent=True)
        # run the models
        t = time.perf_counter()
        smi = dyn.model_iterator.ModelIterator(
            system=c.system,
            all_models=c.all_models,
            settings=c.settings,
            executor=c.executor)
        delt = time.perf_counter()-t
        print(f'Computation time: {delt} seconds = {delt/60} minutes')
    # print all model results
    c.all_models.table.pprint(max_lines=-1, max_width=-1)
    # plot LOSVD
    n_mods = len(c.all_models.table)
    fig, ax = plt.subplots(1, 1)
    # choose an orbit and aperture
    idx_orb = 98
    idx_ap = 25
    # get the appropriate mass-scaling parameter to label the plot
    if 'f_dark_halo' in c.all_models.table.colnames:
        label = 'f_dark_halo'
        leg_title = 'fraction [M_halo/M_*](<R_200)'
    elif 'Mvir_dark_halo' in c.all_models.table.colnames:
        label = 'Mvir_dark_halo'
        leg_title = 'Mvir_halo'
    else:
        raise ValueError('mass-scaling parameter of dark halo not found')
    for i in range(n_mods):
        # get the orbit library for this model
        parset_i = c.all_models.get_parset_from_row(i)
        mod_i = dyn.model.LegacySchwarzschildModel(
            system=c.system,
            settings=c.settings,
            executor=c.executor,
            parspace=c.parspace,
            parset=parset_i)
        orblib_i = mod_i.get_orblib()
        mod_i.orblib.read_losvd_histograms()
        plt.plot(mod_i.orblib.losvd_histograms.x,
                 mod_i.orblib.losvd_histograms.y[idx_orb,:,idx_ap],
                 label=parset_i[label])
    leg = ax.legend()
    leg.set_title(leg_title)
    ax.set_xlabel('$v_{LOS}$ [km/s]')
    ax.set_ylabel('')
    ax.set_title(f'LOSVDs of orbit {idx_orb} in aperture {idx_ap}')
    fig.tight_layout()
    outfile = c.settings.io_settings['output_directory'] + 'losvd_plot.png'
    fig.savefig(outfile)
    plt.close()

    return

if __name__ == '__main__':

    # config files have all parameters fixed except mass-scaling of DM halo
    # the `run_dm_test` function runs 3 models with low, middle, hi DM masses
    # and plots the LOSVD of one orbit in one aperture

    # run NFW
    config_file = 'nfw_config.yaml'
    run_dm_test(config_file=config_file, remove_old_output=False)

    # run gNFW
    config_file = 'gnfw_config.yaml'
    run_dm_test(config_file=config_file, remove_old_output=False)



# end
