#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import numpy as np
import time
import matplotlib.pyplot as plt
from astropy import table
import dynamite as dyn

def plot_losvds(losvd_histogram,
                orb_idx,
                aperture_idx_list,
                ls='-',
                color='k',
                ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    v = losvd_histogram.x
    losvd = losvd_histogram.y[orb_idx, :, :]
    plt.plot(v, np.sum(losvd, 1), color=color, ls=ls, label='total')
    for aperture_idx in aperture_idx_list:
        plt.plot(v,
                 losvd[:, aperture_idx],
                 ls=ls,
                 color=color,
                 label=f'aperture {aperture_idx}')
    plt.gca().set_title(f'LOSVD of orbit {orb_idx}')
    plt.gca().set_xlabel('v [km/s]')
    plt.gca().set_yscale('log')
    plt.gca().legend()
    plt.tight_layout()
    return ax

def run_orbit_losvd_test(make_comparison_losvd=False):

    print('Using DYNAMITE version:', dyn.__version__)
    print('Located at:', dyn.__path__)

    # read configuration
    fname = 'user_test_config.yaml'
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
    plotfile = plotdir + 'orbit_losvds.png'
    if os.path.isfile(plotfile):
        os.remove(plotfile)

    parset = c.parspace.get_parset()
    model = dyn.model.LegacySchwarzschildModel(
        system=c.system,
        settings=c.settings,
        parspace=c.parspace,
        executor=c.executor,
        parset=parset)
    model.setup_directories()
    model.get_model_directory()
    model.get_orblib()
    model.orblib.read_losvd_histograms()

    fname = outdir + 'comparison_losvd.npz'
    if make_comparison_losvd:
        np.savez(fname,
                 xedg=model.orblib.losvd_histograms.xedg,
                 y=model.orblib.losvd_histograms.y)
    tmp = np.load(fname)
    comparison_losvd = dyn.kinematics.Histogram(xedg=tmp['xedg'],
                                                y=tmp['y'],
                                                normalise=False)

    # read orbits and plot them
    if make_comparison_losvd is False:
        orb_idx = 15
        aperture_idx_list = [0, 2, 20, 30]
        ax = plot_losvds(comparison_losvd,
                         orb_idx,
                         aperture_idx_list)
        ax = plot_losvds(model.orblib.losvd_histograms,
                         orb_idx,
                         aperture_idx_list,
                         ls='--',
                         color='r',
                         ax=ax)
        fig = plt.gcf()
        fig.savefig(plotfile)

    print(f'Look at {plotfile}')

    return 0

if __name__ == '__main__':
    run_orbit_losvd_test()

# end
