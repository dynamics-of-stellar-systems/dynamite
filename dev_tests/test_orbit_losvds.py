#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import logging
# import importlib
import numpy as np

# Set matplotlib backend to 'Agg' (compatible when X11 is not running
# e.g., on a cluster). Note that the backend can only be set BEFORE
# matplotlib is used or even submodules are imported!
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
#from astropy import table
import dynamite as dyn

def version_p():
    return sys.version.split()[0]

def version_f():
    v = subprocess.run("gfortran --version", capture_output=True, shell=True, \
        check=True).stdout.decode('utf-8').split(sep='\n')[0].split()[-1]
    return v

def plot_losvds(losvd_histogram,
                orb_idx,
                aperture_idx_list,
                ls='-',
                color='k',
                ax=None,
                title=False,
                seed=None):
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
                 label=f'bin {aperture_idx}')
    if title:
        plt.gca().set_title(f'LOSVD of orbit {orb_idx}\n'
                            f'Python {version_p()}, gfortran {version_f()}\n'
                            f'Random seed: {seed}')
    plt.gca().set_xlabel('v [km/s]')
    plt.gca().set_yscale('log')
    plt.gca().legend(loc='upper right')
    plt.xlim([-1750, 1750])
    plt.tight_layout()
    return ax

def run_orbit_losvd_test(make_comparison_losvd=False):

    logging.info(f'Using DYNAMITE version: {dyn.__version__}')
    logging.info(f'Located at: {dyn.__path__}')

    # read configuration
    fname = 'user_test_config.yaml'
    c = dyn.config_reader.Configuration(fname,
                                        reset_logging=False,
                                        reset_existing_output=True)

    # c.remove_existing_orblibs()
    # c.remove_existing_all_models_file()
    # c.backup_config_file(reset=True)

    plotdir = c.settings.io_settings['plot_directory']
    plotfile = plotdir + f'orbit_losvds-{version_p()}-{version_f()}.png'
    if os.path.isfile(plotfile):
        os.remove(plotfile)

    parset = c.parspace.get_parset()
    model = dyn.model.Model(config=c, parset=parset)
    model.setup_directories()
    orbit_library = model.get_orblib()
    orbit_library.read_losvd_histograms()

    file_dir = os.path.dirname(__file__)
    fname = file_dir if file_dir else '.'
    fname += '/data/comparison_losvd.npz'
    if make_comparison_losvd:
        # create comparison file
        np.savez_compressed(fname,
                            xedg=orbit_library.losvd_histograms[0].xedg,
                            y=orbit_library.losvd_histograms[0].y)

        logging.info(f'Losvd comparison file {fname} created')
    else:
        # read orbits and plot them
        tmp = np.load(fname)
        comparison_losvd = dyn.kinematics.Histogram(xedg=tmp['xedg'],
                                                    y=tmp['y'],
                                                    normalise=False)
        orb_idx = 15
        aperture_idx_list = [0, 2, 20, 30]
        ax = plot_losvds(comparison_losvd,
                         orb_idx,
                         aperture_idx_list)
        ax = plot_losvds(orbit_library.losvd_histograms[0],
                         orb_idx,
                         aperture_idx_list,
                         ls='--',
                         color='r',
                         ax=ax,
                         title=True,
                         seed=c.settings.orblib_settings['random_seed'])
        fig = plt.gcf()
        fig.savefig(plotfile)

        logging.info(f'Look at {plotfile}')
        # we want to print to the console regardless of the logging level
        print(f'Look at {plotfile}')

    return 0

if __name__ == '__main__':

    # For an ultra-minimal logging configuration, not even the
    # logging.basicConfig call is necessary, but here we want to log on the
    # INFO level and delete all other logging settings.
    # Comment out the following logging lines to see warnings only or - if
    # logging was already configured in your kernel - use the existing
    # logging settings.
    ## logging.shutdown()         # Only needed to reset existing logging setup
    ## importlib.reload(logging)  # Only needed to reset existing logging setup
    logging.basicConfig(level=logging.INFO)

    run_orbit_losvd_test()

# end
