#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This script can be called while a DYNAMITE job is running and it will plot
# (kin)chi2 vs model id, chi2 maps, and kinematic maps.
# Call from the directory your config file is in.

import os
import logging
import numpy as np

# Set matplotlib backend to 'Agg' (compatible when X11 is not running
# e.g., on a cluster). Note that the backend can only be set BEFORE
# matplotlib is used or even submodules are imported!
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import dynamite as dyn
import plotter

logger = logging.getLogger()
logger.info(f'Using DYNAMITE version: {dyn.__version__}')
logger.info(f'Located at: {dyn.__path__}')
# print to console anyway...
print('Using DYNAMITE version:', dyn.__version__)
print('Located at:', dyn.__path__)

# read configuration
# if '__file__' in globals():
#     os.chdir(os.path.dirname(__file__))
fname = 'GAMA30890_config.yaml'
c = dyn.config_reader.Configuration(fname, reset_logging=False)

io_settings = c.settings.io_settings
outdir = io_settings['output_directory']
plotdir = outdir + 'plots/'
if not os.path.isdir(plotdir):
    os.mkdir(plotdir)
plotfile_chi2 = plotdir + 'chi2_vs_model_id-in-progress.png'
plotfile_triangle = plotdir + 'chi2_plot-in-progress.png'
plotfile_kinmap = plotdir + 'kinmap-in-progress.png'
for f in [plotfile_chi2, plotfile_triangle, plotfile_kinmap]:
    if os.path.isfile(f):
        os.remove(f)

c.all_models.table.pprint(max_lines=-1, max_width=-1)
text = f'Number of models in all_models table: {len(c.all_models.table)}'
logger.info(text)
print(text)
which_chi2 = c.settings.parameter_space_settings['which_chi2']
models_done = np.where(c.all_models.table['all_done'])
min_chi2 = min(m[which_chi2] for m in c.all_models.table[models_done])
c.all_models.table.add_index(which_chi2)
model_id = c.all_models.table.loc_indices[min_chi2]
model = c.all_models.get_model_from_row(model_id)
text = f'Optimal model parset so far: {model_id}, which_chi2={which_chi2}\n' \
       f'{model.parset}\n{c.all_models.table[model_id]}'
logger.info(text)
print(text)

plt.figure()
plt.plot([i for i in range(len(c.all_models.table))],
         c.all_models.table[which_chi2],
         'rx')
plt.gca().set_title(f'calculated {which_chi2} (red) vs model id')
plt.xlabel('model_id')
plt.xticks(range(len(c.all_models.table)))
plt.ylabel(f'{which_chi2}')
plt.savefig(plotfile_chi2)

theplotter = plotter.Plotter(system=c.system,
                             settings=c.settings,
                             parspace=c.parspace,
                             all_models=c.all_models
                            )
chi2_plot = theplotter.make_chi2_plot()
chi2_plot.savefig(plotfile_triangle)

kinmap = theplotter.plot_kinematic_maps(cbar_lims='data')
kinmap.savefig(plotfile_kinmap)

text = f'Look at {plotfile_chi2}, {plotfile_triangle}, and {plotfile_kinmap}'
logger.info(text)
print(text)

# end