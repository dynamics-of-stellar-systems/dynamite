#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 08:54:48 2021

@author: sabine
"""

import dynamite as dyn
import os
import shutil


print('DYNAMITE')
print('    version', dyn.__version__)
print('    installed at ', dyn.__path__)

fname = 'FCC047_2kin/FCC047_config.yaml'

c = dyn.config_reader.Configuration(fname, reset_logging=True)

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

print(type(c.system))
print(type(c.settings))

parset = c.parspace.get_parset()
print(parset)

nkin = c.system.n_kin
print(f'{nkin} kinematics data sets in system')

model = dyn.model.LegacySchwarzschildModel(
    system=c.system,
    settings=c.settings,
    parspace=c.parspace,
    parset=parset)




model.setup_directories()

orblib=model.get_orblib()

model.get_weights(orblib=orblib)

print(model.chi2)


plotter = dyn.plotter.Plotter(system=c.system,
                              settings=c.settings,
                              parspace=c.parspace,
                              all_models=c.all_models)

for kinset in range(nkin):
    plotfile = f'{plotdir}kin_map{kinset}.png'
    if os.path.isfile(plotfile):
        os.remove(plotfile)
    figure = plotter.plot_kinematic_maps(model,kin_set=kinset+1)
    figure.savefig(plotfile)
    print(f'Look at {plotfile}')
