#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 08:54:48 2021

@author: sabine
"""

import os
import dynamite as dyn
from dynamite import physical_system as physys

print('DYNAMITE')
print('    version', dyn.__version__)
print('    installed at ', dyn.__path__)

fname = 'FCC047_2kin/FCC047_config.yaml'

c = dyn.config_reader.Configuration(fname, reset_logging=True)
c.remove_existing_orblibs()
c.remove_existing_all_models_file()
plotdir = c.settings.io_settings['plot_directory']

print(type(c.system))
print(type(c.settings))

parset = c.parspace.get_parset()
print(parset)

nkin = c.system.n_kin
print(f'{nkin} kinematics data sets in system')

model = dyn.model.Model(
    system=c.system,
    settings=c.settings,
    parspace=c.parspace,
    parset=parset)

model.setup_directories()

orblib=model.get_orblib()

model.get_weights(orblib=orblib)

print(model.chi2)

plotter = dyn.plotter.Plotter(config=c)

stars = c.system.get_component_from_class(physys.TriaxialVisibleComponent)
for kinset in range(nkin):
    plotfile = f'{plotdir}kin_map_{stars.kinematic_data[kinset].name}.png'
    if os.path.isfile(plotfile):
        os.remove(plotfile)
    figure = plotter.plot_kinematic_maps(model,kin_set=kinset)
    figure.savefig(plotfile)
    print(f'Look at {plotfile}')
