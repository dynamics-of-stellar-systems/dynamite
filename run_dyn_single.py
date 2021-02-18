#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 08:54:48 2021

@author: sabine
"""

import dynamite as dyn


print('DYNAMITE')
print('    version', dyn.__version__)
print('    installed at ', dyn.__path__)

fname = 'FCC047_2kin/FCC047_config.yaml'

c = dyn.config_reader.Configuration(fname)

print(type(c.system))
print(type(c.settings))

parset = c.parspace.get_parset()
print(parset)




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

figure = plotter.plot_kinematic_maps(model,kin_set=1)



