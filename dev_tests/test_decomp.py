#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import numpy as np
# import astropy.io
import dynamite as dyn

fname = 'user_test_config_ml.yaml'
# fname = 'FCC047_2kin/FCC047_config.yaml'
c = dyn.config_reader.Configuration(fname,
                                    reset_logging=True,
                                    user_logfile='test_decomp',
                                    reset_existing_output=False)

dyn.model_iterator.ModelIterator(c) # generate models

decomp = dyn.analysis.Decomposition(config=c,
                                    kin_set=0,
                                    decomp_table=True,
                                    comps_weights=True) # do the decomposition

for v_sigma_option in ('moments', 'fit'):
    #select the components and calculate the kinematics for each
    #(this is done with the selection used in Santucci+22)
    #and plot the kinematics
    decomp.plot_decomp(xlim=15,
                       ylim=15,
                       v_sigma_option=v_sigma_option,
                       comps_plot={'thin_d': True, 'thick_d': True,
                                   'disk': True, 'cr_thin_d': True,
                                   'cr_thick_d': True, 'cr_disk': True,
                                   'bulge': True, 'all': True},
                       individual_colorbars=True)
