#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import numpy as np
# import astropy.io
import dynamite as dyn

fname = 'user_test_config_ml.yaml'
#fname = 'FCC047_2kin/FCC047_config.yaml'
c = dyn.config_reader.Configuration(fname,
                                    reset_logging=True,
                                    user_logfile='test_decomp',
                                    reset_existing_output=False)

dyn.model_iterator.ModelIterator(c) # generate models

decomp = dyn.orbit_exploration.Decomposition(config=c,
                                             kin_set=0) # do the decomposition

for conversion in (
                   'gh_expand_around_losvd_mean_and_std',
                   'gh_fit_with_free_v_sigma_params_fortran',
                   'moments',
                   'moments_old',
                   'fit',
                   'gh_fit_with_free_v_sigma_params'
                   ):
    #select the components and calculate the kinematics for each
    #(this is done with the selection used in Santucci+22)
    #and plot the kinematics
    decomp.plot_decomp(xlim=15, ylim=15, conversion=conversion)
