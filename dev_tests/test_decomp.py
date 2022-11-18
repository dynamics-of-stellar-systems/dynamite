#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import numpy as np
# import astropy.io
import dynamite as dyn

fname = 'user_test_config_ml.yaml'
c = dyn.config_reader.Configuration(fname,
                                    reset_logging=True,
                                    user_logfile='test_decomp',
                                    reset_existing_output=False)

decomp = dyn.decomposition.Decomp(c,read_orblib='dynamite')

for conversion in ('gh_fit_with_free_v_sigma_params',
                   'gh_expand_around_losvd_mean_and_std_deviation',
                   'gh_fit_with_free_v_sigma_params_fortran', 'moments'):
#for conversion in ('gh_fit_with_free_v_sigma_params',):
    #select the components and calculate the kinematics for each
    #(this is done with the selection used in Santucci+22)
    #and plot the kinematics
    decomp.plot_decomp(xlim=15, ylim=15, conversion=conversion)
