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

decomp = dyn.decomposition.Decomp(c,read_orblib='decomp')

# for conversion in ('gh_fit', 'losvd_vsig', 'fortran', 'moments'):
for conversion in ('losvd_vsig', 'gh_fit', 'fortran', 'moments'):
    #select the components and calculate the kinematics for each
    #(this is done with the selection used in Santucci+22)
    decomp.comps_aphist(conversion)
    print('Components done')
    #plot the kinematics
    decomp.plot_comps_giu(xlim=15, ylim=15, conversion=conversion)
    print('Plots done')
