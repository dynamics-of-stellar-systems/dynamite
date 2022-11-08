#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import numpy as np
# import astropy.io
from dynamite.decomposition import Decomp

choice = input("Enter 1 for Giulia's galaxy, 2 for NGC6278: ")
if choice == '1':
    in_dir = out_dir = '../../Giulia/decomposition/41059DYN/'
    model_dir = 'bh8.75dc1.15f0.90q0.36p0.990u0.9999/ml7.91/'
elif choice == '2':
    in_dir = 'NGC6278_input/'
    out_dir = 'NGC6278_output/'
    model_dir = 'orblib_000_000/ml1.00/'
else:
    raise ValueError("That's not 1 and not 2...")

decomp = Decomp(input_directory = in_dir,
                output_directory = out_dir,
                model = model_dir)

for conversion in ('gh_fit', 'losvd_vsig', 'fortran', 'moments'):
    #select the components and calculate the kinematics for each
    #(this is done with the selection used in Santucci+22)
    decomp.comps_aphist(conversion)
    print('Components done')
    #plot the kinematics
    decomp.plot_comps_giu(xlim=15, ylim=15, conversion=conversion)
    print('Plots done')
print('**************************')
