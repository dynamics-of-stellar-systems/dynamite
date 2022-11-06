#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import numpy as np
# import astropy.io
from dynamite.decomposition import Decomp

# w_dir='/Users/z5178033/Downloads/' # the working directory where you have all the galaxy catalogue
w_dir='../../Giulia/decomposition/'

# gal_infos=astropy.io.ascii.read(w_dir+'SAMI_gals_cat.dat')
# gals=gal_infos['CATAID']


gal=41059
print('DOING GALAXY ', gal)

# i = np.where(gal_infos['CATAID'] == gal)[0][0]
# Re=gal_infos['Re[arcs]'][i]

#folder where you have the galaxies folders with the data
# w_dir1='/Users/z5178033/Downloads/test_crcut/'
w_dir1=w_dir

#TIM START
# #read the orbits and create the velocity histrogram
# losvd_histograms, proj_mass=run_dec(str(gal), w_dir1, Re)
# print('Orbits read')
# comps=['disk', 'thin_d', 'warm_d', 'bulge', 'all']
#TIM END
gal = str(gal)
decomp = Decomp(galaxy = gal,
                input_directory = w_dir1 + gal + '/',
                output_directory = w_dir1 + gal + '/')


for conversion in ('gh_fit', 'losvd_vsig', 'fortran', 'moments'):
    #select the components and calculate the kinematics for each
    #(this is done with the selection used in Santucci+22)
    decomp.comps_aphist(conversion)
    print('Components done')
    #plot the kinematics
    decomp.plot_comps_giu(gal=gal, xlim=15, ylim=15, conversion=conversion)
    print('Plots done')
print('**************************')
