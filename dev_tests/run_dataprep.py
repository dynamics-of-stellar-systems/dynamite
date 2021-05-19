#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 21:52:29 2021

@author: sabine
"""
import os
import dynamite as dyn
from dynamite.data_prep.generate_kin_input import create_kin_input

dir='Data_prep/Kinematics/'

input=dir+'CALIFA/NGC6278.V1200.rscube_INDOUSv2_SN20_stellar_kin.fits'

outdir='Data_prep/Dyn_input/CALIFA/'
expr=''


#better use
#dyn.model.Model.create_model_directory(output)
if not os.path.exists(outdir):
    os.makedirs(outdir)

create_kin_input('NGC6278',input, outdir, expr='', fit_PA=True, kin_input='CALIFA')

#add PSF to kinfile
gh = dyn.kinematics.GaussHermite() 
gh.add_psf_to_datafile(sigma=[1.06], weight=[1.0], datafile=outdir+'gauss_hermite_kins'+expr+'.ecsv')
