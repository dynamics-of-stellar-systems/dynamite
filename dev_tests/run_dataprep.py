#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""i
Created on Sun May 16 21:52:29 2021

@author: sabine
"""
import os
import dynamite as dyn
from dynamite.data_prep.generate_kin_input import create_kin_input
from dynamite.data_prep.generate_kin_input import read_atlas3d

dir='Data_prep/Kinematics/'

#for Califa only one input needed, Schwarzschild models in https://arxiv.org/pdf/1709.06649.pdf
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



#------------------------------------

#ATLAS3D, Schwarzschild models in https://arxiv.org/pdf/1902.10175.pdf
#original cube and kinematics cube
input=[dir+'ATLAS3D/MS_NGC4570_r1_C2D.fits',
       dir+'ATLAS3D/NGC4570_4moments_ATLAS3d.fits']

outdir='Data_prep/Dyn_input/ATLAS3D/'
expr=''


#better use
#dyn.model.Model.create_model_directory(output)
if not os.path.exists(outdir):
    os.makedirs(outdir)

read_atlas3d(input)
create_kin_input('NGC4570',input, outdir, expr='', fit_PA=True, kin_input='ATLAS3D', ngh=6)

#add PSF to kinfile
#gh = dyn.kinematics.GaussHermite() 
#gh.add_psf_to_datafile(sigma=[1.06], weight=[1.0], datafile=outdir+'gauss_hermite_kins'+expr+'.ecsv')
