#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# first, make sure the paths are set up
# we assume that this script is located and run in the folder dynamite/tests

import os
import sys
import shutil
import numpy as np
from astropy import table
import user_test

#os.chdir('..')
this_dir = os.getcwd()
#print(f'working directory: {this_dir}')
if not this_dir in sys.path:
    sys.path.append(this_dir)

models_folder = this_dir+'/'+'tests/NGC6278/models'
models_file = this_dir+'/'+'tests/NGC6278/all_models.ecsv'

# set n_chi2 to the sample size
n_chi2 = 2
chi2_all = []

# run user_test n_chi2 times
for i in range(n_chi2):
    shutil.rmtree(models_folder, ignore_errors=True)
    if os.path.isfile(models_file):
        os.remove(models_file)
    chi2_all.append([])
    chi2_values = list(user_test.run_user_test()['chi2'])
    for j in range(len(chi2_values)):
        chi2_all[i].append(chi2_values[j])

chi2_all = np.array(chi2_all).T

print(f'chi2 values:\n{chi2_all}')

chi2_average = np.average(chi2_all, axis=1)
chi2_sdev = np.std(chi2_all, axis=1)

t = table.Table()
t['model_id'] = [i for i in range(len(chi2_all))]
t['chi2_average'] = chi2_average
t['chi2_sdev'] = chi2_sdev

print(t)