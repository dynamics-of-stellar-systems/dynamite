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
n_chi2 = 5
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
t['chi2_min'] = np.amin(chi2_all, axis=1)
t['chi2_max'] = np.amax(chi2_all, axis=1)
t['chi2_average'] = chi2_average
t['chi2_sdev'] = chi2_sdev

print(t)

# # backup data from a run with n_chi2 = 5
# t = table.Table()
# t['model_id'] = [i for i in range(6)]
# t['chi2_min'] = [265816.2253, 255570.845, 255249.1577, 256636.7563,
#                  265263.2285, 263122.2601]
# t['chi2_max'] = [266004.54280000005, 256590.6705, 255825.74120000002,
#                  256946.507, 266173.28040000005, 263771.3785]
# t['chi2_average'] = [265922.4676, 255948.46872, 255527.1702, 256833.40252,
#                      265673.03636, 263516.12922]
# t['chi2_sdev'] = [84.22188926028149, 362.50408057701304, 233.47116602014836,
#                   113.23921805508898, 312.363540790302, 234.2891485666512]
# print(t)

t.write('tests/chi2_compare.dat', format='ascii')
