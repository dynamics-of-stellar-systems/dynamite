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

os.chdir('..')
this_dir = os.getcwd()
print(f'working directory: {this_dir}')
if not this_dir in sys.path:
    sys.path.append(this_dir)

models_folder = this_dir+'/tests/NGC6278/models'
models_file = this_dir+'/tests/NGC6278/all_models.ecsv'

# set n_chi2 to the sample size
n_chi2 = 10
chi2_all = []

# run user_test n_chi2 times
for i in range(n_chi2):
    model_results, output_file, n_max = user_test.run_user_test(stat_mode=True)
    chi2_all.append([])
    chi2_values = list(model_results['chi2'])
    for j in range(len(chi2_values)):
        chi2_all[i].append(chi2_values[j])
    if j<n_max:
        chi2_all[i].append(float("NaN"))

chi2_all = np.array(chi2_all).T

t = table.Table()
t['model_id'] = [i for i in range(len(chi2_all))]
t['chi2_min'] = np.nanmin(chi2_all, axis=1)
t['chi2_max'] = np.nanmax(chi2_all, axis=1)
t['chi2_average'] = np.nanmean(chi2_all, axis=1)
t['chi2_sdev'] = np.nanstd(chi2_all, axis=1)

print(t)

np.savetxt(output_file+f'_raw_{n_chi2}.dat', chi2_all)

t.write(output_file+f'_{n_chi2}.dat', format='ascii')
