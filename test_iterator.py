# %load_ext autoreload
# %autoreload 2

import numpy as np
import dynamite_src as dyn
import matplotlib.pyplot as plt

'''
model_example/NGC6278/
    config_example.yaml
    dynamite_script.py
    input_data/

pwd
>model_example/NGC6278/
python dynamite_script.py

model_example/NGC6278/
    config_example.yaml
    dynamite_script.py
    input_data/
    output_data/
'''

fname = 'model_example/config_legacy_example.yaml'
c = dyn.config_reader.ConfigurationReaderYaml(fname, silent=True)
parspace = dyn.parameter_space.ParameterSpace(c.system)
all_models = dyn.schwarzschild.AllModels(from_file=False,
                                         parspace=parspace,
                                         settings=c.settings)

smi = dyn.shwarzschild_model_iterator.SchwarzschildModelIterator(
    system=c.system,
    all_models=all_models,
    settings=c.settings)

all_models.table['kinchi2']
all_models.table['chi2']

chi2 = all_models.table['chi2']+all_models.table['kinchi2']
min_chi2 = np.min(chi2)
del_chi2 = chi2 - min_chi2


np.where(np.array(parspace.par_names)=='f')

parspace.par_names[2]
parspace.par_names[7]

all_models.table['dc']
all_models.table['ml']


# plot
kwscat = dict(cmap=plt.cm.viridis_r,
              s=100)
plt.scatter(all_models.table['dc'],
            all_models.table['ml'],
            c=del_chi2,
            **kwscat)
kwline = dict(color='k', lw=0.5)
plt.axvline(parspace[2].grid_parspace_settings['lo'],
            **kwline)
plt.axvline(parspace[2].grid_parspace_settings['hi'],
            **kwline)
plt.axhline(parspace[7].grid_parspace_settings['lo'],
            **kwline)
plt.axhline(parspace[7].grid_parspace_settings['hi'],
            **kwline)
plt.colorbar()
plt.show()














# end
