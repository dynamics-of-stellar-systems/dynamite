%load_ext autoreload
%autoreload 2

import numpy as np
import dynamite_src as dyn
import matplotlib.pyplot as plt
import os
from astropy.io import ascii

fname = 'pj_model/config_legacy_example.yaml'
c = dyn.config_reader.ConfigurationReaderYaml(fname, silent=True)
parspace = dyn.parameter_space.ParameterSpace(c.system)

# read all existing models, or make empty all_models object if none
all_models = dyn.schwarzschild.AllModels(from_file=True,
                                         parspace=parspace,
                                         settings=c.settings)

all_models.parspace.par_names

parset0 = all_models.get_parset_from_row(0)

mod = dyn.schwarzschild.LegacySchwarzschildModel(
    system=c.system,
    settings=c.settings,
    parset=parset0,
    parspace=parspace,
    execute_run=False
    )
ws = dyn.weight_solvers.LegacyWeightSolver(
    system=c.system,
    mod_dir=mod.directory_noml,
    ml=5.0
    )

mod.directory

chi2, chi2vec, kinchi2 = ws.read_output()

chi2vec

kinchi2

np.sum(chi2vec)

chi2

chi2 + 3.599991E+03

kinchi2


c.system.cmp_list[0].parameters[0].name










# end
