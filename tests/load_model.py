# coding: utf-8
import dynamite as dyn
print('Using DYNAMITE version:', dyn.__version__)
print('Located at:', dyn.__path__)

# read configuration
fname = 'test_reimplement_nnls.yaml'
c = dyn.config_reader.Configuration(fname, silent=True)
parset0 = c.all_models.get_parset_from_row(0)
mod0 = dyn.model.LegacySchwarzschildModel(system=c.system, settings=c.settings, parspace=c.parspace, parset=parset0)
orblib0 = mod0.get_orblib()
ws0 = mod0.get_weights(orblib0)
