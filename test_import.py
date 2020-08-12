%load_ext autoreload
%autoreload 2

import dynamite_src.config_reader as cr

c = cr.ConfigurationReaderYaml('./datafiles/config_example.yaml')

# unpackage the system
sys = c.system
print('System has:')
print(f'   - {sys.n_cmp} componenets')
print(f'   (of which {sys.n_pot} contribute to the potential)')
print(f'   - {sys.n_kin} kinematic datasets')
print(f'   - {sys.n_pop} population datasets')
print('\n')
print('The components are:')
for i in range(sys.n_cmp):
    cmp = sys.cmp_list[i]
    print(f'{i}) {cmp.name}')
    print(f'   with {len(cmp.parameters)} parameters')
    if hasattr(cmp, 'mge'):
        print(f'   with mgefile {cmp.mge.datafile} containing data:')
        print(cmp.mge.data)
    print('\n')












# end
