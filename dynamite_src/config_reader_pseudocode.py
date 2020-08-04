import yaml
import dynamite_src.physical_system as physys
import dynamite_src.parameter_space as parspace
import dynamite_src.kinematics as kinematics

class ConfigurationReaderYaml(object):
    """
    Reads the configuration file
    """
    def __init__(self, filename=None): # instantiate the objects here. instead of the dict, self.system will be a System object
        if filename is not None:
            with open(filename, 'r') as f:
                self.params = yaml.safe_load(f)
        else:
            raise FileNotFoundError('Please specify filename')
        for key, values in self.params.items():
            if key == 'model_components':
                print('model_components:')
                cmp_list = []
                for comp, data_comp in values.items():
                    par_list, kin_list, pop_list = [], [], []
                    for par, data_par in data_comp['parameters']:
                        par = parspace.Parameters(...)
                        par_list += [par]
                    for kin, data_kin in data_comp['kinematics']:
                        # instantiate a kinematics object of type data_kin['type']
                        kin = kinematics.data_kin['type']
                        kin_list += [kin]
                    for pop, data_pop['populations']:
                        # read in populations
                        # instantiate a popuulation object of type data_pop['type']
                        # pop = ...
                        pop_list += [pop]
                    # instantiate a component object of type data['type']
                    component = physys.data['type'](
                        contributes_to_potential=data_comp['contributes_to_potential'],
                        kinematics = kin_list,
                        populations = pop_list,
                        parameters = par_list)
                    cmp_list += [component]

        self.system = physys.System(cmp_list)
