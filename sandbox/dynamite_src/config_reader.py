#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:14:17 2020

@author: maindl
"""

import yaml
import dynamite_src.physical_system as physys
import dynamite_src.parameter_space as parspace

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
        self.system = physys.System()
#        self.__dict__ = par
        for key, values in self.params.items():
            if key == 'model_components':
                print('model_components:')
                for comp, data in values.items():
                    if not data['include']:
                            print('', comp, '  ...ignored')
                            continue
                    if comp == 'black_hole':
                        if data['type'] != 'Plummer':
                            raise ValueError('Black hole must be of type Plummer')
                        print('gas = physys.VisibleComponent(data)', data['contributes_to_potential'])
                        bh = physys.DarkComponent(contributes_to_potential=data['contributes_to_potential'])
                        par = parspace.BlackHoleParameters()
#                            print(par.mass.desc)
                        # instantiate Plummer object here, something like bh = physys.Plummer(data)
                        # add to self.system
                    elif comp == 'gas':
                        if data['type'] != 'VisibleComponent':
                            raise ValueError('Gas must be of type VisibleComponent')
                        else: # code to come
                            print('gas = physys.VisibleComponent(data)', data)
                        par = parspace.GasParameters()
                        # instantiate VisibleComponent object here
                        # add to self.system
                    elif comp == 'stars':
                        if data['type'] != 'VisibleComponent':
                            raise ValueError('Stars must be of type VisibleComponent')
                        else: # code to come
                            print('stars = physys.VisibleComponent(data)', data)
                        par = parspace.StellarParameters()
                        # instantiate VisibleComponent object here
                        # add to self.system
                    # elif... add other components...
                    else:
                        raise ValueError('Unknown component:', comp)

                    for p,v in data['parameters'].items():
                        par.add(name=p,**v)

                    

                    print('', comp, ' ...read')
                    
            elif key == 'other_parameters':
                print('other_parameters exist')
                # add other parameters to system
            elif key == 'orblib_settings':
                print('orblib_settings exist')
                # add orblib setings to config object (for now)
            elif key == 'parameter_space_settings':
                print('parameter_space_settings exist')
                # add perameter space settings to config object (for now)

            else:
                raise ValueError('Unknown configuration key:', key)