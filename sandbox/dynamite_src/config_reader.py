#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:14:17 2020

@author: maindl
"""

import yaml
import dynamite_src.physical_system as physys
import dynamite_src.parameter_space as parspace

class Configuration(object):
    """
    Class that collect misc configuration settings
    """
    def __init__(self):
        self.orblib_settings = {}
        self.parameter_space_settings = {}

    def add(self, kind, values):
        if kind == 'orblib_settings':
            self.orblib_settings = values
        elif kind == 'parameter_space_settings':
            self.parameter_space_settings = values
        else:
            raise ValueError('Config only takes orblib_settings and parameter_space_settings')

    def validate(self):
        if not(self.orblib_settings and self.parameter_space_settings):
            raise ValueError('Config needs orblib_settings and parameter_space_settings')


class ConfigurationReaderYaml(object):
    """
    Reads the configuration file and instantiates the objects
    self.system, ...
    """
    def __init__(self, filename=None): # instantiate the objects here. instead of the dict, self.system will be a System object, etc.
        """
        Reads onfiguration file and instantiates objects. Does some rudimentary
        checks for consistency.

        Parameters
        ----------
        filename : string, needs to refer to an existing file including path

        Raises
        ------
        FileNotFoundError
            DESCRIPTION. If file does not exist or filename is None or not given.

        Returns
        -------
        None.

        """

        if filename is not None:
            with open(filename, 'r') as f:
                self.params = yaml.safe_load(f)
        else:
            raise FileNotFoundError('Please specify filename')

        self.system = physys.System() # instantiate System object
        self.config = Configuration() # instantiate Configuration object
#        self.__dict__ = par

        for key, value in self.params.items(): # walk through the file contents...

            # add components to system

            if key == 'model_components':
                print('model_components:')
                for comp, data in value.items():
                    if not data['include']:
                            print('', comp, '  ...ignored')
                            continue

                    # read black hole data and instantiate Plummer object

                    if comp == 'black_hole':
                        if data['type'] != 'Plummer':
                            raise ValueError('Black hole must be of type Plummer')
                        c = physys.Plummer(contributes_to_potential=data['contributes_to_potential'])
                        c.parameters = parspace.BlackHoleParameters()
                        self.read_parameters(c.parameters, data['parameters'].items())

                    # read gas data and instantiate VisibleComponent object, not implemented yet...

                    elif comp == 'gas': # placeholder, not implemented yet...
                        if data['type'] != 'VisibleComponent':
                            raise ValueError('Gas must be of type VisibleComponent')
                        else: # code to come
                            print('gas = physys.VisibleComponent(data)', data)
                        par = parspace.GasParameters()
                        # instantiate VisibleComponent object here
                        # add to self.system

                    # read stars data and instantiate VisibleComponent object

                    elif comp == 'stars':
                        if data['type'] != 'VisibleComponent':
                            raise ValueError('Stars must be of type VisibleComponent')
                        # assume triaxial symmetry
                        c = physys.VisibleComponent(mge_file = data['mge_file'], symmetry='triax', contributes_to_potential=data['contributes_to_potential'])
                        c.parameters = parspace.StellarParameters()
                        self.read_parameters(c.parameters, data['parameters'].items())

                    # read dark halo data and instantiate NFW object

                    elif comp == 'dark_halo':
                        if data['type'] != 'NFW':
                            raise ValueError('Dark halo must be of type NFW')
                        c = physys.NFW(contributes_to_potential=data['contributes_to_potential'])
                        c.parameters = parspace.DarkHaloParameters()
                        self.read_parameters(c.parameters, data['parameters'].items())

                    # elif... add other components...

                    else:
                        raise ValueError('Unknown component: ' + comp)

                    # fill in common attributes

                    c.contributes_to_potential = data['contributes_to_potential']
                    if 'kinematics' in data:
                        # print('Kinematics!')
                        # print(data['kinematics'])
                        k = []
                        for i, kin in data['kinematics'].items():
                            k.append(kin)
                        # print(k)
                        c.kinematic_data = k
                    if 'population' in data:
                        k = []
                        for i, kin in data['population'].items():
                            k.append(kin)
                        c.population_data = k

                    # add component to system

                    self.system.add_component(c)

                    print('', comp, ' ...read')

            # add other parameters to system

            elif key == 'other_parameters':
                print('other_parameters...')
                for other, data in value.items():
                    if other == 'ml':
                        self.system.ml = parspace.Parameter(**data)
                    elif other == 'distMPc':
                        self.system.distMPc = data
                    elif other == 'galname':
                        self.system.galname = data
                    elif other == 'position_angle':
                        self.system.position_angle = data

            # add orbit library settings to config object

            elif key == 'orblib_settings':
                print('orblib_settings...')
                self.config.add('orblib_settings', value)

            # add parameter space settings to config object

            elif key == 'parameter_space_settings':
                print('parameter_space_settings...')
                self.config.add('parameter_space_settings', value)

            else:
                raise ValueError('Unknown configuration key: ' + key)

        self.system.validate()
        self.config.validate()


    def read_parameters(self, par=None, items=None):
        """
        Will add each key-value pair in items to parameters object par by calling
        par.add(...) and subsequently calls the par.validate() method.

        Parameters
        ----------
        par : ...parameters object, optional
            The default is None.
        items : dictionary, optional
            The default is None.

        Returns
        -------
        None.

        """
        
        for p, v in items:
            par.add(name=p, **v)
        par.validate()
