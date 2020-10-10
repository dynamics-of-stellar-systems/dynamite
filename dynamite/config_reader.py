#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:14:17 2020

@author: maindl
"""

# some tricks to add the current path to sys.path (so the imports below work)

import os.path
import sys

this_dir = os.path.dirname(__file__)
if not this_dir in sys.path:
    sys.path.append(this_dir)

# import required modules/packages

import yaml
import physical_system as physys
import parameter_space as parspace
import kinematics as kinem
import populations as popul
import mges as mge
import model
import executor

class Settings(object):
    """
    Class that collects misc configuration settings
    """
    def __init__(self):
        self.orblib_settings = {}
        self.parameter_space_settings = {}

    def add(self, kind, values):
        if kind == 'orblib_settings':
            self.orblib_settings = values
        elif kind == 'parameter_space_settings':
            self.parameter_space_settings = values
        elif kind == 'legacy_settings':
            self.legacy_settings = values
        elif kind == 'io_settings':
            self.io_settings = values
        elif kind == 'weight_solver_settings':
            self.weight_solver_settings = values
        elif kind == 'executor_settings':
            self.executor_settings = values
        else:
            raise ValueError("""Config only takes orblib_settings
                             and parameter_space_settings
                             and legacy settings
                             and io_settings
                             and weight_solver_settings
                             and executor_settings""")

    def validate(self):
        if not(self.orblib_settings and self.parameter_space_settings and
               self.output_settings and self.weight_solver_settings
               and self.executor_settings):
            raise ValueError("""Config needs orblib_settings
                             and parameter_space_settings
                             and io_settings
                             and weight_solver_settings
                             and executor_settings""")

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.__dict__})')


class Configuration(object):
    """
    Reads the configuration file and instantiates the objects
    self.system, ...
    """
    def __init__(self, filename=None, silent=False):
        """
        Reads configuration file and instantiates objects.
        Does some rudimentary checks for consistency.

        Parameters
        ----------
        filename : string, needs to refer to an existing file including path
        silent : True suppresses output, default=False

        Raises
        ------
        FileNotFoundError
            If file does not exist or filename is None or not given.

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
        self.settings = Settings() # instantiate Settings object

        if not silent: # get paths first
            print('io_settings...')
            print(f' {tuple(self.params["io_settings"].keys())}')
        self.settings.add('io_settings', self.params['io_settings'])
        try:
            for io in ['input', 'output']:
                d = self.settings.io_settings[io+'_directory']
                if len(d) > 0 and d[-1] != '/': # len(d)=0: allow no path, too
                    self.settings.io_settings[io+'_directory'] += '/'
        except:
            raise ValueError('io_settings: check input_directory '
                             'and output_directory')
        # print(self.settings.io_settings)
        # print(self.params['io_settings']['input_directory'])
        # print(self.params['io_settings']['output_directory'])

        for key, value in self.params.items(): # walk through file contents...

            # add components to system

            if key == 'system_components':
                if not silent:
                    print('model_components:')
                for comp, data_comp in value.items():
                    if not data_comp['include']:
                            if not silent:
                                print('', comp, '  ...ignored')
                            continue

                    # instantiate the component

                    if not silent:
                        print(f" {comp}... instantiating {data_comp['type']} "
                              "object")
                    if 'contributes_to_potential' not in data_comp:
                        raise ValueError(f'Component {comp} needs '
                                         'contributes_to_potential attribute')
#                    c = globals()[data_comp['type']](contributes_to_potential
#                        = data_comp['contributes_to_potential'])
                    c = getattr(physys,data_comp['type'])(name = comp,
                            contributes_to_potential \
                            = data_comp['contributes_to_potential'])

                    # initialize the component's paramaters, kinematics,
                    # and populations

                    par_list, kin_list, pop_list = [], [], []

                    # read parameters

                    if 'parameters' not in data_comp:
                        raise ValueError('Component ' + comp + \
                                         ' needs parameters')
                    if not silent:
                        print(f" Has parameters "
                              f"{tuple(data_comp['parameters'].keys())}")
                    for par, data_par in data_comp['parameters'].items():
                        the_parameter = parspace.Parameter(name=par,**data_par)
                        par_list.append(the_parameter)
                    c.parameters = par_list

                    # read kinematics

                    if 'kinematics' in data_comp:
                    # shall we include a check here (e.g., only
                    # VisibleComponent has kinematics?)
                        if not silent:
                            print(f" Has kinematics "
                                  f"{tuple(data_comp['kinematics'].keys())}")
                        for kin, data_kin in data_comp['kinematics'].items():
                            path=self.settings.io_settings['input_directory']
                            kinematics_set = getattr(kinem,data_kin['type'])\
                                                (name=kin,
                                                 input_directory=path,
                                                 **data_kin)
                            kin_list.append(kinematics_set)
                        c.kinematic_data = kin_list

                    # read populations

                    if 'populations' in data_comp:
                    # shall we include a check here (e.g., only
                    # VisibleComponent has populations?)
                        if not silent:
                            print(f" Has populations "
                                  f"{tuple(data_comp['populations'].keys())}")
                        for pop, data_pop in data_comp['populations'].items():
                            populations_set = popul.Populations(name=pop,
                                                                **data_pop)
                            pop_list.append(populations_set)
                        c.population_data = pop_list

                    if 'mge_file' in data_comp:
                        path = self.settings.io_settings['input_directory']
                        c.mge = mge.MGE(input_directory=path,
                                        datafile=data_comp['mge_file'])

                    # add component to system
                    c.validate() # now also adds the right parameter sformat
                    self.system.add_component(c)

            # add system parameters

            elif key == 'system_parameters':
                if not silent:
                    print('system_parameters...')
                    print(f' {tuple(value.keys())}')
                par_list = []
                for other, data in value.items():
                    par_list.append(parspace.Parameter(name=other, **data))
                setattr(self.system, 'parameters', par_list)
                    # if other == 'ml':
                    #     self.system.ml=parspace.Parameter(name=other,**data)
                    # else:
                    #     setattr(self.system, other, data)

            # add system attributes

            elif key == 'system_attributes':
                if not silent:
                    print('system_attributes...')
                    print(f' {tuple(value.keys())}')
                for other, data in value.items():
                    setattr(self.system, other, data)

            # add orbit library settings to Settings object

            elif key == 'orblib_settings':
                if not silent:
                    print('orblib_settings...')
                    print(f' {tuple(value.keys())}')
                self.settings.add('orblib_settings', value)

            # add parameter space settings to Settings object

            elif key == 'parameter_space_settings':
                if not silent:
                    print('parameter_space_settings...')
                    print(f' {tuple(value.keys())}')
                self.settings.add('parameter_space_settings', value)

            # add legacy settings to Settings object

            elif key == 'legacy_settings':
                if not silent:
                    print('legacy_settings...')
                    print(f' {tuple(value.keys())}')
                if value['directory'] == 'default':
                    # this_dir is 'dynamite'
                    value['directory'] = this_dir+'/../legacy_fortran'
                # remove trailing / from path if provided
                if value['directory'][-1]=='/':
                    value['directory'] = value['directory'][:-1]
                self.settings.add('legacy_settings', value)

            # add output settings to Settings object

            elif key == 'io_settings':
                pass # io_settings (paths) have been assigned already...
                # if not silent:
                #     print('io_settings...')
                #     print(f' {tuple(value.keys())}')
                # self.settings.add('io_settings', value)

            # add weight_solver_settings to Settings object

            elif key == 'weight_solver_settings':
                if not silent:
                    print('weight_solver_settings...')
                    print(f' {tuple(value.keys())}')
                self.settings.add('weight_solver_settings', value)

            # add executor_settings to Settings object

            elif key == 'executor_settings':
                if not silent:
                    print('executor_settings...')
                    print(f' {tuple(value.keys())}')
                self.settings.add('executor_settings', value)

            else:
                raise ValueError(f'Unknown configuration key: {key}')

        self.system.validate() # now also adds the right parameter sformat
        if not silent:
            print(f'**** System assembled:\n{self.system}')
            print(f'**** Settings:\n{self.settings}')
        self.validate()
        if not silent:
            print('**** Configuration validated')

        self.parspace = parspace.ParameterSpace(self.system)
        if not silent:
            print('**** Instantiated parameter space')
            print(f'**** Parameter space:\n{self.parspace}')

        self.all_models = model.AllModels(parspace=self.parspace,
                                          settings=self.settings)
        if not silent:
            print('**** Instantiated AllModels object:\n'
                  f'{self.all_models.table}')

        kw_executor = {'system':self.system,
                       'legacy_directory':
                           self.settings.legacy_settings['directory'],
                       'executor_settings':self.settings.executor_settings}
        executor_type = self.settings.executor_settings['type']
        self.executor = getattr(executor, executor_type)(**kw_executor)
        if not silent:
            print(f'**** Instantiated executor object: {executor_type}')

    def validate(self):
        """
        Validates the system and settings. This method is still VERY
        rudimentary and will be adjusted as we add new functionality
        to dynamite. Currently, this method is geared towards legacy mode.

        Returns
        -------
        None.

        """
        if sum(1 for i in self.system.cmp_list \
               if isinstance(i, physys.Plummer)) != 1:
            raise ValueError('System needs to have exactly one Plummer object')
        if sum(1 for i in self.system.cmp_list \
               if isinstance(i, physys.VisibleComponent)) != 1:
            raise ValueError('System needs to have exactly one '
                             'VisibleComponent object')
        if sum(1 for i in self.system.cmp_list \
               if issubclass(type(i), physys.DarkComponent)
               and not isinstance(i, physys.Plummer)) > 1:
            raise ValueError('System needs to have zero or one DM Halo object')
        if not ( 1 < len(self.system.cmp_list) < 4):
            raise ValueError('System needs to comprise exactly one Plummer, '
                             'one VisibleComponent, and zero or one DM Halo '
                             'object(s)')

        for c in self.system.cmp_list:
            if issubclass(type(c), physys.VisibleComponent): # Check vis. comp.
                if c.kinematic_data:
                    for kin_data in c.kinematic_data:
                        if kin_data.type != 'GaussHermite':
                            raise ValueError('VisibleComponent kinematics '
                                             'need GaussHermite type')
                else:
                    raise ValueError('VisibleComponent must have kinematics '
                                     'of type GaussHermite')
                if c.symmetry != 'triax':
                    raise ValueError('Legacy mode: VisibleComponent must be '
                                     'triaxial')
                continue
            if issubclass(type(c), physys.DarkComponent) \
                and not isinstance(c, physys.Plummer):
            # Check allowed dm halos in legacy mode
                if type(c) not in [physys.NFW, physys.Hernquist,
                                   physys.TriaxialCoredLogPotential,
                                   physys.GeneralisedNFW]:
                    raise ValueError(f'DM Halo needs to be of type NFW, '
                                     f'Hernquist, TriaxialCoredLogPotential, '
                                     f'or GeneralisedNFW, not {type(c)}')

        gen_type = self.settings.parameter_space_settings["generator_type"]
        if gen_type != 'GridWalk' and gen_type != 'LegacyGridSearch':
            raise ValueError('Legacy mode: parameter space generator_type '
                             'must be GridWalk or LegacyGridSearch')
        if self.settings.parameter_space_settings["which_chi2"] not in \
            ["chi2", "kinchi2"]:
            raise ValueError('Unknown which_chi2 setting, use chi2 or kinchi2')


    # def read_parameters(self, par=None, items=None):
    #     """
    #     Will add each key-value pair in items to parameters object par by calling
    #     par.add(...) and subsequently calls the par.validate() method.

    #     Parameters
    #     ----------
    #     par : ...parameters object, optional
    #         The default is None.
    #     items : dictionary, optional
    #         The default is None.

    #     Returns
    #     -------
    #     None.

    #     """

    #     for p, v in items:
    #         par.add(name=p, **v)
    #     par.validate()
