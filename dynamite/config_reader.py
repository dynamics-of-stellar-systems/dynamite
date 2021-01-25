
# some tricks to add the current path to sys.path (so the imports below work)

import os.path
import sys
import math

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
        elif kind == 'multiprocessing_settings':
            self.multiprocessing_settings = values
        else:
            raise ValueError("""Config only takes orblib_settings
                             and parameter_space_settings
                             and legacy settings
                             and io_settings
                             and weight_solver_settings
                             and multiprocessing_settings""")

    def validate(self):
        if not(self.orblib_settings and self.parameter_space_settings and
               self.output_settings and self.weight_solver_settings
               and self.multiprocessing_settings):
            raise ValueError("""Config needs orblib_settings
                             and parameter_space_settings
                             and io_settings
                             and weight_solver_settings
                             and multiprocessing_settings""")

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.__dict__})')

class UniqueKeyLoader(yaml.SafeLoader):
    """
    Special yaml loader with duplicate key checking.
    Credits: ErichBSchulz,
    https://stackoverflow.com/questions/33490870/parsing-yaml-in-python-detect-duplicated-keys
    """
    def construct_mapping(self, node, deep=False):
        mapping = []
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            assert key not in mapping, \
                f'Duplicate key in configuration file: {key}'
            mapping.append(key)
        return super().construct_mapping(node, deep)

class Configuration(object):
    """
    Reads the configuration file and instantiates the objects
    self.system, ...
    """

    # Class attributes
    # threshold_delta_chi2 variants for LegacyGridSearch parameter generator
    thresh_chi2_abs = 'threshold_del_chi2_abs'
    thresh_chi2_scaled = 'threshold_del_chi2_as_frac_of_sqrt2nobs'

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
                # self.params = yaml.safe_load(f)
                config_text = f.read()
                self.params = yaml.load(config_text, Loader=UniqueKeyLoader)
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
                        p = par + '_' + comp
                        the_parameter = parspace.Parameter(name=p,**data_par)
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

            # add multiprocessing_settings to Settings object

            elif key == 'multiprocessing_settings':
                if not silent:
                    print('multiprocessing_settings...')
                    print(f' {tuple(value.keys())}')
                # if submitted as slurm script we must add cwd to path
                try: # check if Slurm being using
                    os.environ["SLURM_JOB_CPUS_PER_NODE"]
                    sys.path.append(os.getcwd())
                except KeyError:
                    pass
                if value['ncpus']=='all_available':
                    try:
                        ncpus = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
                    except KeyError:
                        import multiprocessing
                        ncpus = multiprocessing.cpu_count()
                    value['ncpus'] = ncpus
                if not silent:
                    print(f"... using {value['ncpus']} CPUs.")
                self.settings.add('multiprocessing_settings', value)

            else:
                raise ValueError(f'Unknown configuration key: {key}')

        self.system.validate() # now also adds the right parameter sformat
        if not silent:
            print(f'**** System assembled:\n{self.system}')
            print(f'**** Settings:\n{self.settings}')

        self.validate()
        if not silent:
            print('**** Configuration validated')

        if 'generator_settings' in self.settings.parameter_space_settings:
            self.set_threshold_del_chi2( \
                self.settings.parameter_space_settings['generator_settings'])

        self.parspace = parspace.ParameterSpace(self.system)
        if not silent:
            print('**** Instantiated parameter space')
            print(f'**** Parameter space:\n{self.parspace}')

        self.all_models = model.AllModels(parspace=self.parspace,
                                          settings=self.settings)
        if not silent:
            print('**** Instantiated AllModels object:\n'
                  f'{self.all_models.table}')

    def set_threshold_del_chi2(self, generator_settings):
        """
        Sets threshold_del_chi2 depending on scaled or unscaled input. Works
        with the legacy setup only ('stars' component with one set of
        kinematics).

        Parameters
        ----------
        generator_settings : generator_settings dict

        Returns
        -------
        None.

        """
        chi2abs = self.__class__.thresh_chi2_abs
        chi2scaled = self.__class__.thresh_chi2_scaled
        if chi2abs in generator_settings or chi2scaled in generator_settings:
            two_n_obs = self.get_2n_obs()
            if chi2abs in generator_settings:
                thresh = generator_settings[chi2abs]
            else:
                thresh = generator_settings[chi2scaled] * math.sqrt(two_n_obs)
            generator_settings['threshold_del_chi2'] = thresh

    def get_2n_obs(self):
        """
        Returns 2*n_obs = number_GH * number_spatial_bins. Works with the
        legacy setup only ('stars' component with one set of kinematics).

        Returns
        -------
        two_n_obs : 2*n_obs, int

        """
        number_GH = self.settings.weight_solver_settings['number_GH']
        stars = self.system.get_component_from_name('stars')
        two_n_obs = 2 * number_GH * len(stars.kinematic_data[0].data)
        return two_n_obs

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
        chi2abs = self.__class__.thresh_chi2_abs
        chi2scaled = self.__class__.thresh_chi2_scaled
        gen_set=self.settings.parameter_space_settings.get('generator_settings')
        if gen_set != None and (chi2abs in gen_set and chi2scaled in gen_set):
            raise ValueError(f'Only specify one of {chi2abs}, {chi2scaled}, '
                             'not both')


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
