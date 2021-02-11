import os.path
import sys
import math

import yaml
import logging
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
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
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
            text = """Config only takes orblib_settings
                             and parameter_space_settings
                             and legacy settings
                             and io_settings
                             and weight_solver_settings
                             and multiprocessing_settings"""
            self.logger.error(text)
            raise ValueError(text)

    def validate(self):
        if not(self.orblib_settings and self.parameter_space_settings and
               self.output_settings and self.weight_solver_settings
               and self.multiprocessing_settings):
            text = """Config needs orblib_settings
                             and parameter_space_settings
                             and io_settings
                             and weight_solver_settings
                             and multiprocessing_settings"""
            self.logger.error(text)
            raise ValueError(text)

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.__dict__})')

class UniqueKeyLoader(yaml.SafeLoader):
    """
    Special yaml loader with duplicate key checking.
    Credits: ErichBSchulz,
    https://stackoverflow.com/questions/33490870/parsing-yaml-in-python-detect-duplicated-keys
    """

    def construct_mapping(self, node, deep=False):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        mapping = []
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            if key in mapping:
                text = f'Duplicate key in configuration file: {key}'
                self.logger.error(text)
                raise AssertionError(text)
            # assert key not in mapping, \
            #     f'Duplicate key in configuration file: {key}'
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

    def __init__(self, filename=None, silent=None, reset_logging=False):
        """
        Reads configuration file and instantiates objects.
        Does some rudimentary checks for consistency.

        Parameters
        ----------
        filename : string, needs to refer to an existing file including path
        silent : DEPRECATED (diagnostic output handled by logging module)
        reset_logging : bool, if False: use the calling application's logging
                                        settings
                              if True: set logging to Dynamite defaults

        Raises
        ------
        FileNotFoundError
            If file does not exist or filename is None or not given.

        Returns
        -------
        None.

        """
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        logger = self.logger
        if reset_logging is True:
            logger.info('Resetting logging configuration')
            DynamiteLogging()
            logger.debug('Logging set to Dynamite defaults')
        else:
            logger.debug("Dynamite uses the calling application's logging "
                         "settings")
        if silent is not None:
            self.logger.warning("'silent' option is deprecated and ignored")
        
        legacy_dir = \
            os.path.realpath(os.path.dirname(__file__)+'/../legacy_fortran')
            # os.path.dirname(os.path.realpath(__file__))+'/../'legacy_fortran'
        self.logger.debug(f'Legacy Fortran folder: {legacy_dir}')
        
        try:
            with open(filename, 'r') as f:
                # self.params = yaml.safe_load(f)
                config_text = f.read()
        except:
            logger.error(f'Cannot open {filename}, please specify '
                         'existing file')
            raise
        self.params = yaml.load(config_text, Loader=UniqueKeyLoader)
        logger.info(f'Config file {filename} read.')

        self.system = physys.System() # instantiate System object
        self.settings = Settings() # instantiate Settings object

        # get paths first
        logger.info('io_settings...')
        logger.debug(f'Read: {self.params["io_settings"]}')

        self.settings.add('io_settings', self.params['io_settings'])
        try:
            for io in ['input', 'output']:
                d = self.settings.io_settings[io+'_directory']
                if len(d) > 0 and d[-1] != '/': # len(d)=0: allow no path, too
                    self.settings.io_settings[io+'_directory'] += '/'
        except:
            logger.error('io_settings: check input_directory '
                         'and output_directory')
            raise
        logger.debug('io_settings assigned to Settings object')

        for key, value in self.params.items(): # walk through file contents...

            # add components to system

            if key == 'system_components':
                logger.info('model_components...')
                for comp, data_comp in value.items():
                    if not data_comp['include']:
                            logger.info(f'{comp}... ignored')
                            continue

                    # instantiate the component

                    logger.debug(f"{comp}... instantiating {data_comp['type']} "
                              "object")
                    if 'contributes_to_potential' not in data_comp:
                        text = f'Component {comp} needs ' + \
                                'contributes_to_potential attribute'
                        logger.error(text)
                        raise ValueError(text)
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
                        text = f'Component {comp} needs parameters'
                        logger.error(text)
                        raise ValueError(text)
                    logger.debug('Has parameters '
                                 f'{tuple(data_comp["parameters"].keys())}')
                    for par, data_par in data_comp['parameters'].items():
                        p = par + '_' + comp
                        the_parameter = parspace.Parameter(name=p,**data_par)
                        par_list.append(the_parameter)
                    c.parameters = par_list

                    # read kinematics

                    if 'kinematics' in data_comp:
                    # shall we include a check here (e.g., only
                    # VisibleComponent has kinematics?)
                        logger.debug('Has kinematics '
                                    f'{tuple(data_comp["kinematics"].keys())}')
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
                        logger.debug(f'Has populations '
                                f'{tuple(data_comp["populations"].keys())}')
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
                    parset = {p.name[:p.name.rfind(f'_{c.name}')]:p.value \
                              for p in c.parameters}
                    if not c.validate_parset(parset):
                        raise ValueError(f'{c.name}: invalid parameters '
                                         f'{parset}')
                    self.system.add_component(c)

                # once all components added, put all kinematic_data in a list
                self.system.get_all_kinematic_data()

            # add system parameters

            elif key == 'system_parameters':
                logger.info('system_parameters...')
                logger.debug(f'system_parameters: {tuple(value.keys())}')
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
                logger.info('system_attributes...')
                logger.debug(f'system_attributes: {tuple(value.keys())}')
                for other, data in value.items():
                    setattr(self.system, other, data)

            # add orbit library settings to Settings object

            elif key == 'orblib_settings':
                logger.info('orblib_settings...')
                logger.debug(f'orblib_settings: {tuple(value.keys())}')
                self.settings.add('orblib_settings', value)

            # add parameter space settings to Settings object

            elif key == 'parameter_space_settings':
                logger.info('parameter_space_settings...')
                logger.debug(f'parameter_space_settings: {tuple(value.keys())}')
                self.settings.add('parameter_space_settings', value)

            # add legacy settings to Settings object

            elif key == 'legacy_settings':
                logger.info('legacy_settings...')
                logger.debug(f'legacy_settings: {tuple(value.keys())}')
                if value['directory'] == 'default':
                    value['directory'] = legacy_dir
                # remove trailing / from path if provided
                if value['directory'][-1]=='/':
                    value['directory'] = value['directory'][:-1]
                self.settings.add('legacy_settings', value)

            # add output settings to Settings object

            elif key == 'io_settings':
                pass # io_settings (paths) have been assigned already...

            # add weight_solver_settings to Settings object

            elif key == 'weight_solver_settings':
                logger.info('weight_solver_settings...')
                logger.debug(f'weight_solver_settings: {tuple(value.keys())}')
                self.settings.add('weight_solver_settings', value)

            # add multiprocessing_settings to Settings object

            elif key == 'multiprocessing_settings':
                logger.info('multiprocessing_settings...')
                logger.debug(f'multiprocessing_settings: {tuple(value.keys())}')
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
                    logger.info(f"... using {value['ncpus']} CPUs.")
                self.settings.add('multiprocessing_settings', value)

            else:
                text = f'Unknown configuration key: {key}'
                logger.error(text)
                raise ValueError(text)

        self.system.validate() # now also adds the right parameter sformat
        parset = {p.name:p.value for p in self.system.parameters}
        if not self.system.validate_parset(parset):
            raise ValueError(f'Invalid sysetm parameters {parset}')
        logger.info('System assembled')
        logger.debug(f'System: {self.system}')
        logger.debug(f'Settings: {self.settings}')

        self.validate()
        logger.info('Configuration validated')

        if 'generator_settings' in self.settings.parameter_space_settings:
            self.set_threshold_del_chi2( \
                self.settings.parameter_space_settings['generator_settings'])

        self.parspace = parspace.ParameterSpace(self.system)
        logger.info('Instantiated parameter space')
        logger.debug(f'Parameter space: {self.parspace}')

        self.all_models = model.AllModels(parspace=self.parspace,
                                          settings=self.settings,
                                          system=self.system)
        logger.info('Instantiated AllModels object')
        logger.debug(f'AllModels:\n{self.all_models.table}')

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
            self.logger.error('System must have exactly one Plummer object')
            raise ValueError('System must have exactly one Plummer object')
        if sum(1 for i in self.system.cmp_list \
               if isinstance(i, physys.VisibleComponent)) != 1:
            self.logger.error('System needs to have exactly one '
                              'VisibleComponent object')
            raise ValueError('System needs to have exactly one '
                             'VisibleComponent object')
        if sum(1 for i in self.system.cmp_list \
               if issubclass(type(i), physys.DarkComponent)
               and not isinstance(i, physys.Plummer)) > 1:
            self.logger.error('System must have zero or one DM Halo object')
            raise ValueError('System must have zero or one DM Halo object')
        if not ( 1 < len(self.system.cmp_list) < 4):
            self.logger.error('System needs to comprise exactly one Plummer, '
                              'one VisibleComponent, and zero or one DM Halo '
                              'object(s)')
            raise ValueError('System needs to comprise exactly one Plummer, '
                             'one VisibleComponent, and zero or one DM Halo '
                             'object(s)')

        for c in self.system.cmp_list:
            if issubclass(type(c), physys.VisibleComponent): # Check vis. comp.
                if c.kinematic_data:
                    for kin_data in c.kinematic_data:
                        if kin_data.type != 'GaussHermite':
                            self.logger.error('VisibleComponent kinematics '
                                              'need GaussHermite type')
                            raise ValueError('VisibleComponent kinematics '
                                             'need GaussHermite type')
                else:
                    self.logger.error('VisibleComponent must have kinematics '
                                      'of type GaussHermite')
                    raise ValueError('VisibleComponent must have kinematics '
                                     'of type GaussHermite')
                if c.symmetry != 'triax':
                    self.logger.error('Legacy mode: VisibleComponent must be '
                                      'triaxial')
                    raise ValueError('Legacy mode: VisibleComponent must be '
                                     'triaxial')
                continue
            if issubclass(type(c), physys.DarkComponent) \
                and not isinstance(c, physys.Plummer):
            # Check allowed dm halos in legacy mode
                if type(c) not in [physys.NFW, physys.Hernquist,
                                   physys.TriaxialCoredLogPotential,
                                   physys.GeneralisedNFW]:
                    self.logger.error(f'DM Halo needs to be of type NFW, '
                                      f'Hernquist, TriaxialCoredLogPotential, '
                                      f'or GeneralisedNFW, not {type(c)}')
                    raise ValueError(f'DM Halo needs to be of type NFW, '
                                     f'Hernquist, TriaxialCoredLogPotential, '
                                     f'or GeneralisedNFW, not {type(c)}')

        gen_type = self.settings.parameter_space_settings["generator_type"]
        if gen_type != 'GridWalk' and gen_type != 'LegacyGridSearch':
            self.logger.error('Legacy mode: parameter space generator_type '
                              'must be GridWalk or LegacyGridSearch')
            raise ValueError('Legacy mode: parameter space generator_type '
                             'must be GridWalk or LegacyGridSearch')
        chi2abs = self.__class__.thresh_chi2_abs
        chi2scaled = self.__class__.thresh_chi2_scaled
        gen_set=self.settings.parameter_space_settings.get('generator_settings')
        if gen_set != None and (chi2abs in gen_set and chi2scaled in gen_set):
            self.logger.error(f'Only specify one of {chi2abs}, {chi2scaled}, '
                              'not both')
            raise ValueError(f'Only specify one of {chi2abs}, {chi2scaled}, '
                             'not both')


class DynamiteLogging(object):
    """
    Dynamite logging setup. ONLY use if logging has not been configured
    outside of Dynamite.
    If no arguments are give, the logging setup is as follows:
    (1) log to the console with logging level INFO, messages include the level,
        timestamp, class name, and message text
    (2) create a dynamite.log file with logging level DEBUG, messages include
        the level, timestamp, class name, filename:method:line number, and
        message text
    """
    def __init__(self, logfile='dynamite.log', console_level=logging.INFO,
                                               logfile_level=logging.DEBUG,
                                               console_formatter = None,
                                               logfile_formatter = None):
        """
        Initialize Dynamite logging.

        Parameters
        ----------
        logfile : str, optional
            Name of the logfile, logfile=None will not create a logfile.
            The default is 'dynamite.log'.
        console_level : int, optional
            Logfile logging level. The default is logging.INFO.
        logfile_level : int, optional
            Console logging level. The default is logging.DEBUG.
        console_formatter : str, optional
            Format string for console logging. The default is set in the code.
        logfile_formatter : str, optional
            Format string for logfile logging. The default is set in the code.

        Returns
        -------
        None.

        """
        logger = logging.getLogger()       # create logger
        logger.setLevel(logging.DEBUG)     # set level that's lower that wanted

        ch = logging.StreamHandler()       # create console logging handler
        ch.setLevel(console_level)         # set console logging level
        # create formatter
        if console_formatter is None:
            console_formatter = logging.Formatter( \
                '[%(levelname)s] %(asctime)s - %(name)s - '
                '%(message)s', "%H:%M:%S")
        ch.setFormatter(console_formatter) # add the formatter to the handler
        logger.addHandler(ch)              # add the handler to the logger

        if logfile:
            fh = logging.FileHandler(logfile, mode='w') # create handler
            fh.setLevel(logfile_level)             # set file logging level
            # create formatter
            # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            # formatter = logging.Formatter('[%(levelname)s] %(asctime)s.%(msecs)03d - %(name)s - %(message)s', "%Y-%m-%d %H:%M:%S")
            # formatter = logging.Formatter('[%(levelname)s] %(asctime)s - %(filename)s %(funcName)s:%(lineno)d - %(message)s', "%H:%M:%S")
            if logfile_formatter is None:
                logfile_formatter = logging.Formatter( \
                    '[%(levelname)s] %(asctime)s - %(name)s - '
                    '%(filename)s:%(funcName)s:%(lineno)d - '
                    '%(message)s', "%b %d %H:%M:%S" )
                fh.setFormatter(logfile_formatter) # add formatter to handler
            logger.addHandler(fh)                  # add handler to the logger
