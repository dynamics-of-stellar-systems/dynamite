import os
import sys
import shutil
import glob
import fnmatch
import math
import logging
import importlib
import yaml

import dynamite as dyn
from dynamite import physical_system as physys
from dynamite import parameter_space as parspace
from dynamite import kinematics as kinem
from dynamite import populations as popul
from dynamite import mges as mge
from dynamite import model

class Settings(object):
    """
    Class to hold all configuration settings

    Has a dictionary attribute for each entry of the second
    section of the YAML config file (e.g. orblib_settings etc...)

    """
    def __init__(self):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.orblib_settings = {}
        self.parameter_space_settings = {}
        self.legacy_settings = {}
        self.io_settings = {}
        self.weight_solver_settings = {}
        self.multiprocessing_settings = {}

    def add(self, kind, values):
        """Add each setting to the object
        """
        if kind == 'orblib_settings':
            self.orblib_settings = values
        elif kind == 'parameter_space_settings':
            self.parameter_space_settings = values
        elif kind == 'legacy_settings':
            self.legacy_settings = values
        elif kind == 'io_settings':
            try:
                out_dir = values['output_directory']
            except KeyError as e:
                text = 'Output directory not set in config file.'
                self.logger.error(text)
                raise Exception(text) from e
            self.io_settings = values
            self.io_settings['model_directory'] = out_dir + 'models/'
            self.io_settings['plot_directory'] = out_dir + 'plots/'
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
        """Validate that all expected settings are present
        """
        if not(self.orblib_settings and self.parameter_space_settings and
               self.io_settings and self.weight_solver_settings
               and self.multiprocessing_settings):
            text = """Config needs orblib_settings
                             and parameter_space_settings
                             and io_settings
                             and weight_solver_settings
                             and multiprocessing_settings"""
            self.logger.error(text)
            raise ValueError(text)
        self.logger.debug('Settings validated.')

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__dict__})'

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
    Reads configuration file

    Does some rudimentary checks for consistency.
    Builds the output directory tree if it does not exist already
    (does not delete existing data).

    Parameters
    ----------
    filename : string
        needs to refer to an existing file including path
    silent : DEPRECATED
        (diagnostic output handled by logging module)
    reset_logging : bool
        if False: use the calling application's logging settings
        if True: set logging to Dynamite defaults
    reset_existing_output : bool
        if False: do not touch existing data in the output directory tree
        if True: rebuild the output directory tree and delete existing data

    Raises
    ------
    FileNotFoundError
        If file does not exist or filename is None or not given.

    Returns
    -------
    sets attributes:
        - ``self.system``: a ``dyn.physical_system.System`` object
        - ``self.cmp_list``: a list of ``dyn.physical_system.Component`` objects
        - ``self.settings``: a list of ``dyn.config_reader.Settings`` object
        - ``self.settings``: a list of ``dyn.config_reader.Settings`` object

    """

    # Class attributes
    # threshold_delta_chi2 variants for LegacyGridSearch parameter generator
    thresh_chi2_abs = 'threshold_del_chi2_abs'
    thresh_chi2_scaled = 'threshold_del_chi2_as_frac_of_sqrt2nobs'

    def __init__(self, filename=None, silent=None, reset_logging=False,
                 reset_existing_output=False):
        if reset_logging is True:
            DynamiteLogging()
            self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
            self.logger.debug('Logging reset to Dynamite defaults')
        else:
            self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
            self.logger.debug("Dynamite uses the calling application's "
                              "logging settings")
        logger = self.logger
        self.logger.debug(f'This is Python {sys.version.split()[0]}')
        self.logger.debug(f'Using DYNAMITE version {dyn.__version__} '
                          f'located at {dyn.__path__}')

        if silent is not None:
            self.logger.warning("'silent' option is deprecated and ignored")

        legacy_dir = \
            os.path.realpath(os.path.dirname(__file__)+'/../legacy_fortran')
            # os.path.dirname(os.path.realpath(__file__))+'/../'legacy_fortran'
        self.logger.debug(f'Default legacy Fortran directory: {legacy_dir}.')

        self.config_file_name = filename
        try:
            with open(self.config_file_name, 'r') as f:
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

        try:
            for io in ['input', 'output']:
                d = self.params['io_settings'][io+'_directory']
                if len(d) > 0 and d[-1] != '/': # len(d)=0: allow no path, too
                    self.params['io_settings'][io+'_directory'] += '/'
        except:
            logger.error('io_settings: check input_directory '
                         'and output_directory in config file')
            raise
        self.settings.add('io_settings', self.params['io_settings'])
        logger.debug('io_settings assigned to Settings object')
        out_dir = self.settings.io_settings['output_directory']
        if reset_existing_output:
            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir)
                self.logger.info(f'Output directory tree {out_dir} removed.')
        self.make_output_directory_tree()
        self.logger.info(f'Output directory tree: {out_dir}.')

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
                        p = f'{par}-{comp}'
                        the_parameter = parspace.Parameter(name=p,**data_par)
                        par_list.append(the_parameter)
                    c.parameters = par_list

                    # read kinematics

                    if 'kinematics' in data_comp:
                    # shall we include a check here (e.g., only
                    # VisibleComponent has kinematics?)
                        logger.debug('Has kinematics '
                                    f'{list(data_comp["kinematics"].keys())}')
                        for kin, data_kin in data_comp['kinematics'].items():
                            path=self.settings.io_settings['input_directory']
                            kinematics_set = getattr(kinem,data_kin['type'])\
                                                (name=kin,
                                                 input_directory=path,
                                                 **data_kin)
                            kin_list.append(kinematics_set)
                        c.kinematic_data = kin_list

                    # cast hist. values to correct numeric type unless `default`
                    for i, k in enumerate(c.kinematic_data):
                        if (k.hist_width=='default') is False:
                            logger.debug(f'hist_width = {k.hist_width}')
                            k.hist_width = float(k.hist_width)
                        if (k.hist_center=='default') is False:
                            logger.debug(f'hist_center = {k.hist_center}')
                            k.hist_center = float(k.hist_center)
                        if (k.hist_bins=='default') is False:
                            logger.debug(f'hist_bins = {k.hist_bins}')
                            k.hist_bins = int(k.hist_bins)

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

                    if 'mge_pot' in data_comp:
                        path = self.settings.io_settings['input_directory']
                        c.mge_pot = mge.MGE(input_directory=path,
                                        datafile=data_comp['mge_pot'])
                    if 'mge_lum' in data_comp:
                        path = self.settings.io_settings['input_directory']
                        c.mge_lum = mge.MGE(input_directory=path,
                                        datafile=data_comp['mge_lum'])

                    # add component to system
                    c.validate()
                    parset = {c.get_parname(p.name):p.raw_value \
                              for p in c.parameters}
                    if not c.validate_parset(parset):
                        text = f'{c.name}: invalid parameters {parset}'
                        self.logger.error(text)
                        if type(c) is physys.TriaxialVisibleComponent:
                            text = c.suggest_parameter_values()
                            self.logger.error(text)
                        raise ValueError(text)
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
                # set a default value to
                #      orblib_settings --> use_new_mirroring : False
                if 'use_new_mirroring' in value.keys():
                    pass
                else:
                    value.update({'use_new_mirroring':True})
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
                self.logger.debug("Legacy directory set to "
                                  f"{value['directory']}.")

            # add output settings to Settings object

            elif key == 'io_settings':
                pass # io_settings (paths) have been assigned already...

            # add weight_solver_settings to Settings object

            elif key == 'weight_solver_settings':
                logger.info('weight_solver_settings...')
                if 'reattempt_failures' not in value:
                    value['reattempt_failures'] = True
                if value['reattempt_failures']:
                    logger.info('Will attempt to recover partially run models.')
                logger.debug(f'weight_solver_settings: {tuple(value.keys())}')
                self.settings.add('weight_solver_settings', value)

            # add multiprocessing_settings to Settings object

            elif key == 'multiprocessing_settings':
                logger.info('multiprocessing_settings...')
                # if submitted as slurm script we must add cwd to path
                try: # check if Slurm being using
                    os.environ["SLURM_JOB_CPUS_PER_NODE"]
                    sys.path.append(os.getcwd())
                except KeyError:
                    pass
                if 'ncpus' not in value:
                    value['ncpus'] = 'all_available'
                if value['ncpus']=='all_available':
                    value['ncpus'] = self.get_n_cpus()
                logger.info(f"... using {value['ncpus']} CPUs "
                             "for orbit integration.")
                if 'ncpus_weights' not in value:
                    value['ncpus_weights'] = value['ncpus']
                elif value['ncpus_weights'] == 'all_available':
                    value['ncpus_weights'] = self.get_n_cpus()
                logger.info(f"... using {value['ncpus_weights']} CPUs "
                            "for weight solving.")
                if 'modeliterator' not in value:
                    value['modeliterator'] = 'ModelInnerIterator'
                logger.debug(f"... using iterator {value['modeliterator']}.")
                if 'orblibs_in_parallel' not in value:
                    value['orblibs_in_parallel'] = True
                logger.debug("... integrate orblibs in parallel: "
                             f"{value['orblibs_in_parallel']}.")
                logger.debug(f'multiprocessing_settings: {tuple(value.keys())}')
                self.settings.add('multiprocessing_settings', value)

            else:
                text = f'Unknown configuration key: {key}'
                logger.error(text)
                raise ValueError(text)

        self.system.validate() # now also adds the right parameter sformat
        parset = {p.name:p.raw_value for p in self.system.parameters}
        if not self.system.validate_parset(parset):
            text = f'Invalid system parameters {parset}'
            self.logger.error(text)
            raise ValueError(text)
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

        self.all_models = model.AllModels(config=self)
        logger.info('Instantiated AllModels object')
        logger.debug(f'AllModels:\n{self.all_models.table}')

        # self.backup_config_file(reset=False)

    def get_n_cpus(self):
        """"
        Returns the number of avalable CPUs.
        """
        try:
            ncpus = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
        except KeyError:
            import multiprocessing
            ncpus = multiprocessing.cpu_count()
        return ncpus

    def set_threshold_del_chi2(self, generator_settings):
        """
        Sets threshold_del_chi2 depending on scaled or unscaled input.

        Works with the legacy setup only (stars component of class
        TriaxialVisibleComponent with one or more sets of kinematics).

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
        Get 2 * number of kinematic observations

        Used for scaling threshold chi2 values for parameter searches. For
        kinemtic type:

            - ``GaussHermite``, then n_obs = number_GH * number_spatial_bins
            - ``BayesLOSVD``, then n_obs = n_LOSVD_bins * number_spatial_bins

        This returns the sum of (2 * n_obs) for all kinematic sets. Take
        kinematics from the ``TriaxialVisibleComponent`` of the system.

        Returns
        -------
        int
            2 * total number of kinematic observations

        """
        stars = \
          self.system.get_component_from_class(physys.TriaxialVisibleComponent)
        n_obs = 0.
        for k in stars.kinematic_data:
            if k.type == 'GaussHermite':
                number_GH = self.settings.weight_solver_settings['number_GH']
                n_obs += number_GH * len(k.data)
            if k.type == 'BayesLOSVD':
                nvbins = k.data.meta['nvbins']
                n_obs += nvbins * len(k.data)
        two_n_obs = 2 * n_obs
        return two_n_obs

    def remove_existing_orblibs(self):
        """
        Removes the entire model output tree, including all existing orblibs

        Returns
        -------
        None.

        """
        model_dir = self.settings.io_settings['model_directory']
        if os.path.isdir(model_dir):
            shutil.rmtree(model_dir)
            self.logger.info(f'Model output tree {model_dir} removed.')
        else:
            self.logger.warning(f'No model output at {model_dir} to remove.')

    def remove_existing_orbital_weights(self):
        """
        Removes existing orbital weights ('ml' directories).

        Deletes all files matching ``output/*/ml*/``

        Raises
        ------
        Exception if directories cannot be removed.

        Returns
        -------
        None.

        """
        ml_pattern = self.settings.io_settings['model_directory'] + '*/ml*'
        ml_directories = glob.glob(ml_pattern)
        if len(ml_directories) > 0:
            for directory in ml_directories:
                shutil.rmtree(directory)
                self.logger.debug(f'Directory {directory} removed.')
            self.logger.info(f'Orbital weights {ml_pattern} removed.')
        else:
            self.logger.info(f'No orbital weights {ml_pattern} to remove.')

    def remove_existing_plots(self, remove_directory=False):
        """
        Removes existing plots from the plots directory.

        Optionally, the plot directory tree can be removed recursively.

        Parameters
        ----------
        remove_directory : BOOL, optional
            True if the plot directory shall be removed, too. If False,
            only regular files in the plot directory are deleted
            (subdirectories will remain untouched in that case).
            The default is False.

        Raises
        ------
        Exception if files/directories cannot be removed.

        Returns
        -------
        None.

        """
        plot_dir = self.settings.io_settings['plot_directory']
        if os.path.isdir(plot_dir):
            if remove_directory:
                shutil.rmtree(plot_dir)
                self.logger.info(f'Plot directory {plot_dir} deleted.')
            else:
                plot_files = glob.glob(plot_dir + '*')
                for f in plot_files:
                    if os.path.isfile(f):
                        os.remove(f)
                self.logger.info(f'Removed files in {plot_dir}.')
        else:
            self.logger.warning(f'Directory {plot_dir} not found, cannot '
                                'remove plots.')

    def remove_existing_all_models_file(self, wipe_other_files=False):
        """
        Deletes the all models file

        Deletes the all models file if it exists and optionally removes
        all other regular files in the output directory. Additionally
        resets ``self.all_models`` to an empty ``AllModels`` object.

        Parameters
        ----------
        wipe_other_files : Bool, optional
            If True, all regular files in the output directory will be
            deleted and a new backup of the config file will be created.
            If False, only the all models file will be removed.
            The default is False.

        Raises
        ------
        Exception if file(s) cannot be removed.

        Returns
        -------
        None.

        """
        if wipe_other_files:
            output_dir = self.settings.io_settings['output_directory']
            for f in glob.glob(f'{output_dir}*'):
                if os.path.isfile(f):
                    os.remove(f)
            # self.backup_config_file()
            self.logger.info(f'Removed files in {output_dir}.')
        else:
            all_models_file = self.settings.io_settings['output_directory'] \
                              + self.settings.io_settings['all_models_file']
            if os.path.isfile(all_models_file):
                os.remove(all_models_file)
                self.logger.info(f'Deleted existing {all_models_file}.')
        self.all_models = model.AllModels(config=self)
        self.logger.info('Instantiated empty AllModels object')
        self.logger.debug(f'AllModels:\n{self.all_models.table}')

    def remove_all_existing_output(self, wipe_all=False, create_tree=True):
        """
        Removes all existing DYNAMITE output.

        The options determine whether non-DYNAMITE output shall survive in the
        output folders. Also resets ``self.all_models`` to an empty
        ``AllModels`` object.

        Parameters
        ----------
        wipe_all : Bool, optional
            If True, the complete output directory tree will be removed.
            Set to False to keep (a) user files & directories in the output
            directory and (b) user directories in the plots directory.
            The default is False.
        create_tree : Bool, optional
            If True, recreates an empty output directory tree with a
            backup of the config file. Does not recreate directories
            if False. The default is True.

        Raises
        ------
        Exception if wipe_all==True and the output directory tree cannot
        be deleted.

        Returns
        -------
        None.

        """
        if wipe_all:
            out_dir = self.settings.io_settings['output_directory']
            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir)
            else:
                self.logger.warning(f'No output directory {out_dir} to '
                                    'remove.')
            self.logger.info(f'Output directory tree {out_dir} removed.')
        else:
            self.remove_existing_orblibs()
            self.remove_existing_plots(remove_directory=False)
        # Execute in any case to create empty AllModels object:
        self.remove_existing_all_models_file()
        if create_tree:
            self.make_output_directory_tree()
            # self.backup_config_file()

    def make_output_directory_tree(self):
        """
        Create output directory tree. Existing directories not affected.

        Returns
        -------
        None.

        """
        out_dir = self.settings.io_settings['output_directory']
        model_dir = self.settings.io_settings['model_directory']
        plot_dir = self.settings.io_settings['plot_directory']
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
            self.logger.debug(f'Output directory {out_dir} created.')
        else:
            self.logger.debug(f'Using existing output directory {out_dir}.')
        if not os.path.isdir(model_dir):
            os.mkdir(model_dir)
            self.logger.debug(f'Model directory {model_dir} created.')
        else:
            self.logger.debug(f'Using existing model directory {model_dir}.')
        if not os.path.isdir(plot_dir):
            os.mkdir(plot_dir)
            self.logger.debug(f'Plots directory {plot_dir} created.')
        else:
            self.logger.debug(f'Using existing plots directory {plot_dir}.')

    def copy_config_file(self, dest_directory, clean=True):
        """
        Copy config file to dest_directory.

        Creates a copy of the config file, intended to add it to the directory
        holding the model results. The file date will be preserved if possible.

        Parameters
        ----------
        dest_directory : str, mandatory
            The directory the config file will be copied to.
        clean : bool, optional
            If True, all `*`.yaml files in dest_directory will be deleted
            before copying. Default is True.
        """
        if dest_directory[-1] != '/':
            dest_directory += '/'
        if clean:
            del_files = glob.iglob(f'{dest_directory}*.yaml')
            for f_name in del_files:
                if os.path.isfile(f_name):
                    os.remove(f_name)
            self.logger.debug(f'{dest_directory}*.yaml files deleted.')
        shutil.copy2(self.config_file_name, dest_directory)
        self.logger.info('Config file copied to '
                         f'{dest_directory}{self.config_file_name}.')


    def backup_config_file(self, reset=False, keep=None, delete_other=False):
        """
        Copy the config file to the output directory.

        A running index of the format _xxx will be appended to the base file
        name to keep track of earlier config files (config_000.yaml,
        config_001.yaml, config_002.yaml, etc...)

        This method is not used in standard DYNAMITE and provided as a utility.

        Parameters
        ----------
        reset : bool, optional
            If reset==True, all `*`.yaml files in the output directory
            are deleted before the config file is copied.
            The default is False.
        keep : int or NoneType, optional
            If an integer > 0, at most `keep` config files WITH THE SAME
            BASE NAME are kept. Used to control the number of config file
            backups in the output folder. Can be combined with delete_other.
            The current config file will always be backuped. If keep==None,
            nothing is done. The default is None.
        delete_other : bool, optional
            If delete_other==True, all config files WITH A DIFFERENT BASE NAME
            will be deleted. If delete_other==False, nothing is done.
            The default is False.

        Returns
        -------
        None.

        """
        out_dir = self.settings.io_settings['output_directory']
        f_root,f_ext=os.path.splitext(os.path.basename(self.config_file_name))
        if reset:
            del_files = glob.iglob(f'{out_dir}*{f_ext}')
            for fname in del_files:
                if os.path.isfile(fname):
                    os.remove(fname)
            dest_file_name = f'{out_dir}{f_root}_000{f_ext}'
            self.logger.debug('Config file backup reset.')
        else:
            c_pattern = f'{out_dir}{f_root}_[0-9][0-9][0-9]{f_ext}'
            conf_files=glob.iglob(c_pattern)
            conf_roots = [os.path.splitext(i)[0] for i in conf_files]
            indices = [int(i[i.rindex('_')+1:]) for i in conf_roots]
            new_idx = max(indices) + 1 if len(indices)> 0 else 0
            dest_file_name = f'{out_dir}{f_root}_{new_idx:03d}{f_ext}'
            if keep is not None:
                if keep<1 or keep!=int(keep):
                    text = 'Parameter keep must be a positive integer.'
                    self.logger.error(text)
                    raise ValueError(text)
                for i in sorted(indices)[:-keep]:
                    os.remove(f'{out_dir}{f_root}_{i:03d}{f_ext}')
                self.logger.debug(f'{len(indices[:-keep])} config file(s) '
                                  'removed.')
            if delete_other:
                all_conf_files=glob.iglob(f'{out_dir}*_[0-9][0-9][0-9]{f_ext}')
                del_files = [f for f in all_conf_files
                             if not fnmatch.fnmatch(f,c_pattern)]
                for fname in del_files:
                    os.remove(fname)
                self.logger.debug(f'{len(del_files)} other config file(s) '
                                  'removed.')
        shutil.copy(self.config_file_name, dest_file_name)
        self.logger.info(f'Config file backup: {dest_file_name}.')

    def validate(self):
        """
        Validates the system and settings.

        This method is still VERY rudimentary and will be adjusted as we add new
        functionality to dynamite. Currently, this method is geared towards
        legacy mode.

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
        if not 1 < len(self.system.cmp_list) < 4:
            self.logger.error('System needs to comprise exactly one Plummer, '
                              'one VisibleComponent, and zero or one DM Halo '
                              'object(s)')
            raise ValueError('System needs to comprise exactly one Plummer, '
                             'one VisibleComponent, and zero or one DM Halo '
                             'object(s)')

        ws_type = self.settings.weight_solver_settings['type']

        for c in self.system.cmp_list:
            if issubclass(type(c), physys.VisibleComponent): # Check vis. comp.
                if c.kinematic_data:
                    for kin_data in c.kinematic_data:
                        check_gh = (kin_data.type == 'GaussHermite')
                        check_bl = (kin_data.type == 'BayesLOSVD')
                        if (not check_gh) and (not check_bl):
                            self.logger.error('VisibleComponent kinematics type'
                                              'must be GaussHermite or '
                                              'BayesLOSVD')
                            raise ValueError('VisibleComponent kinematics type'
                                             'must be GaussHermite or '
                                             'BayesLOSVD')
                        if check_bl:
                            # check weight solver type
                            if ws_type == 'LegacyWeightSolver':
                                self.logger.error("LegacyWeightSolver can't be "
                                                  "used with BayesLOSVD - use "
                                                  "weight-solver type NNLS")
                                raise ValueError("LegacyWeightSolver can't be "
                                                  "used with BayesLOSVD - use "
                                                  "weight-solver type NNLS")

                else:
                    self.logger.error('VisibleComponent must have kinematics: '
                                      'either GaussHermite or BayesLOSVD')
                    raise ValueError('VisibleComponent must have kinematics: '
                                     'either GaussHermite or BayesLOSVD')
                if c.symmetry != 'triax':
                    self.logger.error('Legacy mode: VisibleComponent must be '
                                      'triaxial')
                    raise ValueError('Legacy mode: VisibleComponent must be '
                                     'triaxial')
                continue
            if issubclass(type(c), physys.DarkComponent) \
                and not isinstance(c, physys.Plummer):
            # Check allowed dm halos in legacy mode
                if type(c) not in [physys.NFW, physys.NFW_m200_c,
                                   physys.Hernquist,
                                   physys.TriaxialCoredLogPotential,
                                   physys.GeneralisedNFW]:
                    text = 'DM Halo needs to be of type NFW, NFW_m200_c, ' \
                           'Hernquist, TriaxialCoredLogPotential, ' \
                           f'or GeneralisedNFW, not {type(c)}'
                    self.logger.error(text)
                    raise ValueError(text)

        gen_type = self.settings.parameter_space_settings["generator_type"]
        allowed_types = ['GridWalk', 'LegacyGridSearch', 'FullGrid']
        if gen_type not in allowed_types:
            text = f'Legacy mode: parameter space generator_type ' \
                   f'must be in {allowed_types}'
            self.logger.error(text)
            raise ValueError(text)
        chi2abs = self.__class__.thresh_chi2_abs
        chi2scaled = self.__class__.thresh_chi2_scaled
        gen_set=self.settings.parameter_space_settings.get('generator_settings')
        if gen_set != None and (chi2abs in gen_set and chi2scaled in gen_set):
            self.logger.error(f'Only specify one of {chi2abs}, {chi2scaled}, '
                              'not both')
            raise ValueError(f'Only specify one of {chi2abs}, {chi2scaled}, '
                             'not both')

        if ws_type == 'LegacyWeightSolver':
            # check velocity histograms settings if LegacyWeightSolver is used.
            # (i) check all velocity histograms have center 0, (ii) force them
            # all to have equal widths and (odd) number of bins
            # these requirements are not needed by orblib_f.f90, but are assumed
            # by the NNLS routine triaxnnl_*.f90 (see 2144-2145 of orblib_f.f90)
            # Therefore this check is based on WeightSolver type.
            stars = self.system.get_component_from_class(
                physys.TriaxialVisibleComponent
                )
            hist_widths = [k.hist_width for k in stars.kinematic_data]
            hist_centers = [k.hist_center for k in stars.kinematic_data]
            hist_bins = [k.hist_bins for k in stars.kinematic_data]
            self.logger.debug('checking all values of hist_center == 0...')
            assert all([x==0 for x in hist_centers]), 'all hist_center values must be 0'
            self.logger.debug('... check passed')
            equal_widths = all([x == hist_widths[0] for x in hist_widths])
            if equal_widths is False:
                max_width = max(hist_widths)
                msg = 'Value of `hist_width` must be the same for all kinematic'
                msg += f' data - defaulting to widest provided i.e. {max_width}'
                self.logger.info(msg)
                for k in stars.kinematic_data:
                    k.hist_width = max_width
            equal_bins = all([x == hist_bins[0] for x in hist_bins])
            if equal_bins is False:
                max_bins = max(hist_bins)
                msg = 'Value of `hist_bins` must be the same for all kinematic'
                msg += f' data - defaulting to largest provided i.e. {max_bins}'
                self.logger.info(msg)
                if max_bins%2 == 0:
                    msg = 'Value of `hist_bins` must be odd: '
                    msg += f'replacing {max_bins} with {max_bins+1}'
                    self.logger.info(msg)
                    max_bins += 1
                for k in stars.kinematic_data:
                    k.hist_bins = max_bins
        self.settings.validate()

        which_chi2 = self.validate_chi2()
        if which_chi2 == 'kinmapchi2' and ws_type != 'LegacyWeightSolver':
            msg = 'kinmapchi2 is only allowed with LegacyWeightSolver'
            self.logger.error(msg)
            raise ValueError(msg)

    def validate_chi2(self, which_chi2=None):
        """
        Validates which_chi2 setting

        Validates the which_chi2 setting in the config file (if argument
        which_chi2 is None) or the string given in the argument.

        Parameters
        ----------
        which_chi2 : str, optional
            If None, the which_chi2 setting in the config file is validated;
            if not None, the string given is validated. The default is None.

        Raises
        ------
        ValueError
            If which_chi2 fails validation.

        Returns
        -------
        which_chi2 : str
            The valid which_chi2 setting: either the value from the config
            file or the string passed as an argument.

        """
        allowed_chi2 = ('chi2', 'kinchi2', 'kinmapchi2')
        if which_chi2 == None:
            which_chi2 = self.settings.parameter_space_settings['which_chi2']
        if which_chi2 not in allowed_chi2:
            text = 'parameter_space_settings: which_chi2 must be one of ' \
                   f'{allowed_chi2}, not {which_chi2}.'
            self.logger.error(text)
            raise ValueError(text)
        return which_chi2


class DynamiteLogging(object):
    """Dynamite logging setup.

    ONLY use if logging has not been configured outside of Dynamite. If no
    arguments are given, the logging setup is as follows:

    -   log to the console with logging level INFO, messages include the
        level, timestamp, class name, and message text
    -   create a dynamite.log file with logging level DEBUG, messages
        include the level, timestamp, class name,
        filename:method:line number, and message text

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

    """
    def __init__(self, logfile='dynamite.log', console_level=logging.INFO,
                                               logfile_level=logging.DEBUG,
                                               console_formatter = None,
                                               logfile_formatter = None):
        logging.shutdown()
        importlib.reload(logging)
        logger = logging.getLogger()       # create logger
        logger.setLevel(logging.DEBUG)     # set level that's lower than wanted
        ch = logging.StreamHandler(stream=sys.stderr) # create console handler
        ch.setLevel(console_level)         # set console logging level
        # create formatter
        if console_formatter is None:
            console_formatter = logging.Formatter( \
                '[%(levelname)s] %(asctime)s - %(name)s - '
                '%(message)s', "%H:%M:%S")
        ch.setFormatter(console_formatter) # add the formatter to the handler
        logger.addHandler(ch)              # add the handler to the logger

        if logfile:
            fh = logging.FileHandler(logfile, mode='w') # create file handler
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
