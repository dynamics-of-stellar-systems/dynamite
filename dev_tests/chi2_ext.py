# Example for a class implementing an external, additive chi2.

# The name of this module needs to be specified in the DYNAMITE config file
# in the chi2_ext component's definition, e.g. ext_module: "chi2_ext".

import logging

class Chi2Ext:
    """Class that implements calculating an external, additive chi2

    It is intended for use with DYNAMITE for components that are completely
    independent of DYNAMITE's orbit and weight calculations,
    such as gas kinematics.

    The class name is give in the DYNAMITE config file (ext_class).

    Parameters
    ----------
    arg1 : type and name as defined in the DYNAMITE config file (ext_class_args)
        _description_
    arg2 : type and name as defined in the DYNAMITE config file (ext_class_args)
        _description_
    """
    def __init__(self, arg1, arg2):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.logger.debug(f'Instantiated with parameters {arg1=}, {arg2=}.\n'
                          f'{type(arg1)=}, {type(arg2)=}')
        # ...
        # Code executed when instantiating the class. This happens once per
        # DYNAMITE run when reading the config file
        #...

    def chi2(self, model_id, config):
        """Method calculating chi2 which is added to all DYNAMITE chi2 values

        The method name is give in the DYNAMITE config file (ext_chi2).

        Parameters
        ----------
        model_id : int
            Model ID in the all_models table.
        config : a ``dyn.config_reader.Configuration`` object

        Returns
        -------
        float
            The external component's chi2 value that is added to all three
            "chi2", "kinchi2", and "kinmapchi2" values after DYNAMITE weight
            solving completes.
        """
        parset = dict(config.all_models.get_parset_from_row(row_id=model_id))
        ml_vscale = config.all_models.get_model_velocity_scaling_factor(
            model_id=model_id)
        self.logger.debug(f'This is chi2() and I got {parset=} from model '
                          f'{model_id} and the configuration for system '
                          f'{config.params["system_attributes"]["name"]}. '
                          f'The velocity scaling factor is {ml_vscale}.')
        # ...
        # Code that calculates the external, additive chi2
        chi2 = 42.
        # ...
        return chi2
