import logging

from dynamite import data


class Populations(data.Integrated):
    """
    Populations class holding attributes and methods pertaining to population
    data


    Parameters
    ----------
    kin_aper : int or None
        If an integer, the index of the kinematic data set (starting with 0)
        which shares its apertures with the population data set.
        If None (the default), the population data set has its own apertures.
    """
    values = []
    def __init__(self,
                 weight=None,
                 hist_width='default',
                 hist_center='default',
                 hist_bins='default',
                 kin_aper=None,
                 pop_cols=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.kin_aper = kin_aper
        if hasattr(self, 'data'):
            if pop_cols is not None:
                self.clean_data(pop_cols)
            self.weight = weight
            self.type = type
            if hist_width=='default':
                self.set_default_hist_width()  # needed for pop data?
            else:
                self.hist_width = float(hist_width)
            if hist_center=='default':
                self.set_default_hist_center()  # needed for pop data?
            else:
                self.hist_center = float(hist_center)
            if hist_bins=='default':
                self.set_default_hist_bins()  # needed for pop data?
            else:
                self.hist_bins = int(hist_bins)
            self.__class__.values = list(self.__dict__.keys())
            self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
            if self.weight is None or self.hist_width is None or \
                    self.hist_center is None or self.hist_bins is None:
                text = 'Populations need (weight, hist_width, hist_center, '\
                   f'hist_bins), but has ({self.weight}, {self.type}, ' \
                   f'{self.hist_width}, {self.hist_center}, {self.hist_bins})'
                self.logger.error(text)
                raise ValueError(text)
            self.n_spatial_bins = len(self.data)

        # self.__class__.values = list(self.__dict__.keys())

    def update(self, **kwargs):
        """
        Updates one or more attributes, including consistency check

        Parameters
        ----------
        **kwargs : attribute=value pairs

        Raises
        ------
        ValueError
            If class attribute does not exist

        Returns
        -------
        None.

        """

        for k, v in kwargs.items():
            if k not in self.__class__.values:
                text = 'Invalid populations key ' + k + '. Allowed keys: ' + \
                       str(tuple(self.__class__.values))
                self.logger.error(text)
                raise ValueError(text)
            setattr(self, k, v)

    def validate(self): # here we can put more validation...
        if sorted(self.__class__.values) != sorted(self.__dict__.keys()):
            text = 'Populations attributes can only be ' + \
                   str(tuple(self.__class__.values)) + ', not ' + \
                   str(tuple(self.__dict__.keys()))
            self.logger.error(text)
            raise ValueError(text)

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__dict__})'

    def get_data(self):
        """Returns the populations data.

        This returns a deep copy of the self.data attribute.

        Returns
        -------
        astropy table
            The populations data

        """
        return self.data.copy(copy_data=True)

    def clean_data(self, pop_cols):
        """
        Removes all data columns except for the index and the populations data.

        Returns
        -------
        None.

        """
        self.data.remove_columns(c
                                 for c in self.data.columns
                                 if c not in pop_cols)
# end
