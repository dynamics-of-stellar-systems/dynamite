import logging

class Populations(object):
    """
    Populations class holding attributes and methods pertaining to population data
    """
    values = []
    def __init__(self,
                 name=None,
                 weight=None,
                 type=None,
                 datafile=None,
                 aperturefile=None,
                 binfile=None,
                 PSF=None,
                 ):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.name = name
        self.weight = weight
        self.type = type
        self.datafile = datafile
        self.aperturefile = aperturefile
        self.binfile = binfile
        self.PSF = PSF
        self.__class__.values = list(self.__dict__.keys())

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
                text = 'Invalid population key ' + k + '. Allowed keys: ' + \
                       str(tuple(self.__class__.values))
                self.logger.error(text)
                raise ValueError(text)
            setattr(self, k, v)

    def validate(self): # here we can put more validation...
        if sorted(self.__class__.values) != sorted(self.__dict__.keys()):
            text = 'Population attributes can only be ' + \
                   str(tuple(self.__class__.values)) + ', not ' + \
                   str(tuple(self.__dict__.keys()))
            self.logger.error(text)
            raise ValueError(text)

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__dict__})'

# end