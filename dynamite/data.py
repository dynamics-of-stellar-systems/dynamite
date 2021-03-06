from astropy.io import ascii
from astropy.table import Table
import logging

class Data(object):

    def __init__(self,
                 name=None,
                 datafile=None,
                 input_directory=None
                 ):
        self.name = name
        if not hasattr(self, 'input_directory'):
            self.datafile = datafile
            self.input_directory = input_directory if input_directory else ''
            if datafile is not None:
                self.data = ascii.read(self.input_directory+self.datafile)
            self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
            self.logger.debug(f'Data {self.name} read from '
                              f'{self.input_directory+self.datafile}')


class Discrete(Data):

    def __init__(self, x, y, values, errors):
        self.x = x
        self.y = y
        self.values = values


class Integrated(Data):

    def __init__(self,
                 aperturefile=None,
                 binfile=None,
                 maskfile=None,
                 **kwargs
                 ):
        self.aperturefile = aperturefile
        self.binfile = binfile
        self.maskfile = maskfile
        super().__init__(**kwargs)
        if hasattr(self, 'data'):
            self.PSF = self.data.meta['PSF']
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')

    def add_psf_to_datafile(self,
                            sigma=[1.],
                            weight=[1.],
                            datafile='datafile.ecsv'):
        """Write PSF into the datafile

        Quantities are written as meta-data under 'PSF'. PSF represented by a
        sequence of (Gaussian) sigma and weights.

        Parameters
        ----------
        sigma : list
            Gaussian std dev of PSF components
        weight : list
            weights of PSF components
        datafile : string
            Description of parameter `datafile`.

        """
        assert type(sigma) is list
        assert type(weight) is list
        assert isinstance(datafile, str)
        if hasattr(self, 'PSF'):
            self.logger.warning('Warning: this dataset already has an ' + \
                                'associated PSF! Possibly overwriting an ' + \
                                'existing PSF in the datafile')
        old_table = ascii.read(datafile)
        meta = old_table.meta
        psf = {'sigma':sigma, 'weight':weight}
        meta.update({'PSF':psf})
        new_table = Table(old_table, meta=meta)
        new_table.write(datafile, format='ascii.ecsv', overwrite=True)






# end
