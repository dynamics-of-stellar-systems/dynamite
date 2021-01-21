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
        pass

    def add_psf_to_datafile(self,
                            sigma=[1.],
                            weight=[1.],
                            datafile='datafile.ecsv'):
        logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        assert type(sigma) is list
        assert type(sigma) is list
        assert isinstance(datafile, str)
        if hasattr(self, 'PSF'):
            logger.warning('Warning: this dataset already has an ' + \
                           'associated PSF')
            logger.warning('Possibly overwriting an existing PSF in the ' + \
                           'datafile')
            # print('Warning: this dataset already has an associated PSF')
            # print('Possibly overwriting an existing PSF in the datafile')
        psf = {'sigma':sigma, 'weight':weight}
        meta = {'PSF':psf}
        old_table = ascii.read(datafile)
        new_table = Table(old_table, meta=meta)
        new_table.write(datafile, format='ascii.ecsv', overwrite=True)






# end
