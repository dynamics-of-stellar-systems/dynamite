import data

import numpy as np
from scipy import special, stats
from astropy import table
from astropy.io import ascii

# TODO: move some of the kwargs from the init of 'Kinematics' to the init of
# higher level data classes, e.g. all Integrated objects will need aperturefile,
# binfile, maskfile, PSF.

# QUESTION: what does values list do...? Why is it outside of any method?

class Kinematics(data.Data):
    """
    Kinematics class holding attributes and methods pertaining to kinematics data
    """
    values = []
    def __init__(self,
                 name=None,
                 weight=None,
                 type=None,
                 parametrization=None,
                 datafile=None,
                 aperturefile=None,
                 binfile=None,
                 maskfile=None,
                 PSF=None,
                 ):
        self.name = name
        self.weight = weight
        self.type = type
        self.parametrization = parametrization
        self.datafile = datafile
        self.aperturefile = aperturefile
        self.binfile = binfile
        self.maskfile = maskfile
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
                raise ValueError('Invalid kinematics key ' + k + '. Allowed keys: ' + str(tuple(self.__class__.values)))
            setattr(self, k, v)

    def validate(self): # here we can put more validation...
        if sorted(self.__class__.values) != sorted(self.__dict__.keys()):
            raise ValueError('Kinematics attributes can only be ' + str(tuple(self.__class__.values)) + ', not ' + str(tuple(self.__dict__.keys())))

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__dict__})'


class GaussHermite(Kinematics, data.Integrated):

    def __init__(self, **kwargs):
        # super goes left to right, i.e. first calls "Kinematics" __init__, then
        # calls data.Integrated's __init__
        super(GaussHermite, self).__init__(**kwargs)
        # if filename is not None:
        data = self.read_file(self.datafile)
        self.vbin_id = data['vbin_id']
        self.v = data['v']
        self.dv = data['dv']
        self.sigma = data['sigma']
        self.sigma = data['dsigma']
        self.h3 = data['h3']
        self.dh3 = data['dh3']
        self.h4 = data['h4']
        self.dh4 = data['dh4']

    def read_file_old_format(self, filename):
        f = open(filename)
        header = f.readline()
        f.close()
        header = header.split()
        n_vbins, n_gh = header
        n_vbins, n_gh = int(n_vbins), int(n_gh)
        self.n_vbins = n_vbins
        self.n_gh = n_gh
        # n_gh read from the file header is currently unused and incorrect
        # TODO: fix/remove it... for now hardcode n_gh = 4
        self.n_gh = 4
        # TODO: move names/dtypes to the file header
        names = ['vbin_id', 'v', 'dv', 'sigma', 'dsigma', 'n_gh']
        for i_gh in range(3, self.n_gh+1):
            names += [f'h{i_gh}', f'dh{i_gh}']
        ncols = len(names)
        dtype = [float for i in range(ncols)]
        dtype[0] = int
        dtype[5] = int
        data = np.genfromtxt(filename,
                             skip_header=1,
                             names=names,
                             dtype=dtype)
        # 6th column of kin_data.dat is == (unused and incorrect?) n_gh
        # TODO: remove it
        return data

    def convert_file_from_old_format(self,
                                     filename_old_format,
                                     filename_new_format):
        data = self.read_file_old_format(filename_old_format)
        data = table.Table(data)
        data.remove_column('n_gh')
        data.write(filename_new_format, format='ascii.ecsv')
        return

    def read_file(self, filename):
        data = ascii.read(filename)
        return data


# end
