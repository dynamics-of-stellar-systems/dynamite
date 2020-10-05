from astropy.io import ascii


class Data(object):

    def __init__(self,
                 name=None,
                 datafile=None,
                 input_directory=None
                 ):
        self.name = name
        if not hasattr(self, input_directory):
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
        self.PSF = self.data.meta['PSF']
        pass





# end
