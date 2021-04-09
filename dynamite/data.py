from astropy.io import ascii
from astropy.table import Table
import logging
import numpy as np
from plotbin import display_pixels

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
                              f'{self.input_directory}{self.datafile}')


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
        if self.aperturefile is not None and self.binfile is not None:
            self.read_aperture_and_bin_files()
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

    def read_aperture_and_bin_files(self):
        '''read aperture and bin files and store the required data for plotting
        with plotbin.display_pixels to the dictionary self.dp_args'''
        # read aperture file
        aperture_fname = self.input_directory+self.aperturefile
        lines = [line.rstrip('\n').split() for line in open(aperture_fname)]
        strhead = lines[0]
        minx = float(lines[1][0])
        miny = float(lines[1][1])
        sx = float(lines[2][0])
        sy = float(lines[2][1])
        maxx = sx + minx
        sy = sy + miny
        angle_deg = float(lines[3][0])
        nx = int(lines[4][0])
        ny = int(lines[4][1])
        dx = sx / nx
        xr = np.arange(nx, dtype=float) * dx + minx + 0.5 * dx
        yc = np.arange(ny, dtype=float) * dx + miny + 0.5 * dx
        xi = np.outer(xr, (yc * 0 + 1))
        xt = xi.T.flatten()
        yi = np.outer((xr * 0 + 1), yc)
        yt = yi.T.flatten()
        radeg = 57.2958
        xi = xt
        yi = yt
        # read bin file
        bin_fname = self.input_directory+self.binfile
        lines_bins = [line.rstrip('\n').split() for line in open(bin_fname)]
        i = 0
        str_head = []
        i_var = []
        grid = []
        while i < len(lines_bins):
            for x in lines_bins[i]:
                if i == 0:
                    str_head.append(str(x))
                if i == 1:
                    i_var.append(int(x))
                if i > 1:
                    grid.append(int(x))
            i += 1
        str_head = str(str_head[0])
        i_var = int(i_var[0])
        grid = np.ravel(np.array(grid))
        # bins start counting at 1 in fortran and at 0 in idl:
        grid = grid - 1
        # Only select the pixels that have a bin associated with them.
        s = np.ravel(np.where((grid >= 0)))
        x, y = xi[s], yi[s]
        # store the arguments needed to use `plotbin.display_pixels`
        dp_args = {'x':x,
                   'y':y,
                   'dx':dx,
                   'idx_bin_to_pix':grid[s],
                   'angle':angle_deg}
        self.dp_args = dp_args

    def get_map_plotter(self):
        """Get a function which can plot maps of this dataset

        Returns
        -------
        function
            a function `map plotter` which plots kinematic maps

        """
        def map_plotter(bin_data, **kw_display_pixels):
            """plot maps with the spatial binning of this dataset

            This is a wrapper around plotbin.display_pixel

            Parameters
            ----------
            bin_data : array (nbins,)
                the binned data to be plotted with entries ordered as the bins
                are ordered in the table self.data
            **kw_display_pixels : type
                keyword arguments passed directly to plotbin.display_pixel

            Returns
            -------
            None

            """
            pix_data = bin_data[self.dp_args['idx_bin_to_pix']]
            display_pixels.display_pixels(self.dp_args['x'],
                                          self.dp_args['y'],
                                          pix_data,
                                          pixelsize=self.dp_args['dx'],
                                          angle=self.dp_args['angle'],
                                          **kw_display_pixels)
            return
        return map_plotter


# end
