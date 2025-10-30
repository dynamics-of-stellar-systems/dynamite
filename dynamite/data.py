import logging
from astropy.io import ascii
from astropy.table import Table
import numpy as np
from plotbin import display_pixels

class Data(object):
    """Abstract class for data in Astropy ECSV files

    The data is stored in the Astropy table at ``self.data``

    Parameters
    ----------
    name : string
        Descriptve name of this data set
    datafile : string
        name of the Astropy ECSV datafile
    input_directory : string, or None
        location of the data file

    """

    def __init__(self,
                 name=None,
                 datafile=None,
                 input_directory=None
                 ):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.name = name
        if not hasattr(self, 'input_directory'):
            self.datafile = datafile
            self.input_directory = input_directory if input_directory else ''
            if datafile is not None:
                self.data = ascii.read(self.input_directory+self.datafile)
                self.logger.debug(f'Data {self.name} read from '
                                  f'{self.input_directory}{self.datafile}')
                data_array = np.lib.recfunctions.structured_to_unstructured(
                    self.data.as_array())
                if np.isnan(data_array).any():
                    txt=f'Input file {self.input_directory}{datafile} has nans'
                    self.logger.error(f'{txt} at: '
                                      f'{np.argwhere(np.isnan(data_array))}.')
                    raise ValueError(txt)


class Discrete(Data):
    """Class for discrete data

    # TODO: make this!

    """
    def __init__(self, x, y, values, errors):
        self.x = x
        self.y = y
        self.values = values


class Integrated(Data):
    """Abstract class for integrated data

    e.g. kinematic maps, population maps. TODO: deal with mask files!

    Parameters
    ----------
    aperturefile : string
        the name of the ``aperture.dat`` file of this dataset
    binfile : string
        the name of the ``bins.dat`` file of this dataset
    maskfile : string
        the name of the maskfile of this dataset
    **kwargs :
        other keyword arguments
    """
    def __init__(self,
                 aperturefile=None,
                 binfile=None,
                 **kwargs
                 ):
        self.aperturefile = aperturefile
        self.binfile = binfile
        super().__init__(**kwargs)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if hasattr(self, 'data'):
            self.PSF = self.data.meta['PSF']
            if abs(sum(self.PSF['weight'])-1.0) > 1e-8:
                txt = f"PSF weights add up to {sum(self.PSF['weight'])}, " + \
                      "not 1.0."
                if hasattr(self, 'datafile'):
                    txt += ' Check input data in '
                    if hasattr(self, 'input_directory'):
                        txt += f'{self.input_directory}'
                    txt += f'{self.datafile}.'
                self.logger.error(txt)
                raise ValueError(txt)
        if self.aperturefile is not None and self.binfile is not None:
            self.read_aperture_and_bin_files()

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
            output filename

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
        """read aperture and bin files

        Read the two files and store the required data for plotting
        with plotbin.display_pixels to the dictionary self.dp_args'''

        Returns
        -------
        Sets the attribute ``self.dp_args``

        """
        # read aperture file
        aperture_fname = self.input_directory + self.aperturefile
        lines = [line.rstrip('\n').split() for line in open(aperture_fname)
                                           if line[0] != '#']
        minx = float(lines[0][0])
        miny = float(lines[0][1])
        sx = float(lines[1][0])
        sy = float(lines[1][1])
        sy = sy + miny
        angle_deg = float(lines[2][0])
        nx = int(lines[3][0])
        ny = int(lines[3][1])
        dx = sx / nx
        xr = np.arange(nx, dtype=float) * dx + minx + 0.5 * dx
        yc = np.arange(ny, dtype=float) * dx + miny + 0.5 * dx
        xi = np.outer(xr, (yc * 0 + 1))
        xt = xi.T.flatten()
        yi = np.outer((xr * 0 + 1), yc)
        yt = yi.T.flatten()
        xi = xt
        yi = yt
        # read bin file
        bin_fname = self.input_directory + self.binfile
        lines_bins = [line.rstrip('\n').split() for line in open(bin_fname)
                                                if line[0] != '#']
        i = 0
        i_var = []
        grid = []
        while i < len(lines_bins):
            for x in lines_bins[i]:
                if i == 0:
                    i_var.append(int(x))
                if i > 0:
                    grid.append(int(x))
            i += 1
        i_var = int(i_var[0])
        grid = np.ravel(np.array(grid))
        if not (nx * ny == i_var == len(grid)):
            txt = f'Numbers of apertures do not match: {nx}*{ny}={nx*ny} ' \
                  f'({self.aperturefile}), {i_var} ({self.binfile}), ' \
                  f'{len(grid)} (bin data points in {self.binfile}).'
            self.logger.error(txt)
            raise ValueError(txt)
        self.logger.debug(f'{self.aperturefile} and {self.binfile} read.')
        n_bins_kinem = self.data[-1][0] + (1 if self.data[0][0] == 0 else 0)
        if not (n_bins_kinem == len(self.data) == max(grid)):
            txt = f'Numbers of kinematic bins do not match: {len(self.data)}'\
                  f' (length of {self.datafile}), {n_bins_kinem} (last id in '\
                  f'{self.datafile}), max number {max(grid)} in {self.binfile}.'
            self.logger.error(txt)
            raise ValueError(txt)
        self.logger.debug(f'Number of vbins in {self.datafile}, '
                          f'{self.binfile} validated.')
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
                   'angle':angle_deg} # Angle in degrees measured counter
                                      # clockwise from the galaxy major axis
                                      # to the X-axis of the input data
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

    def convert_to_plot_coords(self, x, y):
        """Convert a set of co-ordinates (x,y) into those used for plotting.

        i.e. rotated so major-axis lies along x-axis

        Parameters
        ----------
        x : float, or array like
            x-coordinates
        y : float, or array like
            x-coordinates

        Returns
        -------
        tuple
            rotated (x,y) co-ordinates which can be used on plots

        """
        ang = np.radians(self.dp_args['angle'])
        x, y = x*np.cos(ang) - y*np.sin(ang), x*np.sin(ang) + y*np.cos(ang)
        return x, y

# end
