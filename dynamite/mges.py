import logging
import numpy as np
from astropy import table
from dynamite import data

class MGE(data.Data):
    """Multi Gaussian Expansions"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.validate_q_values()

    def validate_q_values(self):
        """Validates the mge's q values

        Any q 'too close to 1' will be set to q=NINES for numerical stability.
        If any changes are made, a warning message will be logged. Note that
        the 'closeness' to 1 might be machine dependent - for us 'four nines',
        0.9999, worked...

        Returns
        -------
        None. Any changes are applied to ``self.data``.

        """
        NINES = 0.9999
        new_mge = False
        for r in self.data:
            if r['q'] > NINES:
                self.logger.warning(f'changing q={r["q"]} to q={NINES} for '
                                    'numerical stability.')
                r['q'] = NINES
                new_mge = True
        if new_mge:
            self.logger.warning(f'New mge:\n{self.data}')

    def read_file_old_format(self, filename):
        """read the MGE data from a text file

        old format = text file, 4 columns (I, sigma, q, PA_twist), with one
        header line

        Parameters
        ----------
        filename : string
            name of file

        Returns
        -------
        ndarray
            MGE data in structrued numpy array

        """
        with open(filename) as fp:
            n_cmp_mge = fp.readline()
        n_cmp_mge = int(n_cmp_mge)
        dat = np.genfromtxt(filename,
                            skip_header=1,
                            max_rows=n_cmp_mge,
                            names=['I', 'sigma', 'q', 'PA_twist'])
        if np.isnan(dat).any():
            txt = f'Input file {filename} has nans'
            self.logger.error(f'{txt} at: {np.argwhere(np.isnan(dat))}.')
            raise ValueError(txt)
        return dat

    def convert_file_from_old_format(self,
                                     filename_old_format,
                                     filename_new_format):
        """convert old mge file to ECSV file

        Parameters
        ----------
        filename_old_format : string
            old filename
        filename_new_format : string
            new filename

        Returns
        -------
        None
            saves file with name ``filename_new_format``

        """
        data = self.read_file_old_format(filename_old_format)
        data = table.Table(data)
        data.write(filename_new_format, format='ascii.ecsv')

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__dict__})'

    def get_projected_masses(self, parset, apertures):
        # TODO:
        # calculate the mass of the mge in observed 2D apertures given the
        # parameter set containing intrinsic axis ratios (p, q, u)
        # for now, use legacy implementation below which reads from file
        pass

    def get_projected_masses_from_file(self, directory_noml):
        """read mge projected masses from ``mass_aper.dat``

        Parameters
        ----------
        directory_noml : string
            name of model directory exclusing the ``ml/`` extension

        Returns
        -------
        array
            array of aperture masses of the MGE

        """
        fname = f'{directory_noml}datfil/mass_aper.dat'
        aperture_masses = np.loadtxt(fname, skiprows=1)
        # remove first column (aperture index)
        aperture_masses = aperture_masses[:,1]
        return aperture_masses

    def get_intrinsic_masses(self, parset, grid):
        # TODO: reimplement
        # calculate the mass of the mge in observed 3D grid given the
        # parameter set containing intrinsic axis ratios (p, q, u)
        # for now, use legacy implementation below which reads from file
        pass

    def get_intrinsic_masses_from_file(self, directory_noml):
        """read mge intrinsic masses from ``mass_qgrid.dat``

        Parameters
        ----------
        directory_noml : string
            name of model directory exclusing the ``ml/`` extension

        Returns
        -------
        array
            3D intrinsic_masses masses of the MGE in a polar grid with sizes
            (n_r, n_theta, n_phi) which are defined in the config file.
            Their defaults are (6,6,10)

        """
        fname = f'{directory_noml}datfil/mass_qgrid.dat'
        shape = np.loadtxt(fname, max_rows=1, dtype=int)
        intrinsic_masses = np.loadtxt(fname, skiprows=1)
        intrinsic_masses = np.reshape(intrinsic_masses, shape)
        return intrinsic_masses

    def __add__(self,other):
        """Concatenate two MGEs, preserving row order.

        The input_directory and filename attributes are inherited from the first MGE.

        Parameters
        ----------
        other : Object of type MGE
            the MGE to add to this one

        Returns
        -------
        new_mge : Object of type MGE

        """
        mge1_data = self.data
        mge2_data = other.data
        mge1_data['row_merge_ID'] = list(range(1,len(mge1_data)+1))
        mge2_data['row_merge_ID'] = list(range(len(mge1_data)+1,len(mge1_data)+len(mge2_data)+1))

        new_data = table.join(mge1_data, mge2_data, join_type='outer')
        new_data.sort('row_merge_ID')
        new_data.remove_columns('row_merge_ID')

        new_mge = MGE(input_directory=self.input_directory, datafile=self.datafile)
        new_mge.data = new_data

        return new_mge



# end
