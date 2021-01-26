import data
import numpy as np
from astropy import table

class MGE(data.Data):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def read_file_old_format(self, filename):
        with open(filename) as fp:
            n_cmp_mge = fp.readline()
        n_cmp_mge = int(n_cmp_mge)
        dat = np.genfromtxt(filename,
                            skip_header=1,
                            max_rows=n_cmp_mge,
                            names=['I', 'sigma', 'q', 'PA_twist'])
        return dat

    def convert_file_from_old_format(self,
                                     filename_old_format,
                                     filename_new_format):
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
        fname = f'{directory_noml}datfil/mass_qgrid.dat'
        shape = np.loadtxt(fname, max_rows=1, dtype=int)
        intrinsic_masses = np.loadtxt(fname, skiprows=1)
        intrinsic_masses = np.reshape(intrinsic_masses, shape)
        return intrinsic_masses

# TODO:
# class MGE_from_image(MGE):
#
#     def __init__(self, image_filename):
#         # # img = read_image(image_filename)
#         # # tmp = self.fit_mge_to_image(img)
#         # super(MGE_from_image, self).__init__(n = tmp['n'],
#         #                                      I = tmp['I'],
#         #                                      sigma = tmp['sigma'],
#         #                                      q = tmp['q'])
#         pass
#
#     def fit_mge_to_image():
#         # e.g. Capellari code here
#         pass


# TODO:
# class intrinsic_MGE(MGE):
#     # intrinsic 3D density
#
#     def __init__(self, n, I, sigma, q):
#         super(intrinsic_MGE, self).__init__(n = tmp['n'],
#                                             I = tmp['I'],
#                                             sigma = tmp['sigma'],
#                                             q = tmp['q'])


# TODO:
# class intrinsic_MGE_from_xyz_grid(intrinsic_MGE):
#     # intrinsic 3D density
#
#     def __init__(self, xyz_grid, rho):
#         tmp = self.fit_mge()
#         super(intrinsic_MGE, self).__init__(n = tmp['n'],
#                                             I = tmp['I'],
#                                             sigma = tmp['sigma'],
#                                             q = tmp['q'])
#
#     def fit_mge(self, xyz_grid, rho):
#         # code to fit MGE given density evaluated at grid of (x,y,z) points
#         return




# end
