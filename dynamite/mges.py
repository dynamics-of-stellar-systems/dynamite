import data
import numpy as np
from astropy import table
from astropy.io import ascii


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
