class MGE(object):

    def __init__(self):
        pass
    # def __init__(self, n, I, sigma, q):
    #     self.n = n
    #     self.I = I				    # intensity e.g. L/pc^2, kpc^-2, M/m^2
    #     self.sigma = sigma			# angular sigma
    #     self.q = q				    # flattenings

    # def surface_brightness(self, x, y):
    #     # function to calculate surface_brightness at (x,y)
    #     pass

    # def deproject(self, inc):
    #     # function to deproject MGE given an
    #     return intrinsic_MGE()


class MGE_from_file(MGE):

    def __init__(self, mge_filename):
        # tmp = read_file(mge_filename)
        # super(MGE_from_file, self).__init__(n = tmp['n'],
        #                                     I = tmp['I'],
        #                                     sigma = tmp['sigma'],
        #                                     q = tmp['q'])
        pass

# class MGE_from_image(MGE):

#     def __init__(self, image_filename):
#         # # img = read_image(image_filename)
#         # # tmp = self.fit_mge_to_image(img)
#         # super(MGE_from_image, self).__init__(n = tmp['n'],
#         #                                      I = tmp['I'],
#         #                                      sigma = tmp['sigma'],
#         #                                      q = tmp['q'])
#         pass

#     def fit_mge_to_image():
#         # e.g. Capellari code here
#         pass


# class intrinsic_MGE(MGE):
#     # intrinsic 3D density

#     def __init__(self, n, I, sigma, q):
#         super(intrinsic_MGE, self).__init__(n = tmp['n'],
#                                             I = tmp['I'],
#                                             sigma = tmp['sigma'],
#                                             q = tmp['q'])


# class intrinsic_MGE_from_xyz_grid(intrinsic_MGE):
#     # intrinsic 3D density

#     def __init__(self, xyz_grid, rho):
#         tmp = self.fit_mge()
#         super(intrinsic_MGE, self).__init__(n = tmp['n'],
#                                             I = tmp['I'],
#                                             sigma = tmp['sigma'],
#                                             q = tmp['q'])

#     def fit_mge(self, xyz_grid, rho):
#         # code to fit MGE given density evaluated at grid of (x,y,z) points
#         return




# end
