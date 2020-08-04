# import numpy as np

class Data(object):

    def __init__(self):
        pass


class Discrete(Data):

    def __init__(self, x, y, values, errors):
        self.x = x
        self.y = y
        self.values = values


class Integrated(Data):

    def __init__(self, **kwargs):#aperture, values, errors):
        # self.aperture = aperture
        # self.values = values
        # self.errors = errors
        pass


# class Kinematic(Data):
#
#     def __init__(self):
#         pass
#
#     # def transform_orbits_to_observables(orb_lib):
#     #     # placeholder for code to transform orbits
#     #     return observables
#
#
# class GaussHermite(Integrated, Kinematic):
#
#     def __init__(self,
#                  filename=None):
#         pass
#
#     def transform_orbits_to_observables(orb_lib):
#         # actual code to transform orbits to GH coefficients
#         return gauss_hermite_coefficients


# class BSplines(Integrated, Kinematic):

#     def __init__(self,
#                  filename=None):
#         pass

#     def transform_orbits_to_observables(orb_lib):
#         # actual code to transform orbits to B-spline coefficients
#         return b_spline_coefficients


# class DiscreteLOS(Discrete, Kinematic):

#     def __init__(self, filename=None):
#         pass

#     def transform_orbits_to_observables(orb_lib):
#         # actual code to transform orbits to pdf evaluated for orbit libs
#         return pdf


# class Population(Data):

#     def __init__(self):
#         pass

#     def colouring_recipe(self, orb_lib, orb_wts):
#         # code to find orbit labels (colours) given observed self.values and an
#         # orbit library and orbit weights
#         # return orb_labels
#         pass


# class PopulationMap(Integrated, Population):

#     def __init__(self, filename=None):
#         pass

#     def colouring_recipe(self, orb_lib, orb_wts):
#         pass


# class DiscretePopulation(Discrete, Population):

#     def __init__(self, filename=None):
#         pass

#     def colouring_recipe(self, orb_lib, orb_wts):
#         pass




# end
