"""

V1.0.0: Written. Michele Cappellari, Oxford, 04 April 2014
V1.0.1: Included angle keyword and kwargs.
    MC, Campos do Jordao, Brazil, 23 November 2015
V1.0.2: Updated documentation and included usage warning,
    MC, Oxford, 11 January 2015
V1.0.3: Use for loop with large arrays to reduce memory usage.
    MC, Oxford, 5 July 2017
V1.0.4: Changed imports for plotbin as a package. MC, Oxford, 17 April 2018    

"""

import warnings
import numpy as np

from plotbin.cap_display_pixels import display_pixels

def display_bins_generators(xBin, yBin, velBin, x, y, angle=None, **kwargs):
    """
    Displays a Voronoi binned map starting from the original coordinates of the pixels
    and the coordinates of the *generators* (not the centroids!) of the Voronoi
    tessellation, as provided in output e.g. by my voronoi_2d_binning routine.

    NB: When possible, instead of this routine, one should use the more general display_bins
    routine which uses the binNumber of every spaxel instead of the Voronoi generators.

    :param xBin: coordinates of the *generators* of the Voronoi tessellation
    :param yBin:
    :param velBin:
    :param x: coordinates of the original spaxels
    :param y:
    :param angle:
    :param kwargs:
    :return: image
    """

    warnings.warn('When possible, usage of the routine display_bins is preferred to display_bins_generators')

    assert xBin.size == yBin.size == velBin.size, 'The vectors (XBIN, YBIN, VEL) must have the same size'
    assert x.size == y.size, 'The vectors (X, Y) must have the same size'
    assert x.size >= xBin.size, 'The vectors (X, Y) cannot be smaller than (XBIN, YBIN)'
    
    # Perform a Voronoi tessellation starting from the coordinates
    # of the generators and the coordinates of the original pixels
    #
    if x.size < 1e4:
        binNum = np.argmin((x[:, None] - xBin)**2 + (y[:, None] - yBin)**2, axis=1)
    else:  # use for loop to reduce memory usage
        binNum = np.zeros(x.size, dtype=int)
        for j, (xj, yj) in enumerate(zip(x, y)):
            binNum[j] = np.argmin((xj - xBin)**2 + (yj - yBin)**2)

    f = display_pixels(x, y, velBin[binNum], angle=angle, **kwargs)
    
    return f
    

