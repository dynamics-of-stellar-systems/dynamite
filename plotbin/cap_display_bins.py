"""
NAME:
    display_bins(x, y, bin_num, vel_bin)
    
AUTHOR:
    Michele Cappellari, University of Oxford
    E-mail: michele.cappellari_at_physics.ox.ac.uk

PURPOSE:
    This simple routine illustrates how to display a Voronoi binned map.
    Keyword parameters will be passed to display_pixels() (see help there).
    
INPUTS:
    (x, y): (length npix) Coordinates of the original spaxels before binning;
    bin_num: (length npix) Integer bin number corresponding to each (x, y) pair,
            as provided in output by the voronoi_2d_binning() routine;
            The index goes from zero to nbins-1.
    vel_bin: (length nbins) Quantity associated to each bin, resulting
            e.g. from the kinematic extraction from the binned spectra.
          
MODIFICATION HISTORY:
    V1.0.0: Michele Cappellari, Oxford, 15 January 2015
    V1.0.1: Further input checks. MC, Oxford, 15 July 2015
    V1.0.2: Raise an error if bin_num is not integer. Pass kwargs to display_pixels.
        Thanks to Rebekka Schupp (MPIA) for the feedback.
        MC, Oxford, 31 July 2017
    V1.0.3: Changed imports for plotbin as a package. MC, Oxford, 17 April 2018
    
"""

import numpy as np

from plotbin.cap_display_pixels import display_pixels

def display_bins(x, y, bin_num, vel_bin, **kwargs):

    assert bin_num.dtype.kind == 'i', "bin_num must be integer"
    assert x.size == y.size == bin_num.size, "The vectors (x, y, bin_num) must have the same size"
    assert np.unique(bin_num).size == vel_bin.size, "vel_bin size does not match number of bins"
        
    img = display_pixels(x, y, vel_bin[bin_num], **kwargs)
    
    return img

