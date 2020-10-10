#######################################################################
#
# Copyright (C) 2004-2017, Michele Cappellari
# E-mail: michele.cappellari_at_physics.ox.ac.uk
#
# This software is provided as is without any warranty whatsoever.
# Permission to use, for non-commercial purposes is granted.
# Permission to modify for personal or internal use is granted,
# provided this copyright and disclaimer are included unchanged
# at the beginning of the file. All other rights are reserved.
#
#######################################################################
#
# NAME:
#   symmetrize_velfield()
#
# PURPOSE:
#   This routine generates a bi-symmetric ('axisymmetric') or point-symmetric
#   version of a given set of kinematic measurements.
#   PA: is the angle in degrees, measured counter-clockwise,
#       from the vertical axis (Y axis) to the galaxy major axis.
#       This is *not* necessary and not used when SYM = 3 or 4.
#   SYM: by-simmetry: is 1 for (V, h3, h5) and 2 for (sigma, h4, h6)
#     point-simmetry: is 3 for (V, h3, h5) and 4 for (sigma, h4, h6)
#
# HISTORY:
#   V1.0.0: Michele Cappellari, Vicenza, 21 May 2004
#   V1.0.1: Added MISSING keyword to TRIGRID call. Flipped velocity sign.
#       Written basic documentation. MC, Leiden, 25 May 2004
#   V1.1.0: Included point-symmetric case. Remco van den Bosch, Leiden, 18 January 2005
#   V1.1.1: Minor code revisions. MC, Leiden, 23 May 2005
#   V1.1.2: Important: changed definition of PA to be measured counterclockwise
#       with respect to the positive Y axis, as in astronomical convention and
#       consistently with my FIND_GALAXY routine. MC, Leiden, 1 June 2005
#   V1.1.3: Added optional keyword TRIANG. Corrected rare situation with w=-1.
#       MC, Leiden, 2 June 2005
#   V1.1.4: Added prefix SYMM_ to internal functions to prevent conflicts
#       with external functions with the same name. MC, Oxford, 11 May 2007
#   V2.0.0 : Completely rewritten without any loop. MC, Oxford, 8 October 2013
#   V2.0.1: Uses TOLERANCE keyword of TRIANGULATE to try to avoid IDL error
#       "TRIANGULATE: Points are co-linear, no solution." MC, Oxford, 2 December 2013
#   V3.0.0: Translated from IDL into Python. MC, Oxford, 14 February 2014
#   V3.0.1: Fixed rare case where interpolated value at boundary becomes
#       NaN due to numerical accuracy. MC, Oxford, 20 May 2015
#   V3.1.0: Re-implemented point-symmetric case. MC, Oxford, 1 November 2017
#
#######################################################################

import numpy as np
from scipy import interpolate

#----------------------------------------------------------------------
#     Michele cappellari, Paranal, 10 November 2013

def _rotate_points(x, y, ang):
    """
    Rotates points counter-clockwise by an angle ANG-90 in degrees.
    
    """    
    theta = np.radians(ang - 90.)
    xNew = x*np.cos(theta) - y*np.sin(theta)
    yNew = x*np.sin(theta) + y*np.cos(theta)

    return xNew, yNew
    
#----------------------------------------------------------------------
    
def symmetrize_velfield(xbin, ybin, vel_bin, sym=2, pa=90.):
    """
    This routine generates a bi-symmetric ('axisymmetric') of point-symmetric
    version of a given set of kinematical measurements.
    PA: is the angle in degrees, measured counter-clockwise,
      from the vertical axis (Y axis) to the galaxy major axis.
    SYM: by-simmetry: is 1 for (V, h3, h5) and 2 for (sigma, h4, h6)
      point-simmetry: is 3 for (V, h3, h5) and 4 for (sigma, h4, h6)
    
    """        
    xbin, ybin, vel_bin = map(np.asarray, [xbin, ybin, vel_bin])

    assert xbin.size == ybin.size == vel_bin.size, \
        "The vectors (xbin, ybin, velBin) must have the same size"
    assert isinstance(sym, int), "sym must be integer"
    assert 1 <= sym <= 4, "must be 1 <= sym <= 4"

    if sym < 3:
        x, y = _rotate_points(xbin, ybin, -pa)  # Negative PA for counter-clockwise
        xout = np.hstack([x,-x, x,-x])
        yout = np.hstack([y, y,-y,-y])
        vel_out = interpolate.griddata((x, y), vel_bin, (xout, yout))
        vel_out = vel_out.reshape(4, xbin.size)
        vel_out[0, :] = vel_bin  # see V3.0.1
        if sym == 1:
            vel_out[[1, 3], :] *= -1.
    else:
        vel_out = interpolate.griddata((xbin, ybin), vel_bin, (-xbin, -ybin))
        if sym == 3:
            vel_out = -vel_out
        vel_out = np.row_stack([vel_bin, vel_out])

    vel_sym = np.nanmean(vel_out, axis=0)
    
    return vel_sym

#----------------------------------------------------------------------
