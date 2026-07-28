import numpy as np

GRAV_CONST_KM = 6.67428e-11*1.98892e30/1e9  # unit km**3 / Msun / s**2
PARSEC_KM = 1.4959787068e8*(648.000e3/np.pi)  # unit km

weight_file = 'orbit_weights.ecsv'  # weights file
p_masses_file = 'mass_aper.ecsv'  # projected masses file

def ARC_KPC(distance):
    """Returns the conversion factor from arcseconds to kiloparsecs.
    The distance is in MPc, and the result is in kpc/arcsec."""
    return distance * np.pi / 648

def ARC_KM(distance):
    """Returns the conversion factor from arcseconds to km.
    The distance is in MPc, and the result is in km/arcsec."""
    return distance * 1e6 * np.tan(np.pi / 648e3) * PARSEC_KM
