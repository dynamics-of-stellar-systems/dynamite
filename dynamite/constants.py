import numpy as np

GRAV_CONST_KM = 6.67428e-11*1.98892e30/1e9
PARSEC_KM = 1.4959787068e8*(648.000e3/np.pi)
RHO_CRIT = (3.*((7.3000e-5)/PARSEC_KM)**2)/(8.*np.pi*GRAV_CONST_KM)

weight_file = 'orbit_weights.ecsv'