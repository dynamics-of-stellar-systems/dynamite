import numpy as np

grav_const_km = 6.67428e-11*1.98892e30/1e9
parsec_km = 1.4959787068e8*(648.000e3/np.pi)
rho_crit = (3.*((7.3000e-5)/parsec_km)**2)/(8.*np.pi*grav_const_km)
