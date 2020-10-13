import numpy as np

G = 6.67e-11
Msun = 1.989e30
kpc = 3.086e19
kms = 1000.
Guse = G * Msun / kpc
deg = 360./2./np.pi
arcsec = 360./2./np.pi * 60. * 60.
arcmin = arcsec / 60.
Msunkpc3toGeVcm3 = 3.7973132915271756e-08
kpccm = 1e3*3.0856775*1e18
rhocrit = 135.05
oden = 200
h = 0.7
year = 365.0*24.0*60.0*60.0
