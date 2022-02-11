import numpy as np
from constants import * 
from functions import * 

#This file contains all the code options and choices for 
#running a given model. Throughout, -1 means auto-calculate.

#Data files and output base filename:
whichgal = 'Ocen'
dgal_kpc = 5.4
infile = output_base+whichgal+'/dgal_%.1f/' % (dgal_kpc)+whichgal
outdirbase = output_base+whichgal+'/dgal_%.1f/' % (dgal_kpc)
print('Using Omega cen distance: %.1f kpc' % (dgal_kpc))

#Plot ranges and sample points [-1 means auto-calculate]:
rplot_inner = 1e-5
rplot_outer = 5.0
rplot_pnts = 2500
y_sigLOSmax = 25
ymin_Sigstar = 1e-3
ymax_Sigstar = 1e5
yMlow = 1e2
yMhigh = 1e9
yrholow = 1e3
yrhohigh = 1e12
alp3sig = 0.0
sigmlow = 1e-3
sigmhigh = 5.0

#Code options:
propermotion = 'yes'
virialshape = 'no'

#Overplot true solution (for mock data). If 
#yes, then the true solutions should be passed
#in: ranal,betatrue(ranal),betatruestar(ranal),
#truemass(ranal),trueden(ranal),truedlnrhodlnr(ranal).
overtrue = 'no'

#Radial grid range for Jeans calculation:
rmin = 5.0e-5
rmax = 5.0

#Galaxy properties. Assume here that the baryonic mass
#has the same radial profile as the tracer stars. If this
#is not the case, you should set Mstar_rad and Mstar_prof 
#here. The variables barrad_min, barrad_max and bar_pnts 
#set the radial range and sampling of the baryonic mass model.
Mstar = 3.0e6
Mstar_err = Mstar*0.75
baryonmass_follows_tracer = 'yes'
barrad_min = 1.0e-5
barrad_max = 5.0
bar_pnts = 2500


###########################################################
#Priors

#For surface density fit tracertol = [0,1] sets the spread 
#around the best-fit value from the binulator.
tracertol = 0.75

#Cosmology priors on the coreNFWtides model. mWDM(keV) is
#the mass of a thermal relic; <0 means CDM; sig_c200 is 
#the scatter of c200 in log10 space. If the cosmo_cprior
#is set, then we include a Gaussian spread in M200-c200 in
#the likelihood. Without this, M200-c200 enters only if 
#used to set the priors, below.
cosmo_cprior = 'yes'
sig_c200 = 0.1
mWDM = -1
if (mWDM > 0):
    cosmo_cfunc = lambda M200,h : \
        cosmo_cfunc_WDM(M200,h,OmegaM,rhocrit,mWDM)

#Velocity anisotropy priors:
betr0min = -3
betr0max = 0.0
betnmin = 1.0
betnmax = 3.0
bet0min = -0.4
bet0max = 0.4
betinfmin = -0.4
betinfmax = 1.0

#CoreNFWtides priors:
logM200low = 8.5
logM200high = 10.5
clow = 5.0
chigh = 50.0
rclow = 1e-3
rchigh = 2.0
logrclow = np.log10(rclow)
logrchigh = np.log10(rchigh)
nlow = 0.0
nhigh = 1.0

#Use here Jacobi radius from:
#https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1127K/abstract
rtlow = 160.0 / 1000.0
rthigh = 200.0 / 1000.0
logrtlow = np.log10(rtlow)
logrthigh = np.log10(rthigh)
dellow = 4.0
delhigh = 8.0

if (cosmo_cprior == 'yes'):
    clow = 1.0
    chigh = 100.0

#Priors on central dark mass [set logMcenlow/high very negative
#to switch this off. Mcen is the mass in Msun; acen is the
#scale length in kpc, usually assumed smaller than Rhalf
#to avoid degeneracies with the stellar mass]:
logMcenlow = 1.0
logMcenhigh = 6.0
acenlow = 1.0e-6
acenhigh = 1.0e-3

#Priors on rotation [set Arotlow = Arothigh = 0 to switch off.
#Arot defined s.t. vphimean^2 / (2 sigr^2) = Arot(r/Rhalf)
#which yields linear rotation with radius. Arot = 0.5 means
#an equal balance of rotation and pressure support at Rhalf].
#Note that:
#https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.1460S/abstract
#find a mean rotation amplitude of 4.27+/-0.52 km/s. At Rhalf,
#the dispersion of Ocen is 12-14km/s. This yields Arot in the
#range: 0.035 < Arot < 0.08:
#Arotlow = 0.035
#Arothigh = 0.08

#Switch off rotation, however, in favour of simply adding it
#in quadratude to the dispersion profile (usual trick):
Arotlow = 0.0
Arothigh = 1.0e-12

#Priors on distance [set drangelow = drangehigh = 1.0 to switch
#off. True distance follows as dgal_kpc * drange s.t. we
#usually want drangelow < 1.0 and drangehigh > 1.0]:

#Priors here based on other distance determinations, e.g.:
#https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.5957B/abstract
drangelow = 0.925
drangehigh = 1.05

###########################################################
#Post processing options:

#For calculating J-factors:
calc_Jfac = 'no'
alpha_Jfac_deg = 0.5
