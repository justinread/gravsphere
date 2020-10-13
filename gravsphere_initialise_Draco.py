import numpy as np
from constants import * 
from functions import * 

#This file contains all the code options and choices for 
#running a given model. Throughout, -1 means auto-calculate.

#Data files and output base filename:
whichgal = 'Draco'
infile = './Output/'+whichgal+'/'+whichgal
outdirbase = './Output/'+whichgal+'/'

#Plot ranges and sample points [-1 means auto-calculate]:
rplot_inner = -1
rplot_outer = -1
rplot_pnts = 100
y_sigLOSmax = 15
ymin_Sigstar = 1e-6
ymax_Sigstar = 100
yMlow = 1e4
yMhigh = 1e10
yrholow = 1e7
yrhohigh = 1e9
alp3sig = 0.0

#Code options:
propermotion = 'no'
virialshape = 'yes'

#Overplot true solution (for mock data). If 
#yes, then the true solutions should be passed
#in: ranal,betatrue(ranal),betatruestar(ranal),
#truemass(ranal),trueden(ranal),truedlnrhodlnr(ranal).
overtrue = 'no'

#Radial grid range for Jeans calculation:
rmin = -1.0
rmax = -1.0

#Galaxy properties. Assume here that the baryonic mass
#has the same radial profile as the tracer stars. If this
#is not the case, you should set Mstar_rad and Mstar_prof 
#here. The variables barrad_min, barrad_max and bar_pnts 
#set the radial range and sampling of the baryonic mass model.
Mstar = 0.29e6
Mstar_err = Mstar * 0.25
baryonmass_follows_tracer = 'yes'
barrad_min = 0.0
barrad_max = 10.0
bar_pnts = 250


###########################################################
#Priors

#For surface density fit tracertol = [0,1] sets the spread 
#around the best-fit value from the binulator.
tracertol = 0.5

#Cosmology priors on the coreNFWtides model:
cosmo_cprior = 'no'

#Velocity anisotropy priors:
betr0min = -2
betr0max = 0.0
betnmin = 1.0
betnmax = 10.0
bet0min = -1.0
bet0max = 1.0
betinfmin = -1.0
betinfmax = 1.0

#CoreNFWtides priors:
logM200low = 8.5
logM200high = 10.5
clow = cosmo_cfunc(10.0**logM200high,h)
logclow = np.log10(clow)-0.1
clow = 10.0**logclow
chigh = cosmo_cfunc(10.0**logM200low,h)*1.4
logchigh = np.log10(chigh)+0.2
chigh = 10.0**logchigh
rclow = 1e-2
rchigh = 10.0**0.5
logrclow = np.log10(rclow)
logrchigh = np.log10(rchigh)
sigmlow = 1e-3
sigmhigh = 5.0
nlow = 0.0
nhigh = 1.0
rtlow = 1.0
rthigh = 10.0
logrtlow = np.log10(rtlow)
logrthigh = np.log10(rthigh)
dellow = 3.01
delhigh = 5.0


###########################################################
#Post processing options:

#For calculating J-factors:
get_Juse = get_J
calc_Jfac = 'no'
alpha_Jfac_deg = 0.5 

  





