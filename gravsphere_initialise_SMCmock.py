import numpy as np
from constants import * 
from functions import * 

#This file contains all the code options and choices for 
#running a given model. Throughout, -1 means auto-calculate.

#Data files and output base filename:
whichgal = 'SMCmock_3kpc'
#whichgal = 'SMCmock'
infile = output_base+whichgal+'/'+whichgal
outdirbase = output_base+whichgal+'/'

#Plot ranges and sample points [-1 means auto-calculate]:
rplot_inner = 1e-1
rplot_outer = 10.0
rplot_pnts = 100
y_sigLOSmax = 25
ymin_Sigstar = 1e-4
ymax_Sigstar = 1e-1
yMlow = 1e4
yMhigh = 1e10
yrholow = 1e3
yrhohigh = 1e8
alp3sig = 0.0
sigmlow = 1e-3
sigmhigh = 5.0

#Code options:
propermotion = 'no'
virialshape = 'yes'

#Overplot true solution (for mock data). If 
#yes, then the true solutions should be passed
#in: ranal,betatrue(ranal),betatruestar(ranal),
#truemass(ranal),trueden(ranal),truedlnrhodlnr(ranal).
#Pass a zero array if one of these is not available.
overtrue = 'yes'
data = \
    np.genfromtxt('../Data/SMC_mock/SMC_den.txt',\
                  dtype='f8')
ranal = data[:,0]
trueden = data[:,1]
data = \
    np.genfromtxt('../Data/SMC_mock/SMC_betastar.txt',\
                  dtype='f8')
betatruestar = np.interp(ranal,data[:,0],data[:,1])
betatrue = 2.0*betatruestar/(1.0+betatruestar)
truemass = np.zeros(len(ranal))
truedlnrhodlnr = np.zeros(len(ranal))

#Radial grid range for Jeans calculation:
rmin = -1
rmax = -1

#Galaxy properties. Assume here that the baryonic mass
#has the same radial profile as the tracer stars. If this
#is not the case, you should set Mstar_rad and Mstar_prof 
#here. The variables barrad_min, barrad_max and bar_pnts 
#set the radial range and sampling of the baryonic mass model.
Mstar = 0.0
Mstar_err = 1.0
baryonmass_follows_tracer = 'yes'
barrad_min = 0.0
barrad_max = 10.0
bar_pnts = 250


###########################################################
#Priors

#For surface density fit tracertol = [0,1] sets the spread 
#around the best-fit value from the binulator.
tracertol = 1.0e-3

#Cosmology priors on the coreNFWtides model. mWDM(keV) is
#the mass of a thermal relic; <0 means CDM; sig_c200 is 
#the scatter of c200 in log10 space. If the cosmo_cprior
#is set, then we include a Gaussian spread in M200-c200 in
#the likelihood. Without this, M200-c200 enters only if 
#used to set the priors, below.
cosmo_cprior = 'no'
sig_c200 = 0.1
mWDM = -1
if (mWDM > 0):
    cosmo_cfunc = lambda M200,h : \
        cosmo_cfunc_WDM(M200,h,OmegaM,rhocrit,mWDM)

#Velocity anisotropy priors:
betr0min = -2
betr0max = 1.0
betnmin = 1.0
betnmax = 3.0
bet0min = -0.1
bet0max = 0.1
betinfmin = -0.1
betinfmax = 1.0
#bet0min = -0.01
#bet0max = 0.01
#betinfmin = -0.1
#betinfmax = 1.0

#CoreNFWtides priors:
logM200low = 5.5
logM200high = 11.5
clow = 1.0
chigh = 100.0
rclow = 1e-2
rchigh = 1e2
logrclow = np.log10(rclow)
logrchigh = np.log10(rchigh)
nlow = -1.0
nhigh = 1.0
rtlow = 0.1
rthigh = 10.0
logrtlow = np.log10(rtlow)
logrthigh = np.log10(rthigh)
dellow = 3.01
delhigh = 8.0

if (cosmo_cprior == 'yes'):
    clow = 1.0
    chigh = 100.0

#Priors on central dark mass [set logMcenlow/high very negative
#to switch this off. Mcen is the mass in Msun; acen is the
#scale length in kpc, usually assumed smaller than Rhalf
#to avoid degeneracies with the stellar mass]:
logMcenlow = -4
logMcenhigh = -3
acenlow = 1e-5
acenhigh = 1e-2

#Priors on rotation [Arot defined as:
#vphimean^2 / (2 sigr^2) = Arot(r/Rhalf) which yields linear
#rotation with radius. (Arot = 0.5 means an equal balance of
#rotation and pressure support at Rhalf.)]:
Arotlow = 0.0
Arothigh = 1.0e-12

#Priors on distance [True distance follows as:
#dgal_kpc * drange s.t. we usually want drangelow < 1.0 and
#drangehigh > 1.0]:
dgal_kpc = 60.0
drangelow = 0.99999
drangehigh = 1.00001
    
###########################################################
#Post processing options:

#For calculating D+J-factors:
calc_Jfac = 'no'
alpha_Jfac_deg = 0.5
calc_Dfac = 'no'
alpha_Dfac_deg = 0.5
