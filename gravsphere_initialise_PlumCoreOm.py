import numpy as np
from constants import * 
from functions import * 

#This file contains all the code options and choices for 
#running a given model. Throughout, -1 means auto-calculate.

#Data files and output base filename:
nstars = 1000
if (nstars == 1000):
    whichgal = 'PlumCoreOm'
elif (nstars == 100):
    whichgal = 'PlumCoreOm100'
else:
    whichgal = 'PlumCoreOm10000'
infile = output_base+'GCmock/'+whichgal+'/'+whichgal
outdirbase = output_base+'GCmock/'+whichgal+'/'

#Plot ranges and sample points [-1 means auto-calculate]:
rplot_inner = 1e-2
rplot_outer = 5.0
rplot_pnts = 50
y_sigLOSmax = 25
ymin_Sigstar = 1e-6
ymax_Sigstar = 100
yMlow = 1e4
yMhigh = 1e10
yrholow = 1e5
yrhohigh = 1e9
alp3sig = 0.2
sigmlow = 1e-3
sigmhigh = 5.0

#Code options:
propermotion = 'yes'
virialshape = 'yes'

#Overplot true solution (for mock data). If 
#yes, then the true solutions should be passed
#in: ranal,betatrue(ranal),betatruestar(ranal),
#truemass(ranal),trueden(ranal),truedlnrhodlnr(ranal).
overtrue = 'yes'
rho0,r0,alp,bet,gam,rstar,ra = \
    np.array([400./1000. * 1000.**3.,1.0,1.0,\
              3.0,0.0,25./100.,100./100.*25./100.])
ranal = np.logspace(-3,1,np.int(250))
betatrue = ranal**2./(ranal**2. + ra**2.)
betatruestar = betatrue/(2.0-betatrue)
truemass = alpbetgammass(ranal,rho0,r0,alp,bet,gam)
trueden = alpbetgamden(ranal,rho0,r0,alp,bet,gam)
truedlnrhodlnr = alpbetgamdlnrhodlnr(ranal,rho0,r0,alp,bet,gam)

#Radial grid range for Jeans calculation:
rmin = -1
rmax = -1

#Galaxy properties. Assume here that the baryonic mass
#has the same radial profile as the tracer stars. If this
#is not the case, you should set Mstar_rad and Mstar_prof 
#here. The variables barrad_min, barrad_max and bar_pnts 
#set the radial range and sampling of the baryonic mass model.
Mstar = -1.0
Mstar_err = 1.0
baryonmass_follows_tracer = 'yes'
barrad_min = 0.0
barrad_max = 10.0
bar_pnts = 250


###########################################################
#Priors

#For surface density fit tracertol = [0,1] sets the spread 
#around the best-fit value from the binulator.
tracertol = 0.1

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
betr0max = 0.0
betnmin = 1.0
betnmax = 3.0
bet0min = -0.01
bet0max = 0.01
betinfmin = -0.1
betinfmax = 1.0

#CoreNFWtides priors:
logM200low = 7.5
logM200high = 11.5
clow = 1.0
chigh = 100.0
rclow = 1e-2
rchigh = 10.0
logrclow = np.log10(rclow)
logrchigh = np.log10(rchigh)
nlow = 0.0
nhigh = 1.0
rtlow = 1.0
rthigh = 20.0
logrtlow = np.log10(rtlow)
logrthigh = np.log10(rthigh)
dellow = 3.01
delhigh = 5.0

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
dgal_kpc = 100.0
drangelow = 0.99999
drangehigh = 1.00001
    
###########################################################
#Post processing options:

#For calculating D+J-factors:
calc_Jfac = 'no'
alpha_Jfac_deg = 0.5 
calc_Dfac = 'no'
alpha_Dfac_deg = 0.5






