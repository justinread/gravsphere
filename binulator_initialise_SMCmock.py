import numpy as np
from binulator_apis import *
from constants import * 

#Data files and output base filename:
whichgal = 'SMCmock'
infile_kin = './Data/SMC_mock/New_Radial_Vel_3sig.dat'
infile_phot = './Data/SMC_mock/Sims_Star_counts_full.dat'
outfile = './Output/'+whichgal+'/'+whichgal

#Plot ranges:
xpltmin = 1e-2
xpltmax = 10.0
surfpltmin = 1e-6
surfpltmax = 100
vztwopltmin = 0
vztwopltmax = 15
vzfourpltmin = 1e2
vzfourpltmax = 1e6

#Number of stars per bin [-1 indicates that
#binning was already done elsewhere]:
Nbin = -1
Nbinkin = 150

#Priors for surface density fit. Array values are:
#[M1,M2,M3,a1,a2,a3] where M,a are the Plummer mass
#and scale length. [-1 means use full radial range].
p0in_min = np.array([1e-4,1e-4,1e-4,0.01,0.01,0.01])
p0in_max = np.array([1e2,1e2,1e2,5.0,5.0,10.0])
Rfitmin = -1
Rfitmax = -1

#Priors for binulator velocity dispersion calculation. 
#Array values are: [vzmean,alp,bet,backamp,backmean,backsig], 
#where alp is ~the dispersion, bet=[0.1,10] is a shape parameter,
#and "back" is a Gaussian of amplitude "backamp", describing 
#some background. [0 means use full radial range].
p0vin_min = np.array([-10,1.0,1.0,1e-4,-10,5.0])
p0vin_max = np.array([10,10.0,3.0,1.0,10.0,15.0])
vfitmin = -20.0
vfitmax = 20.0

#Priors for VSP calculation. This is the range of
#powerlaw index for the assumed fall-off of <vlos^4>
#beyond the outermost data point.
alpmin = 0.0
alpmax = 3.0

#Convert input data to binulator format (see APIs, above).
#Note that we also calculate Rhalf directly from the data here.
#This is superseded later by the calculation from the fit but
#we include it here as a sanity check.
R, surfden, surfdenerr, Rhalf, Rkin, vz, vzerr, mskin = \
    smc_api(infile_phot,infile_kin)
