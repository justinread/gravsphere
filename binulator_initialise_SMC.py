import numpy as np
from binulator_apis import *
from constants import * 

#Data files and output base filename:
whichgal = 'SMC'
infile_kin = './Data/SMC/Obs_Radial_Vel_raw.dat'
infile_phot = './Data/SMC/Zaritsky_counts.dat'
infile_prop = './Data/SMC/Obs_Prop_Mot_raw.dat'
outfile = output_base+whichgal+'/'+whichgal

#Plot ranges:
xpltmin = 1e-2
xpltmax = 10.0
surfpltmin = 1e-6
surfpltmax = 100
vztwopltmin = 0
vztwopltmax = 35
vzfourpltmin = 1e5
vzfourpltmax = 1e7

#Number of stars per bin [-1 indicates that
#binning was already done elsewhere]:
Nbin = -1
Nbinkin = 50

#Priors for surface density fit. Array values are:
#[M1,M2,M3,a1,a2,a3] where M,a are the Plummer mass
#and scale length. [-1 means use full radial range].
p0in_min = np.array([1e-3,1e-3,-0.5,0.1,0.1,0.1])
p0in_max = np.array([5,5,5,2.5,2.5,2.5])
Rfitmin = -1
Rfitmax = 2.0

#Priors for binulator velocity dispersion calculation. 
#Array values are: [vzmean,alp,bet,backamp,backmean,backsig], 
#where alp is ~the dispersion, bet=[0.1,10] is a shape parameter,
#and "back" is a Gaussian of amplitude "backamp", describing 
#some background. [0 means use full radial range].
p0vin_min = np.array([-20,1.0,1.0,1e-4,-20,40.0])
p0vin_max = np.array([20,40.0,5.0,1.0,20.0,150.0])
vfitmin = 0
vfitmax = 0
Rfitvmin = -1
Rfitvmax = -1

#Convert input data to binulator format (see APIs, above).
#Note that we also calculate Rhalf directly from the data here.
#If use_dataRhalf = 'yes', then we will use this data Rhalf
#instead of the fitted Rhalf. This can be useful if the 
#SB falls off very steeply, which cannot be captured easily
#by the sum over Plummer spheres that binulator/gravsphere
#assumes.
R, surfden, surfdenerr, Rhalf, Rkin, vz, vzerr, mskin = \
    smc_api(infile_phot,infile_kin)
use_dataRhalf = 'yes'

#Propermotions:
propermotion = 'yes'
dgal_kpc = 60.0
if (propermotion == 'yes'):
    Nbinkin_prop = Nbinkin
    x, y, vx, vxerr, vy, vyerr, msprop = \
        smc_prop_api(infile_prop,dgal_kpc)
