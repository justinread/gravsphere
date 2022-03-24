import numpy as np
from binulator_apis import *
from constants import * 

#Run on multiproc:
nprocs = 10

#Data files and output base filename:
whichgal = 'Ocen'
infile_kin = '../Data/Ocen/Ocen_NoyoloVDV_RVs_N.txt'
infile_phot = '../Data/Ocen/Ocen_surfden.txt'
infile_prop = '../Data/Ocen/Ocen_HST_Gaia_PMs.txt'
dgal_kpc = 5.5
outfile = output_base+whichgal+'/dgal_%.1f/' % (dgal_kpc)+whichgal
print('Using Omega cen distance: %.1f kpc' % (dgal_kpc))

#Plot ranges:
xpltmin = 1e-6
xpltmax = 5.0
surfpltmin = 1e-3
surfpltmax = 1e5
vztwopltmin = 0
vztwopltmax = 30
vzfourpltmin = 1e5
vzfourpltmax = 1e7

#Number of stars per bin [-1 indicates that
#binning was already done elsewhere]:
Nbin = -1
Nbinkin = -1

#Priors for surface density fit. Array values are:
#[M1,M2,M3,a1,a2,a3] where M,a are the Plummer mass
#and scale length. [-1 means use full radial range].
p0in_min = np.array([0,0,0,0.001,0.001,0.001])
p0in_max = np.array([9,9,9,0.05,0.05,0.05])
Rfitmin = -1
Rfitmax = -1

#Priors for binulator velocity dispersion calculation. 
#Array values are: [vzmean,alp,bet,backamp,backmean,backsig], 
#where alp is ~the dispersion, bet=[0.1,10] is a shape parameter,
#and "back" is a Gaussian of amplitude "backamp", describing 
#some background. [0 means use full radial range].
p0vin_min = np.array([-50.0,10.0,1.0,1e-4,-250.0,60.0])
p0vin_max = np.array([50.0,60.0,5.0,1.0,250.0,150.0])
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

#Ocen data is prebinned so most of these variables will be
#just dummy values so that binulator outputs the right
#files for gravsphere to read.
R, surfden, surfdenerr, Rhalf,\
    rbin,vzmeanbin,vzmeanbinlo,vzmeanbinhi,\
    vztwobin,vztwobinlo,vztwobinhi,\
    vzfourbin,vzfourbinlo,vzfourbinhi,\
    backampbin,backampbinlo,backampbinhi,\
    backmeanbin,backmeanbinlo,backmeanbinhi,\
    backsigbin,backsigbinlo,backsigbinhi,\
    vsp1,vsp1lo,vsp1hi,\
    vsp2,vsp2lo,vsp2hi,\
    vsp1store,\
    vsp2store,\
    rbinpt,vphimeanbin,vphimeanbinlo,vphimeanbinhi,\
    vphitwobin,vphitwobinlo,vphitwobinhi,\
    vphifourbin,vphifourbinlo,vphifourbinhi,\
    backptampbin,backptampbinlo,backptampbinhi,\
    backptmeanbin,backptmeanbinlo,backptmeanbinhi,\
    backptsigbin,backptsigbinlo,backptsigbinhi,\
    rbinpR,vRmeanbin,vRmeanbinlo,vRmeanbinhi,\
    vRtwobin,vRtwobinlo,vRtwobinhi,\
    vRfourbin,vRfourbinlo,vRfourbinhi,\
    backpRampbin,backpRampbinlo,backpRampbinhi,\
    backpRmeanbin,backpRmeanbinlo,backpRmeanbinhi,\
    backpRsigbin,backpRsigbinlo,backpRsigbinhi = \
        ocen_prebin_api(dgal_kpc,infile_phot,\
                        infile_kin,infile_prop)
use_dataRhalf = 'no'
propermotion = 'yes'
