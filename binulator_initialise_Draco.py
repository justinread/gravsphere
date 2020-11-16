import numpy as np
from binulator_apis import *
from constants import * 

#Data files and output base filename:
whichgal = 'Draco'
infile_kin = './Data/Walker_dwarfs/dra_justin1_spec.dat'
infile_phot = './Data/Walker_dwarfs/dra_justin1_phot.dat'
infile_jardel = './Data/Walker_dwarfs/jardel_virus.txt'
outfile = output_base+whichgal+'/'+whichgal

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
Nbin = 30
Nbinkin = 30

#Priors for surface density fit. Array values are:
#[M1,M2,M3,a1,a2,a3] where M,a are the Plummer mass
#and scale length. [-1 means use full radial range].
p0in_min = np.array([1e-4,1e-4,1e-4,0.01,0.01,0.01])
p0in_max = np.array([1e2,1e2,1e2,1.0,1.0,1.0])
Rfitmin = -1
Rfitmax = -1

#Priors for binulator velocity dispersion calculation. 
#Array values are: [vzmean,alp,bet,backamp,backmean,backsig], 
#where alp is ~the dispersion, bet=[0.1,10] is a shape parameter,
#and "back" is a Gaussian of amplitude "backamp", describing 
#some background. [0 means use full radial range].
p0vin_min = np.array([-0.1,1.0,1.0,1e-5,-300.0,25.0])
p0vin_max = np.array([0.1,25.0,5.0,1e-4,300.0,300.0])
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
dgal_kpc = 76.0
R, surfden, surfdenerr, Rhalf, \
    Rkin, vz, vzerr, mskin, vsys = \
    walker_api(infile_phot,infile_kin,dgal_kpc,Nbin)
use_dataRhalf = 'no'

#Include the Virus-W inner data from:
#https://ui.adsabs.harvard.edu/abs/2013ApJ...763...91J/abstract
include_jardel = 'no'
if (include_jardel == 'yes'):
    print('Including Virus-W data from Jardel et al. 2013')
    data_jar = np.genfromtxt(infile_jardel,dtype='f8')
    RA_DRA_deg = 260.0516667  #McConnachie rev.
    DEC_DRA_deg = 57.9152778  #
    vsys_DRA_walker = -291.0  #

    rjar = np.sqrt(((data_jar[:,0] - RA_DRA_deg)*\
                np.cos(DEC_DRA_deg/360.0*2.0*np.pi))**2.0+\
                (data_jar[:,1] - DEC_DRA_deg)**2.0)/\
                360.0*2.0*np.pi*dgal_kpc
    vzjar = data_jar[:,2]
    vzjarerr = data_jar[:,3]
    msjar = np.zeros(len(vzjar))+1.0    
    print('Jardel systemic velocity:',\
          np.sum(vzjar*msjar)/np.sum(msjar))
    vzjar = vzjar - vsys

    Rkin = np.concatenate((rjar,Rkin))
    vz = np.concatenate((vzjar,vz))
    vzerr = np.concatenate((vzjarerr,vzerr))
    mskin = np.concatenate((msjar,mskin))
    print('Updated effective no. of tracers (kinematic):',\
          np.sum(mskin))
