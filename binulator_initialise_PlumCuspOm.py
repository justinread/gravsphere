import numpy as np
from binulator_apis import *
from constants import * 

#Data files and output base filename:
nstars = 10000
if (nstars == 1000):
    whichgal = 'PlumCuspOm'
    data_file_kin = \
        './Data/GC_mock/PlumCuspOm/gs010_bs050_rcrs010_rarc100_cusp_0064mpc3_df_1000_0_err.dat'
    #Number of stars per bin [-1 indicates that
    #binning was already done elsewhere]:
    Nbinkin = 30
elif (nstars == 100):
    whichgal = 'PlumCuspOm100'
    data_file_kin = \
        './Data/GC_mock/PlumCuspOm/gs010_bs050_rcrs010_rarc100_cusp_0064mpc3_df_100_0_err.dat'
    #Number of stars per bin [-1 indicates that
    #binning was already done elsewhere]:
    Nbinkin = 25 
else:
    whichgal = 'PlumCuspOm10000'
    data_file_kin = \
        './Data/GC_mock/PlumCuspOm/gs010_bs050_rcrs010_rarc100_cusp_0064mpc3_df_10000_0_err.dat'
    #Number of stars per bin [-1 indicates that
    #binning was already done elsewhere]:
    Nbinkin = 100
Nbin = 100
data_file_phot = \
    './Data/GC_mock/PlumCuspOm/gs010_bs050_rcrs010_rarc100_cusp_0064mpc3_df_10000_0_err.dat'
outfile = output_base+'GCmock/'+whichgal+'/'+whichgal

#Plot ranges:
xpltmin = 1e-2
xpltmax = 10.0
surfpltmin = 1e-6
surfpltmax = 100
vztwopltmin = 0
vztwopltmax = 15
vzfourpltmin = 1e2
vzfourpltmax = 1e6

#Priors for surface density fit. Array values are:
#[M1,M2,M3,a1,a2,a3] where M,a are the Plummer mass
#and scale length. [-1 means use full radial range].
p0in_min = np.array([1e-4,1e-4,1e-4,0.01,0.2,0.2])
p0in_max = np.array([1e2,1e2,1e2,1.0,1.0,1.0])
Rfitmin = -1
Rfitmax = -1

#Priors for binulator velocity dispersion calculation. 
#Array values are: [vzmean,alp,bet,backamp,backmean,backsig], 
#where alp is ~the dispersion, bet=[0.1,10] is a shape parameter,
#and "back" is a Gaussian of amplitude "backamp", describing 
#some background. [0 means use full radial range].
p0vin_min = np.array([-0.1,1.0,0.5,1e-5,-300.0,25.0])
p0vin_max = np.array([0.1,25.0,2.0,1e-4,300.0,300.0])
vfitmin = 0
vfitmax = 0

#Priors for VSP calculation. This is the range of
#powerlaw index for the assumed fall-off of <vlos^4>
#beyond the outermost data point.
alpmin = 0.0
alpmax = 1.0

#Convert input data to binulator format (see APIs, above).
#Note that we also calculate Rhalf directly from the data here.
#If use_dataRhalf = 'yes', then we will use this data Rhalf
#instead of the fitted Rhalf. This can be useful if the 
#SB falls off very steeply, which cannot be captured easily
#by the sum over Plummer spheres that binulator/gravsphere
#assumes.
R, surfden, surfdenerr, Rhalf, Rkin, vz, vzerr, mskin = \
    gc_api(data_file_phot,data_file_kin,Nbin)
use_dataRhalf = 'no'

#Calculate true VSPs and output for comparison:
rho0s,r0s,alps,bets,gams = np.array([1.0,0.1,2,5,0.1])
rho0,r0,alp,bet,gam,rstar,ra = \
    np.array([64./1000. * 1000.**3.,1.0,1.0,\
              3.0,1.0,10./100.,100./100.*10./100.])
vsp1, vsp2, zeta_A, zeta_B = \
    alpbetgamvsp(rho0s,r0s,alps,bets,gams,\
                 rho0,r0,alp,bet,gam,ra)
print('True VSPs:',vsp1,vsp2)

#Richardson / Fairbairn VSPs:
vsp1_RF, vsp2_RF = richfair_vsp(vz,Rkin,mskin)
print('Richardson/Fairbairn VSPs:',vsp1_RF,vsp2_RF)
