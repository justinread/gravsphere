import numpy as np
from binulator_apis import *
from constants import * 

#Data files and output base filename:
whichgal = 'PlumCoreOm'
data_file = \
    './Data/GC_mock/PlumCoreOm/gs010_bs050_rcrs025_rarc100_core_0400mpc3_df_1000_0_err.dat'
outfile = './Output/GCmock/'+whichgal+'/'+whichgal

#Plot ranges:
xpltmin = 1e-2
xpltmax = 10.0
surfpltmin = 1e-6
surfpltmax = 100
vztwopltmin = 0
vztwopltmax = 25
vzfourpltmin = 1e3
vzfourpltmax = 1e7

#Number of stars per bin [-1 indicates that
#binning was already done elsewhere]:
Nbin = 15
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
p0vin_min = np.array([-0.1,1.0,1.0,1e-5,-300.0,30.0])
p0vin_max = np.array([0.1,30.0,5.0,1e-4,300.0,300.0])
vfitmin = 0
vfitmax = 0

#Priors for VSP calculation. This is the range of
#powerlaw index for the assumed fall-off of <vlos^4>
#beyond the outermost data point.
alpmin = 0.0
alpmax = 1.0

#Convert input data to binulator format (see APIs, above).
#Note that we also calculate Rhalf directly from the data here.
#This is superseded later by the calculation from the fit but
#we include it here as a sanity check.
R, surfden, surfdenerr, Rhalf, Rkin, vz, vzerr, mskin = \
    gc_api(data_file,Nbin)
    
#Calculate true VSPs and output for comparison:
rho0s,r0s,alps,bets,gams = np.array([1.0,0.25,2,5,0.1])
rho0,r0,alp,bet,gam,rstar,ra = \
    np.array([400./1000. * 1000.**3.,1.0,1.0,\
              3.0,0.0,25./100.,100./100.*25./100.])
vsp1, vsp2, zeta_A, zeta_B = \
    alpbetgamvsp(rho0s,r0s,alps,bets,gams,\
                 rho0,r0,alp,bet,gam,ra)
print('True VSPs:',vsp1,vsp2)

#Richardson / Fairbairn VSPs:
vsp1_RF, vsp2_RF = richfair_vsp(vz,Rkin,mskin)
print('Richardson/Fairbairn VSPs:',vsp1_RF,vsp2_RF)


