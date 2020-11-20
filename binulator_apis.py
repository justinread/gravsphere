import numpy as np
from scipy.integrate.quadrature import simps as integrator
from constants import *
from functions import * 

###########################################################
#APIs for loading in data. These should convert input data
#to the "binulator" format. This is:
#   1. Surface density as: R[kpc]; surfden[N kpc^-2], 
#      normalised such that its integral = 1; and surfdenerr,  
#      a 68% confidence interval, assumed to be Gaussian 
#      symmetric. 
#   2. Radial positions (Rkin[kpc]) and line-of-sight velocities
#      for tracers (vz[km/s]) with errors (vzerr[km/s]).
#      Each star can have a membership probability, mskin=[0,1].
#      Velocity errors are assumed to be Gaussian symmetric.
#      The mean velocity should be subtracted. (This is
#      not vital per-se, but makes it easier to set more general
#      priors on the the velocity pdf fits.) 
#
#Some APIs additionally calculate Rhalf directly from the data, 
#though this is not required. In this case, setting
#use_dataRhalf = 'yes' in the initalisation script will use 
#this data Rhalf instead of the fitted Rhalf. This can be useful if 
#the SB falls off very steeply, which cannot be captured easily
#by the sum over Plummer spheres that binulator/gravsphere
#assumes.

def walker_api(infile_phot,infile_kin,dgal_kpc,Nbin):    
    #First surface density:
    data_phot = np.genfromtxt(infile_phot,dtype='f8')
    R = data_phot[:,4]*dgal_kpc/arcmin
    ms = data_phot[:,10]
    print('Total effective no. of tracers (photometric):', np.sum(ms))
    rbin, surfden, surfdenerr, Rhalf = \
        binthedata(R,ms,Nbin)
    print('Data Rhalf:', Rhalf)

    #Now velocity data:
    data_kin_vs = np.genfromtxt(infile_kin,dtype='f8')
    gotvz = np.where(data_kin_vs[:,6])[0]
    Rkin = data_kin_vs[gotvz,4]*dgal_kpc/arcmin
    vz = data_kin_vs[gotvz,6]
    vzerr = data_kin_vs[gotvz,7]
    mskin = data_kin_vs[gotvz,15]
    vsys = np.sum(vz*mskin)/np.sum(mskin)
    print('Systemic velocity:',vsys)
    vz = vz - vsys
    print('Total effective no. of tracers (kinematic):', np.sum(mskin))

    return rbin, surfden, surfdenerr, \
        Rhalf, Rkin, vz, vzerr, mskin, vsys

def collins_api(infile_phot,infile_kin,Nbin):    
    #First surface density:
    data_phot = np.genfromtxt(infile_phot,dtype='f8')
    R = data_phot[:,0]/1000.0
    ms = data_phot[:,1]
    print('Total effective no. of tracers (photometric):', np.sum(ms))
    rbin, surfden, surfdenerr, Rhalf = \
        binthedata(R,ms,Nbin)
    print('Data Rhalf:', Rhalf)

    #Now velocity data:
    data_kin_vs = np.genfromtxt(infile_kin,dtype='f8')
    Rkin = data_kin_vs[:,0]/1000.0
    vz = data_kin_vs[:,1]
    vzerr = data_kin_vs[:,2]
    mskin = data_kin_vs[:,3]
    vz = vz - np.sum(vz*mskin)/np.sum(mskin)
    print('Total effective no. of tracers (kinematic):', np.sum(mskin))

    return rbin, surfden, surfdenerr, Rhalf, Rkin, vz, vzerr, mskin
    
def smc_api(infile_phot,infile_kin):    
    #First surface density:
    data_phot = np.genfromtxt(infile_phot,dtype='f8')
    
    #This is already binned:
    rbin = data_phot[:,2]
    surfden = data_phot[:,0]
    surfdenerr = data_phot[:,1]

    #Normalise the surface density and calc.
    #the half light radius from the data:
    Rhalf, Menc_tot = surf_renorm(rbin,surfden)
    surfden = surfden / Menc_tot
    surfdenerr = surfdenerr / Menc_tot
    print('Data Rhalf:', Rhalf)

    #Now velocity data:
    data_kin = np.genfromtxt(infile_kin,dtype='f8')
    Rkin = data_kin[:,2]
    vz = data_kin[:,0]
    vzerr = data_kin[:,1]
    mskin = np.zeros(len(Rkin))+1.0
    vz = vz - np.sum(vz*mskin)/np.sum(mskin)
    print('Total effective no. of tracers (kinematic):', np.sum(mskin))

    return rbin, surfden, surfdenerr, Rhalf, Rkin, vz, vzerr, mskin
    
def gc_api(data_file_phot,data_file_kin,Nbin):
    #Read in the data:
    data = np.genfromtxt(data_file_phot,dtype='f8')

    #Surface density:
    R = np.sqrt(data[:,0]**2.0 + data[:,1]**2.0)
    ms = np.zeros(len(R)) + 1.0
    print('Total effective no. of tracers (photometric):', np.sum(ms))
    rbin, surfden, surfdenerr, Rhalf = \
        binthedata(R,ms,Nbin)
    print('Data Rhalf:', Rhalf)
    
    #Velocity data:
    data = np.genfromtxt(data_file_kin,dtype='f8')
    Rkin = np.sqrt(data[:,0]**2.0 + data[:,1]**2.0)
    vz = data[:,5]
    vzerr = np.zeros(len(vz))+2.0
    mskin = np.zeros(len(Rkin)) + 1.0
    vz = vz - np.sum(vz*mskin)/np.sum(mskin)
    print('Total effective no. of tracers (kinematic):', np.sum(mskin))

    return rbin, surfden, surfdenerr, Rhalf, Rkin, vz, vzerr, mskin

def gc_prop_api(data_file_kin):
    #Propermotion data:
    data = np.genfromtxt(data_file_kin,dtype='f8')
    x = data[:,0]
    y = data[:,1]
    vx = data[:,3]
    vxerr = np.zeros(len(vx))+2.0
    vy = data[:,4]
    vyerr = np.zeros(len(vy))+2.0
    mskin = np.zeros(len(x))+1.0
    vx = vx - np.sum(vx*mskin)/np.sum(mskin)
    vy = vy - np.sum(vy*mskin)/np.sum(mskin)
    
    return x, y, vx, vxerr, vy, vyerr, mskin

def smc_prop_api(data_file_kin,dgal_kpc):
    #Propermotion data:
    data = np.genfromtxt(data_file_kin,dtype='f8')
    pmRA = data[:,0]
    pmDEC = data[:,1]
    errpmRA = data[:,2]
    errpmDEC = data[:,3]
    R = data[:,4]
    angle = data[:,5]
    mskin = np.zeros(len(pmRA))+1.0

    #Calc x-y in kpc:
    y = R*np.cos(angle*np.pi/180.0)
    x = R*np.sin(angle*np.pi/180.0)
    
    #Convert to km/s:
    vx = pmRA*4.74*dgal_kpc
    vxerr = errpmRA*4.74*dgal_kpc
    vy = pmDEC*4.74*dgal_kpc
    vyerr = errpmDEC*4.74*dgal_kpc
    
    return x, y, vx, vxerr, vy, vyerr, mskin
