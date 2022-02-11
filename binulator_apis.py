import numpy as np
from scipy.integrate import simps as integrator
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

def SegI_api(infile_phot,infile_kin,dgal_kpc):
    #First surface density:
    data_phot = np.genfromtxt(infile_phot,dtype='f8')
    rbin = data_phot[:,0]*dgal_kpc/arcmin
    surfden = data_phot[:,1]
    surfdenerr = (-data_phot[:,2]+data_phot[:,3])/2.0

    #Normalise the surface density and calc.
    #the half light radius from the data:
    Rhalf, Menc_tot = surf_renorm(rbin,surfden)
    surfden = surfden / Menc_tot
    surfdenerr = surfdenerr / Menc_tot
    print('Data Rhalf:', Rhalf)
    
    #Now velocity data:
    data_kin_vs = np.genfromtxt(infile_kin,dtype='f8')
    Rkin = data_kin_vs[:,3]*dgal_kpc/arcmin
    vz = data_kin_vs[:,1]
    vzerr = data_kin_vs[:,2]
    mskin = data_kin_vs[:,6]

    #Remove stars with very large uncertainty:
    ecut = 10.0
    print('Cutting velocity error on %f km/s' % (ecut))
    print('As compared to min/max error: %f, %f km/s' % \
          (np.min(vzerr),np.max(vzerr)))
    Ruset = Rkin[vzerr < ecut]
    vzuset = vz[vzerr < ecut]
    vzerruset = vzerr[vzerr < ecut]
    msuset = mskin[vzerr < ecut]
    
    #Remove non-members:
    pcut = 0.9
    Ruse = Ruset[msuset > pcut]
    vzuse = vzuset[msuset > pcut]
    vzerruse = vzerruset[msuset > pcut]
    msuse = msuset[msuset > pcut]
    
    vsys = np.sum(vzuse*msuse)/np.sum(msuse)
    print('Systemic velocity:',vsys)
    vzuse = vzuse - vsys
    print('Total effective no. of tracers (kinematic):', np.sum(msuse))
    
    return rbin, surfden, surfdenerr, \
        Rhalf, Ruse, vzuse, vzerruse, msuse, vsys

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

    #Cut out stars with large uncertainty:
    ecut = 0.1
    errpmRA = data[:,2]
    esel = errpmRA < ecut
    pmRAt = data[esel,0]
    pmDECt = data[esel,1]
    errpmRAt = data[esel,2]
    errpmDECt = data[esel,3]
    Rt = data[esel,4]
    anglet = data[esel,5]
    esel = errpmDECt < ecut
    pmRA = pmRAt[esel]
    pmDEC = pmDECt[esel]
    errpmRA = errpmRAt[esel]
    errpmDEC = errpmDECt[esel]
    R = Rt[esel]
    angle = anglet[esel]
    mskin = np.zeros(len(pmRA))+1.0
    print('Min/max proper motion error:',\
          np.min(errpmRA),np.max(errpmRA),\
          np.min(errpmDEC),np.max(errpmDEC))
    
    #Calc x-y in kpc:
    y = R*np.cos(angle*np.pi/180.0)
    x = R*np.sin(angle*np.pi/180.0)
    
    #Convert to km/s:
    vx = pmRA*4.74*dgal_kpc
    vxerr = errpmRA*4.74*dgal_kpc
    vy = pmDEC*4.74*dgal_kpc
    vyerr = errpmDEC*4.74*dgal_kpc
    vx = vx - np.sum(vx*mskin)/np.sum(mskin)
    vy = vy - np.sum(vy*mskin)/np.sum(mskin)
    print('Min/max tangential velocity error:',\
          np.min(vxerr),np.max(vxerr),\
          np.min(vyerr),np.max(vyerr))
    
    return x, y, vx, vxerr, vy, vyerr, mskin

def ocen_prebin_api(dgal_kpc,infile_phot,infile_kin,infile_prop):
    #Surface density:
    data_phot = np.genfromtxt(infile_phot,dtype='f8')
    rbin_phot = data_phot[:,0]*dgal_kpc/arcsec
    surfden = data_phot[:,1]
    surfdenerr = (data_phot[:,2]+data_phot[:,3])/2.0

    #Normalise the surface density and calc.
    #the half light radius from the data:
    Rhalf, Menc_tot = surf_renorm(rbin_phot,surfden)
    surfden = surfden / Menc_tot
    surfdenerr = surfdenerr / Menc_tot
    print('Data Rhalf:', Rhalf)
    
    #Now the kinematic data:
    print('Reading RV data ... ',infile_kin)
    data_kin = np.genfromtxt(infile_kin,dtype='f8')
    print('Reading proper motion data ... ',infile_prop)
    data_pm = np.genfromtxt(infile_prop,dtype='f8')

    rbin_kin = data_kin[:,0]*dgal_kpc/arcsec
    sigpmean = data_kin[:,1]
    sigperr = data_kin[:,2]

    rbin_kinp = data_pm[:,0]*dgal_kpc/arcsec
    sigpmr = data_pm[:,1]*dgal_kpc/arcsec*1.0e-3/year*kpc/kms
    sigpmrerr = data_pm[:,2]*dgal_kpc/arcsec*1.0e-3/year*kpc/kms
    sigpmt = data_pm[:,3]*dgal_kpc/arcsec*1.0e-3/year*kpc/kms
    sigpmterr = data_pm[:,4]*dgal_kpc/arcsec*1.0e-3/year*kpc/kms

    #Set up dummy data:
    vzmeanbin = np.zeros(len(rbin_kin))
    vzmeanbinlo = np.zeros(len(rbin_kin))
    vzmeanbinhi = np.zeros(len(rbin_kin))
    vztwobin = sigpmean
    vztwobinlo = sigpmean - sigperr
    vztwobinhi = sigpmean + sigperr
    vzfourbin = np.zeros(len(rbin_kin))
    vzfourbinlo = np.zeros(len(rbin_kin))
    vzfourbinhi = np.zeros(len(rbin_kin))
    backampbin = np.zeros(len(rbin_kin))
    backampbinlo = np.zeros(len(rbin_kin))
    backampbinhi = np.zeros(len(rbin_kin))
    backmeanbin = np.zeros(len(rbin_kin))
    backmeanbinlo = np.zeros(len(rbin_kin))
    backmeanbinhi = np.zeros(len(rbin_kin))
    backsigbin = np.zeros(len(rbin_kin))
    backsigbinlo = np.zeros(len(rbin_kin))
    backsigbinhi = np.zeros(len(rbin_kin))

    vsp1 = 0.0
    vsp1lo = 0.0
    vsp1hi = 0.0
    vsp2 = 0.0
    vsp2lo = 0.0
    vsp2hi = 0.0
    vsp1store = np.zeros(10)
    vsp2store = np.zeros(10)

    vphimeanbin = np.zeros(len(rbin_kinp))
    vphimeanbinlo = np.zeros(len(rbin_kinp))
    vphimeanbinhi = np.zeros(len(rbin_kinp))
    vphitwobin = sigpmt
    vphitwobinlo = sigpmt - sigpmterr
    vphitwobinhi = sigpmt + sigpmterr
    vphifourbin = np.zeros(len(rbin_kinp))
    vphifourbinlo = np.zeros(len(rbin_kinp))
    vphifourbinhi = np.zeros(len(rbin_kinp))
    backptampbin = np.zeros(len(rbin_kinp))
    backptampbinlo = np.zeros(len(rbin_kinp))
    backptampbinhi = np.zeros(len(rbin_kinp))
    backptmeanbin = np.zeros(len(rbin_kinp))
    backptmeanbinlo = np.zeros(len(rbin_kinp))
    backptmeanbinhi = np.zeros(len(rbin_kinp))
    backptsigbin = np.zeros(len(rbin_kinp))
    backptsigbinlo = np.zeros(len(rbin_kinp))
    backptsigbinhi = np.zeros(len(rbin_kinp))
    
    vRmeanbin = np.zeros(len(rbin_kinp))
    vRmeanbinlo = np.zeros(len(rbin_kinp))
    vRmeanbinhi = np.zeros(len(rbin_kinp))
    vRtwobin = sigpmr
    vRtwobinlo = sigpmr - sigpmrerr
    vRtwobinhi = sigpmr + sigpmrerr
    vRfourbin = np.zeros(len(rbin_kinp))
    vRfourbinlo = np.zeros(len(rbin_kinp))
    vRfourbinhi = np.zeros(len(rbin_kinp))
    backpRampbin = np.zeros(len(rbin_kinp))
    backpRampbinlo = np.zeros(len(rbin_kinp))
    backpRampbinhi = np.zeros(len(rbin_kinp))
    backpRmeanbin = np.zeros(len(rbin_kinp))
    backpRmeanbinlo = np.zeros(len(rbin_kinp))
    backpRmeanbinhi = np.zeros(len(rbin_kinp))
    backpRsigbin = np.zeros(len(rbin_kinp))
    backpRsigbinlo = np.zeros(len(rbin_kinp))
    backpRsigbinhi = np.zeros(len(rbin_kinp))

    return rbin_phot, surfden, surfdenerr, Rhalf,\
        rbin_kin,vzmeanbin,vzmeanbinlo,vzmeanbinhi,\
        vztwobin,vztwobinlo,vztwobinhi,\
        vzfourbin,vzfourbinlo,vzfourbinhi,\
        backampbin,backampbinlo,backampbinhi,\
        backmeanbin,backmeanbinlo,backmeanbinhi,\
        backsigbin,backsigbinlo,backsigbinhi,\
        vsp1,vsp1lo,vsp1hi,\
        vsp2,vsp2lo,vsp2hi,\
        vsp1store,\
        vsp2store,\
        rbin_kinp,vphimeanbin,vphimeanbinlo,vphimeanbinhi,\
        vphitwobin,vphitwobinlo,vphitwobinhi,\
        vphifourbin,vphifourbinlo,vphifourbinhi,\
        backptampbin,backptampbinlo,backptampbinhi,\
        backptmeanbin,backptmeanbinlo,backptmeanbinhi,\
        backptsigbin,backptsigbinlo,backptsigbinhi,\
        rbin_kinp,vRmeanbin,vRmeanbinlo,vRmeanbinhi,\
        vRtwobin,vRtwobinlo,vRtwobinhi,\
        vRfourbin,vRfourbinlo,vRfourbinhi,\
        backpRampbin,backpRampbinlo,backpRampbinhi,\
        backpRmeanbin,backpRmeanbinlo,backpRmeanbinhi,\
        backpRsigbin,backpRsigbinlo,backpRsigbinhi

def ocen_api(dgal_kpc,infile_phot,infile_kin,infile_prop):
    #Surface density:
    data_phot = np.genfromtxt(infile_phot,dtype='f8')
    rbin_phot = data_phot[:,0]*dgal_kpc/arcsec
    surfden = data_phot[:,1]
    surfdenerr = (data_phot[:,2]+data_phot[:,3])/2.0
    
    #Normalise the surface density and calc.
    #the half light radius from the data:
    Rhalf, Menc_tot = surf_renorm(rbin_phot,surfden)
    surfden = surfden / Menc_tot
    surfdenerr = surfdenerr / Menc_tot
    print('Data Rhalf:', Rhalf)
    
    #Now the kinematic data:
    print('Reading unbinned RV data ... ',infile_kin)
    data_kin = np.genfromtxt(infile_kin,dtype='f8')
    print('Reading unbinned proper motion data ... ',infile_prop)
    data_pm = np.genfromtxt(infile_prop,dtype='f8')

    Rkin = data_kin[:,0]*dgal_kpc/arcsec
    vz = data_kin[:,1]
    vzerr = data_kin[:,2]
    mskin = np.zeros(len(vzerr))+1.0

    x = data_pm[:,0]*dgal_kpc/arcsec
    y = data_pm[:,1]*dgal_kpc/arcsec
    vx = data_pm[:,2]*dgal_kpc/arcsec*1.0e-3/year*kpc/kms
    vxerr = data_pm[:,3]*dgal_kpc/arcsec*1.0e-3/year*kpc/kms
    vy = data_pm[:,4]*dgal_kpc/arcsec*1.0e-3/year*kpc/kms
    vyerr = data_pm[:,5]*dgal_kpc/arcsec*1.0e-3/year*kpc/kms
    msprop = np.zeros(len(vxerr))+1.0
    
    return rbin_phot, surfden, surfdenerr, Rhalf, \
        Rkin, vz, vzerr, mskin, \
        x, y, vx, vxerr, vy, vyerr, msprop
