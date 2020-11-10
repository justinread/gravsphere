###########################################################
#Binulator
###########################################################

#Python programme to bin discrete dwarf galaxy
#data for Jeans modelling with GravSphere. To
#work with your own data, you should write/add an
#"API" to put it in "binulator" format (see below).
#Binulator will then output all the files that 
#GravSphere needs to run its models.

###########################################################
#Main code:

#Imports & dependencies:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import emcee
from scipy.integrate.quadrature import simps as integrator
from functions import *
from constants import *
from binulator_surffuncs import * 
from binulator_velfuncs import * 
from figures import * 
import sys

#Suppress warning output:
import warnings
warnings.simplefilter("ignore")

#Welcome blurb: 
print('###### BINULATOR VERSION 1.0 ######\n')

###########################################################
#Input data selection here:

#Initialise the data. This bit sets up everything for
#binulator, returning "outfile", rbin,surfden,surfdenerr,
#Nbin, Nbinkin, Rkin, vz, vzerr, mskin etc. (see e.g.
#binulator_initialise_Draco.py for details).

#MW satellites:
#from binulator_initialise_Draco import *
#from binulator_initialise_Fornax import *
#from binulator_initialise_Fornax_tides import *
#from binulator_initialise_SMC import *

#Mocks:
#from binulator_initialise_SMCmock import *
#from binulator_initialise_PlumCoreOm import *
from binulator_initialise_PlumCuspOm import *

#M31 satellites:
#from binulator_initialise_And21 import *

#Code options:
#quicktestSB => set to 'yes' if you want to fit just the 
#SB alone. This is worthwhile to make sure you have a 
#good fit before running the fit to the velocity data
#that takes longer to complete.
quicktestSB = 'no'

#Some output about parameter choices:
print('Doing galaxy:',whichgal)
if (Nbin > 0):
    print('Number of stars per photometric bin:', Nbin)
else:   
    print('Photometric data pre-binned')
if (Nbinkin > 0):
    print('Number of stars per kinematic bin:', Nbinkin)
else:
    print('Kinematic data pre-binned')
print('Priors allow for kurtosis in the range: %.2f < k < %.2f' % \
      (kurt_calc(p0vin_max),kurt_calc(p0vin_min)))
if (Rfitmin > 0):
    print('Fitting R >',Rfitmin)
if (Rfitmax > 0):
    print('Fitting R <',Rfitmax)
if (vfitmin != 0):
    print('Fitting v >',vfitmin)
if (vfitmax != 0):
    print('Fitting v <',vfitmax)    

#Fit surface density:
if (quicktestSB == 'yes'):
    print('### Only fitting the surface brightness profile ###')
print('Fitting the surface density profile:\n')

#Cut back to fit range:
if (Rfitmin > 0):
    Rfit_t = R[R > Rfitmin]
    surfdenfit_t = surfden[R > Rfitmin]
    surfdenerrfit_t = surfdenerr[R > Rfitmin]
else:
    Rfit_t = R
    surfdenfit_t = surfden
    surfdenerrfit_t = surfdenerr
if (Rfitmax > 0):
    Rfit = Rfit_t[Rfit_t < Rfitmax]
    surfdenfit = surfdenfit_t[Rfit_t < Rfitmax]
    surfdenerrfit = surfdenerrfit_t[Rfit_t < Rfitmax]
else:
    Rfit = Rfit_t
    surfdenfit = surfdenfit_t
    surfdenerrfit = surfdenerrfit_t
    
Rplot,surf_int,Rhalf_int,p0best = \
    tracerfit(Rfit,surfdenfit,surfdenerrfit,\
              p0in_min,p0in_max)
print('Fitted Rhalf: %f -%f +%f' % \
    (Rhalf_int[0],Rhalf_int[0]-Rhalf_int[1],\
     Rhalf_int[2]-Rhalf_int[0]))
if (use_dataRhalf == 'yes'):
    Rhalf_int[0] = Rhalf

#Bin the velocity data and calculate the VSPs with uncertainties:
if (quicktestSB == 'no' and Nbinkin > 0):
    print('Fitting the velocity data:\n')
    nsamples = 2500
    
    #Cut back to fit range:
    if (Rfitmin > 0):
        Rf_t = Rkin[Rkin > Rfitmin]
        vzfit_t = vz[Rkin > Rfitmin]
        vzerrfit_t = vzerr[Rkin > Rfitmin]
        msfit_t = mskin[Rkin > Rfitmin]
    else:
        Rf_t = Rkin
        vzfit_t = vz
        vzerrfit_t = vzerr
        msfit_t = mskin
    if (Rfitmax > 0):
        Rf = Rf_t[Rf_t < Rfitmax]
        vzfit = vzfit_t[Rf_t < Rfitmax]
        vzerrfit = vzerrfit_t[Rf_t < Rfitmax]
        msfit = msfit_t[Rf_t < Rfitmax]
    else:
        Rf = Rf_t
        vzfit = vzfit_t
        vzerrfit = vzerrfit_t
        msfit = msfit_t
    
    #Do the fit:
    rbin,vzmeanbin,vzmeanbinlo,vzmeanbinhi,\
        vztwobin,vztwobinlo,vztwobinhi,\
        vzfourbin,vzfourbinlo,vzfourbinhi,\
        backampbin,backampbinlo,backampbinhi,\
        backmeanbin,backmeanbinlo,backmeanbinhi,\
        backsigbin,backsigbinlo,backsigbinhi,\
        vsp1,vsp1lo,vsp1hi,vsp2,vsp2lo,vsp2hi,\
        ranal,vzfourstore,vsp1store,vsp2store = \
            velfit(Rf,vzfit,vzerrfit,msfit,Nbinkin,\
                vfitmin,vfitmax,\
                p0vin_min,p0vin_max,p0best,\
                alpmin,alpmax,nsamples,outfile)
    print('Fitted VSP1: %f+%f-%f' % (vsp1,vsp1hi-vsp1,vsp1-vsp1lo))
    print('Fitted VSP2: %f+%f-%f' % (vsp2,vsp2hi-vsp2,vsp2-vsp2lo))
    

    ###########################################################
    #Store the surface density, velocity dispersions, VSPs,
    #best fit surfden parameters and Rhalf for GravSphere:
    f = open(outfile+'_surfden.txt','w')
    for i in range(len(Rfit)):
        f.write('%f %f %f\n' % \
            (Rfit[i], surfdenfit[i], surfdenerrfit[i]))
    f.close()
    f = open(outfile+'_p0best.txt','w')
    for i in range(len(p0best)):
        f.write('%f\n' % \
            (p0best[i]))
    f.close()
    f = open(outfile+'_Rhalf.txt','w')
    for i in range(len(Rhalf_int)):
        f.write('%f\n' % (Rhalf_int[i]))
    f.close()
    f = open(outfile+'_vel.txt','w')
    for i in range(len(rbin)):
        f.write('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n' % \
            (rbin[i],vzmeanbin[i],vzmeanbinlo[i],vzmeanbinhi[i],\
             vztwobin[i],vztwobinlo[i],vztwobinhi[i],\
             vzfourbin[i],vzfourbinlo[i],vzfourbinhi[i],\
             backampbin[i],backampbinlo[i],backampbinhi[i],\
             backmeanbin[i],backmeanbinlo[i],backmeanbinhi[i],\
             backsigbin[i],backsigbinlo[i],backsigbinhi[i]))
    f.close()
    f = open(outfile+'_vsps.txt','w')
    f.write('%f %f %f\n' % (vsp1,vsp1lo,vsp1hi))
    f.write('%f %f %f\n' % (vsp2,vsp2lo,vsp2hi))
    f.close()
    f = open(outfile+'_vsp1full.txt','w')
    for i in range(len(vsp1store)):
        f.write('%f\n' % (vsp1store[i]))
    f.close()
    f = open(outfile+'_vsp2full.txt','w')
    for i in range(len(vsp2store)):
        f.write('%f\n' % (vsp2store[i]))
    f.close()


###########################################################
#Plots:

##### Surface density fit #####
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(mylinewidth)
ax.minorticks_on()
ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
plt.xticks(fontsize=myfontsize)
plt.yticks(fontsize=myfontsize)

plt.xlabel(r'Radius [kpc]',fontsize=myfontsize)
plt.ylabel(r'Surface density [${\rm kpc}^{-2}$]',fontsize=myfontsize)

plt.loglog()

plt.errorbar(R,surfden,surfdenerr,color='black',ecolor='black',\
             fmt='o',\
             linewidth=mylinewidth,marker='o',markersize=5,\
             alpha=0.5)
plt.fill_between(Rplot,surf_int[5,:],surf_int[6,:],\
            facecolor='blue',alpha=0.25,\
            edgecolor='none')
plt.fill_between(Rplot,surf_int[3,:],surf_int[4,:],\
            facecolor='blue',alpha=0.33,\
            edgecolor='none')
plt.fill_between(Rplot,surf_int[1,:],surf_int[2,:],\
            facecolor='blue',alpha=0.66,\
            edgecolor='none')
            
#Best fit model:
surfbest = threeplumsurf(Rplot,p0best[0],p0best[1],p0best[2],\
                         p0best[3],p0best[4],p0best[5])
plt.plot(Rplot,surfbest,color='red',linewidth=2)

plt.axvline(x=Rhalf_int[0],color='blue',alpha=0.5,\
            linewidth=mylinewidth)

plt.xlim([xpltmin,xpltmax])
plt.ylim([surfpltmin,surfpltmax])

plt.savefig(outfile+'_surfden.pdf',bbox_inches='tight')

##### Tracer 3D density profile (best fit) #####
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(mylinewidth)
ax.minorticks_on()
ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
plt.xticks(fontsize=myfontsize)
plt.yticks(fontsize=myfontsize)

plt.xlabel(r'Radius [kpc]',fontsize=myfontsize)
plt.ylabel(r'Density [${\rm kpc}^{-3}$]',fontsize=myfontsize)

plt.loglog()

denbest = threeplumden(Rplot,p0best[0],p0best[1],p0best[2],\
                         p0best[3],p0best[4],p0best[5])
plt.plot(Rplot,denbest,color='red',linewidth=2)

plt.axvline(x=Rhalf_int[0],color='blue',alpha=0.5,\
            linewidth=mylinewidth)

plt.xlim([xpltmin,xpltmax])
plt.ylim([surfpltmin,surfpltmax])

plt.savefig(outfile+'_tracerden.pdf',bbox_inches='tight')

##### Tracer 3D mass profile (best fit) #####
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(mylinewidth)
ax.minorticks_on()
ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
plt.xticks(fontsize=myfontsize)
plt.yticks(fontsize=myfontsize)

plt.xlabel(r'Radius [kpc]',fontsize=myfontsize)
plt.ylabel(r'Normalised Mass',fontsize=myfontsize)

plt.loglog()

massbest = threeplummass(Rplot,p0best[0],p0best[1],p0best[2],\
                         p0best[3],p0best[4],p0best[5])
plt.plot(Rplot,massbest,color='red',linewidth=2)

plt.axvline(x=Rhalf_int[0],color='blue',alpha=0.5,\
            linewidth=mylinewidth)

plt.xlim([xpltmin,xpltmax])
plt.ylim([1e-4,10.0])

plt.savefig(outfile+'_tracermass.pdf',bbox_inches='tight')

##### Plot velocity bin fits ##### 
if (quicktestSB == 'no'):
    ##### Velocity dispersion profile ##### 
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    ax.minorticks_on()
    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)

    plt.xlabel(r'Radius [kpc]',fontsize=myfontsize)
    plt.ylabel(r'Velocity dispersion [${\rm km\,s}^{-1}$]',\
           fontsize=myfontsize)

    plt.semilogx()

    plt.errorbar(rbin,vztwobin,\
             yerr=[vztwobin-vztwobinlo,\
                   vztwobinhi-vztwobin],\
             color='black',ecolor='black',\
             fmt='o',\
             linewidth=mylinewidth,marker='o',\
             markersize=5,alpha=0.5)

    plt.axvline(x=Rhalf_int[0],color='blue',alpha=0.5,\
            linewidth=mylinewidth)

    plt.xlim([xpltmin,xpltmax])
    plt.ylim([vztwopltmin,vztwopltmax])

    plt.savefig(outfile+'_vztwo.pdf',bbox_inches='tight')

    ##### Fourth velocity moment ##### 
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    ax.minorticks_on()
    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)

    plt.xlabel(r'Radius [kpc]',fontsize=myfontsize)
    plt.ylabel(r'Fourth velocity moment [${\rm km}^4\,{\rm s}^{-4}$]',\
            fontsize=myfontsize)

    plt.loglog()

    plt.errorbar(rbin,vzfourbin,\
             yerr=[vzfourbin-vzfourbinlo,\
                   vzfourbinhi-vzfourbin],\
             color='black',ecolor='black',\
             fmt='o',\
             linewidth=mylinewidth,marker='o',\
             markersize=5,alpha=0.5)

    #Plot individual fits + falloff:
    for i in range(nsamples):
        plt.plot(ranal,vzfourstore[i,:],linewidth=1,alpha=0.25)

    plt.axvline(x=Rhalf_int[0],color='blue',alpha=0.5,\
            linewidth=mylinewidth)

    plt.xlim([xpltmin,xpltmax])
    plt.ylim([vzfourpltmin,vzfourpltmax])

    plt.savefig(outfile+'_vzfour.pdf',bbox_inches='tight')

    ###### VSP1 and VSP2 histograms #####
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)

    nbin = 25
    n, bins, patches = plt.hist(vsp1store,bins=nbin,\
        range=(np.min(vsp1store),\
               np.max(vsp1store)),\
        facecolor='b', \
        histtype='bar',alpha=0.5, \
        label='vs_1')
 
    plt.xlabel(r'$v_{s1}\,[{\rm km}^4\,{\rm s}^{-4}]$',\
           fontsize=myfontsize)
    plt.ylabel(r'$N$',fontsize=myfontsize)
    plt.ylim([0,np.max(n)])
    plt.savefig(outfile+'output_vsp1.pdf',bbox_inches='tight')

    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)

    nbin = 25
    n, bins, patches = plt.hist(vsp2store,bins=nbin,\
        range=(np.min(vsp2store),\
               np.max(vsp2store)),\
        facecolor='b', \
        histtype='bar',alpha=0.5, \
        label='vs_2')
 
    plt.xlabel(r'$v_{s2}\,[{\rm km}^4\,{\rm s}^{-4}]$',\
           fontsize=myfontsize)
    plt.ylabel(r'$N$',fontsize=myfontsize)
    plt.ylim([0,np.max(n)])
    plt.savefig(outfile+'output_vsp2.pdf',bbox_inches='tight')


###########################################################
#Exit:
print('\nThank you for using the Binulator! Have a nice day.\n')
