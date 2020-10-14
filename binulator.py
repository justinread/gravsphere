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

#Welcome blurb: 
print('###### BINULATOR VERSION 1.0 ######\n')

#Initialise the data. This bit sets up everything for
#binulator, returning "outfile", rbin,surfden,surfdenerr,
#Nbin, Nbinkin, Rkin, vz, vzerr, mskin etc. (see e.g.
#binulator_initialise_Draco.py for details).

#from binulator_initialise_Draco import *
#from binulator_initialise_SMCmock import *
#from binulator_initialise_PlumCoreOm import *
from binulator_initialise_PlumCuspOm import *

#Some output about parameter choices:
if (Nbin > 0):
    print('Number of stars per photometric bin:', Nbin)
else:   
    print('Photometric data pre-binned')
if (Nbinkin > 0):
    print('Number of stars per kinematic bin:', Nbinkin)
else:
    print('Kinematic data pre-binned')

#Fit surface density:
print('Fitting the surface density profile:\n')
Rplot,surf_int,Rhalf_int,p0best = \
    tracerfit(R,surfden,surfdenerr,\
              Rfitmin,Rfitmax,\
              p0in_min,p0in_max)
print('Fitted Rhalf: %f -%f +%f' % \
    (Rhalf_int[0],Rhalf_int[0]-Rhalf_int[1],Rhalf_int[2]-Rhalf_int[0]))

#Bin the velocity data and calculate the VSPs with uncertainties:
print('Fitting the velocity data:\n')
nsamples = 500
rbin,vzmeanbin,vzmeanbinerr,vztwobin,vztwobinerr,\
        vzfourbin,vzfourbinerr,\
        backampbin,backampbinerr,\
        backmeanbin,backmeanbinerr,backsigbin,backsigbinerr,\
        vsp1,vsp1err,vsp2,vsp2err,\
        ranal,vzfourstore = \
    velfit(Rkin,vz,vzerr,mskin,Nbinkin,\
           vfitmin,vfitmax,\
           p0vin_min,p0vin_max,p0best,\
           alpmin,alpmax,nsamples,outfile)
print('Fitted VSP1: %f+/-%f' % (vsp1,vsp1err))
print('Fitted VSP2: %f+/-%f' % (vsp2,vsp2err))


###########################################################
#Store the surface density, velocity dispersions, VSPs,
#best fit surfden parameters and Rhalf for GravSphere:
f = open(outfile+'_surfden.txt','w')
for i in range(len(R)):
    f.write('%f %f %f\n' % \
            (R[i], surfden[i], surfdenerr[i]))
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
    f.write('%f %f %f %f %f %f %f %f %f %f %f %f %f\n' % \
            (rbin[i],vzmeanbin[i],vzmeanbinerr[i],vztwobin[i],\
             vztwobinerr[i],\
             vzfourbin[i],vzfourbinerr[i],\
             backampbin[i],backampbinerr[i],\
             backmeanbin[i],backmeanbinerr[i],backsigbin[i],\
             backsigbinerr[i]))
f.close()
f = open(outfile+'_vsps.txt','w')
f.write('%f %f\n' % (vsp1,vsp1err))
f.write('%f %f\n' % (vsp2,vsp2err))
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
             linewidth=mylinewidth,marker='o',markersize=5,alpha=0.75)
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

plt.axvline(x=Rhalf,color='black',alpha=0.5,\
            linewidth=mylinewidth)
plt.axvline(x=Rhalf_int[0],color='blue',alpha=0.5,\
            linewidth=mylinewidth)

plt.xlim([xpltmin,xpltmax])
plt.ylim([surfpltmin,surfpltmax])

plt.savefig(outfile+'_surfden.pdf',bbox_inches='tight')

##### Velocity dispersion fit ##### 
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

plt.errorbar(rbin,vztwobin,vztwobinerr,color='black',ecolor='black',\
             fmt='o',\
             linewidth=mylinewidth,marker='o',markersize=5,alpha=0.75)

plt.axvline(x=Rhalf,color='black',alpha=0.5,\
            linewidth=mylinewidth)
plt.axvline(x=Rhalf_int[0],color='blue',alpha=0.5,\
            linewidth=mylinewidth)

plt.xlim([xpltmin,xpltmax])
plt.ylim([vztwopltmin,vztwopltmax])

plt.savefig(outfile+'_vztwo.pdf',bbox_inches='tight')

##### Fourth velocity moment fit ##### 
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

plt.errorbar(rbin,vzfourbin,vzfourbinerr,color='black',ecolor='black',\
             fmt='o',\
             linewidth=mylinewidth,marker='o',markersize=5,alpha=0.75)

#Plot individual fits + falloff:
for i in range(nsamples):
    plt.plot(ranal,vzfourstore[i,:],linewidth=1,alpha=0.25)

plt.axvline(x=Rhalf,color='black',alpha=0.5,\
            linewidth=mylinewidth)
plt.axvline(x=Rhalf_int[0],color='blue',alpha=0.5,\
            linewidth=mylinewidth)

plt.xlim([xpltmin,xpltmax])
plt.ylim([vzfourpltmin,vzfourpltmax])

plt.savefig(outfile+'_vzfour.pdf',bbox_inches='tight')


###########################################################
#Exit:
print('\nThank you for using the Binulator! Have a nice day.\n')
