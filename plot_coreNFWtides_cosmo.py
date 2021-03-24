###########################################################
#plot_coreNFWtides
###########################################################

#Python programme to plot the coreNFWtides model in CDM
#WDM and SIDM.

###########################################################
#Main code:

#Imports & dependencies:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.optimize import curve_fit
from functions import *
from constants import *
from figures import *
import sys

#Set LCDM density profile (WLM-like):
ranal = np.logspace(-3,1,np.int(250))
M200 = 1.0e10
rc_kpc = 1.0
n = 0.0
rt_kpc = 100.0
delta = 5.0
c200 = cosmo_cfunc(M200,h)
rho_LCDM = corenfw_tides_den(ranal,M200,c200,rc_kpc,n,rt_kpc,delta)

#And now for WDM:
mWDM_keV = 1.0
c200_WDM = cosmo_cfunc_WDM(M200,h,OmegaM,rhocrit,mWDM_keV)
rho_LWDM = corenfw_tides_den(ranal,M200,c200_WDM,rc_kpc,n,rt_kpc,delta)

#And for SIDM calculate the cross section assuming complete
#core formation (n=1):
rho_SIDM = corenfw_tides_den(ranal,M200,c200,rc_kpc,1.0,rt_kpc,delta)
sigm = sidm_novel(rc_kpc,M200,c200,oden,rhocrit)
print('SIDM cross section for this model: %f cm^2/g' % (sigm))

###########################################################
#Plot:

fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(mylinewidth)
ax.minorticks_on()
ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
plt.xticks(fontsize=myfontsize)
plt.yticks(fontsize=myfontsize)

plt.loglog()

plt.xlabel(r'Radius\,[kpc]',fontsize=myfontsize)
plt.ylabel(r'Density\,[${\rm M}_\odot$\,kpc$^{-3}$]',fontsize=myfontsize)

plt.plot(ranal,rho_LCDM,linewidth=mylinewidth,color='blue',\
         label='LCDM',alpha=0.5)
plt.plot(ranal,rho_LWDM,linewidth=mylinewidth,color='red',\
         label='LWDM (1\,keV)',alpha=0.5)
plt.plot(ranal,rho_SIDM,linewidth=mylinewidth,color='green',\
         label='LSIDM (%.1f\,cm$^2$g$^{-1}$)' % (sigm),alpha=0.5)

plt.legend(fontsize=mylegendfontsize)
plt.xlim([1e-2,5])
plt.ylim([1e6,1e10])
plt.legend(loc='upper right',fontsize=18)
plt.savefig(output_base+'coreNFWtides_cosmo.pdf',bbox_inches='tight')
