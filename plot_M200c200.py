###########################################################
#plot_M200c200
###########################################################

#Python programme to plot the M200c200 relation in CDM
#versus WDM. This is to test cosmology priors on the 
#gravsphere mass models.

###########################################################
#Main code:

#Imports & dependencies:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from functions import *
from constants import *
from figures import *
import sys

#Set up mass array:
M200 = np.logspace(5,12,np.int(500))

#Set up range of mWDM(keV) to explore:
mWDM = np.array([0.1,0.5,1,5,20,50])

###########################################################
#Plot:

##### M200c200 #####
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

plt.xlabel(r'M$_{200}$ [M$_\odot$]',fontsize=myfontsize)
plt.ylabel(r'c$_{200}$',fontsize=myfontsize)

for i in range(len(mWDM)):
    c200 = cosmo_cfunc_WDM(M200,h,OmegaM,rhocrit,mWDM[i])
    plt.plot(M200,c200,linewidth=mylinewidth,\
            label=r'm$_X$=%.2f' % (mWDM[i]))

plt.legend(fontsize=mylegendfontsize)
plt.savefig(output_base+'plot_M200c200.pdf',bbox_inches='tight')

