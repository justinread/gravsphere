###########################################################
#plot_velpdf
###########################################################

#Python programme to plot the velpdf function used by
#binulator to fit each velocity bin.

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

#Set up example vel PDF plot. Parameters are:
#[vzmean,alp,bet,backamp,backmean,backsig]
beta_range = np.array([1.0,2.0,5.0])
mycol = ['black','blue','green']
pars_beta_range = np.array([0.0,15.0,0.0,0.0,0.0,50.0])
vz = np.linspace(-50,50,np.int(1e3))
vzerr = np.zeros(len(vz)) + 2.0

###########################################################
#Plots:

#Sequence of velpdf showing its kurtosis range for the
#true PDF convolved with Gaussian velocity errors, and
#a fast analytic approximation. Includes a stacked error
#residual plot.

fig = plt.figure(figsize=(figx,figy))
plt.rcParams['font.size'] = myfontsize

ax = fig.add_subplot(211)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(mylinewidth)
ax.minorticks_on()
ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.set_ylabel(r'frequency',fontsize=myfontsize)

for i in range(len(beta_range)):
    pars_beta_range[2] = beta_range[i]
    pdfnoerr = velpdf_noerr(vz,pars_beta_range)
    pdf = velpdf(vz,vzerr,pars_beta_range)
    pdffast = velpdffast(vz,vzerr,pars_beta_range)
    kurt = kurt_calc(pars_beta_range)
    plt.plot(vz,pdffast,linewidth=mylinewidth,\
             color=mycol[i],linestyle='dashed')
    plt.plot(vz,pdf,linewidth=mylinewidth,color=mycol[i],\
             label=r'$\beta_v,\kappa = %.1f,%.1f$' % \
             (beta_range[i],kurt))

ax.set_xlim([-50,50])
plt.legend(fontsize=16,loc='upper left')

ax2 = fig.add_subplot(212)
for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(mylinewidth)
ax2.minorticks_on()
ax2.tick_params('both', length=10, width=2, which='major')
ax2.tick_params('both', length=5, width=1, which='minor')

ax2.set_xlabel(r'v$_{\rm los}$ (km/s)',fontsize=myfontsize)
ax2.set_ylabel(r'Residual (\%)',fontsize=myfontsize)

for i in range(len(beta_range)):
    pars_beta_range[2] = beta_range[i]
    pdf = velpdf(vz,vzerr,pars_beta_range)
    pdffast = velpdffast(vz,vzerr,pars_beta_range)
    residual = (pdf-pdffast) / np.max(pdf) * 100.0
    
    plt.plot(vz,residual,linewidth=mylinewidth,\
             color=mycol[i])

ax2.set_xlim([-50,50])
#ax2.set_ylim([-20,20])

plt.savefig(output_base+'velpdf_kurtosis_example.pdf',\
            bbox_inches='tight')
