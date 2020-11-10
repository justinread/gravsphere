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

#Set up examples [vzmean,alp,bet,backamp,backmean,backsig]:
pars = np.zeros(6)
pars[0] = -15.0
pars[1] = 10.0
pars[2] = 1.5
pars[3] = 0.5
pars[4] = 35.0
pars[5] = 50.0
beta_range = np.linspace(0.75,9.0,5)
pars_beta_range = np.array([0.0,10.0,0.0,0.0,0.0,50.0])
vz = np.linspace(-75,75,np.int(1e3))
vzerr = np.zeros(len(vz)) + 25.0
pdf = velpdf(vz,vzerr,pars)
pdfmonte = velpdfmonte(vz,vzerr,pars)
pdffast = velpdffast(vz,vzerr,pars)
pdfnoerr = velpdf_noerr(vz,pars)

###########################################################
#Plots:

##### Velpdf #####
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(mylinewidth)
ax.minorticks_on()
ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
plt.xticks(fontsize=myfontsize)
plt.yticks(fontsize=myfontsize)

plt.xlabel(r'v$_z$ [km/s]',fontsize=myfontsize)
plt.ylabel(r'frequency',fontsize=myfontsize)

plt.plot(vz,pdf,linewidth=mylinewidth,color='blue',label='PDF',alpha=0.5)
plt.plot(vz,pdfmonte,linewidth=mylinewidth,color='black',label='PDFmonte',alpha=0.5)
plt.plot(vz,pdffast,linewidth=mylinewidth,color='green',\
    label='PDFfast',alpha=0.5)
plt.plot(vz,pdfnoerr,linewidth=mylinewidth,color='red',alpha=0.5,\
    label='PDF no error')

plt.legend(fontsize=mylegendfontsize)
plt.savefig(output_base+'velpdf_example.pdf',bbox_inches='tight')

#And now plot a sequence of velpdf showing its kurtosis range:
fig = plt.figure(figsize=(figx,figy))
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(mylinewidth)
ax.minorticks_on()
ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
plt.xticks(fontsize=myfontsize)
plt.yticks(fontsize=myfontsize)

plt.xlabel(r'v$_z$ [km/s]',fontsize=myfontsize)
plt.ylabel(r'frequency',fontsize=myfontsize)

for i in range(len(beta_range)):
    pars_beta_range[2] = beta_range[i]
    pdfnoerr = velpdf_noerr(vz,pars_beta_range)
    kurt = kurt_calc(pars_beta_range)
    plt.plot(vz,pdfnoerr,linewidth=mylinewidth,alpha=0.5,\
             label=r'$\beta_v = %.1f \mid \kappa = %.1f$' % \
             (beta_range[i],kurt))

plt.legend(fontsize=15,loc='upper left')
plt.savefig(output_base+'velpdf_kurtosis_example.pdf',bbox_inches='tight')
