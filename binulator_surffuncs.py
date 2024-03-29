import numpy as np
from scipy.integrate import simps as integrator
from scipy.special import gamma
from constants import *
from functions import * 
import emcee
from multiprocessing import Pool
from multiprocessing import cpu_count

#Functions:
def lnprob_surf(theta, x, y, yerr, p0in_min, p0in_max):
    lp = lnprior_set_surf(theta,p0in_min,p0in_max)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike_surf(theta, x, y, yerr)

def lnlike_surf(theta, x, y, yerr):
    if (theta[0] < 0 or theta[1] < 0 or theta[2] < 0):
        #If using neg. Plummer component, add some more
        #"data points" at large & small radii bounded on
        #zero and the outermost data point. This
        #will disfavour models with globally
        #negative tracer density.
        x,y,yerr = Sig_addpnts(x,y,yerr)
        
    model = threeplumsurf(x,theta[0],theta[1],theta[2],\
               theta[3],theta[4],theta[5])
        
    #If using negative Plummer components,
    #shrink the error to disfavour globally
    #negative models. The shrinkage amount
    #is designed to ensure the likelihood is
    #penalised so we pick the smallest error
    #on any point, divide it by the total
    #number of points and then divide that
    #by 1e3 for good measure.
    if (theta[0] < 0 or theta[1] < 0 or theta[2] < 0):
        if (np.min(model) < 0):
            yerr[np.where(model < 0)] = \
                 np.min(yerr)/np.float(len(x))/1.0e3

    inv_sigma2 = 1.0/(yerr)**2.0
    lnlike_out = -0.5*(np.sum((y-model)**2*inv_sigma2))
    
    if (lnlike_out != lnlike_out):
        lnlike_out = -np.inf
            
    return lnlike_out

def lnprior_set_surf(theta,p0in_min,p0in_max):
    ndims = len(theta)
    minarr = np.zeros(ndims)
    maxarr = np.zeros(ndims)
    for i in range(ndims):
        minarr[i] = p0in_min[i]
        maxarr[i] = p0in_max[i]
    if all(minarr < thetau < maxarr for \
           minarr,thetau,maxarr in \
           zip(minarr,theta,maxarr)):
        return 0.0
    return -np.inf

#Functions for emcee (surface density):
def tracerfit(R,surfden,surfdenerr,p0in_min,p0in_max,nprocs):
    #Code to fit the surface density profile with
    #a three-Plummer model:

    #Emcee parameters:
    nwalkers = 500
    nmodels = 10000

    #Starting guess:
    ndims = len(p0in_min)
    pos = np.zeros((nwalkers, ndims), dtype='float')
    p0in_startmin = p0in_min
    p0in_startmax = p0in_max
    for i in range(ndims):
        pos[:,i] = np.random.uniform(p0in_startmin[i],\
            p0in_startmax[i],nwalkers)

    #Run chains:
    with Pool(processes = nprocs) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndims, \
                    lnprob_surf, \
                    args=(R, surfden, surfdenerr, \
		          p0in_min, p0in_max),pool=pool)
        sampler.run_mcmc(pos, nmodels, progress=True)

    #Extract results + errors:
    burn = np.int(0.75*nmodels)
    chisq = -2.0 * \
            sampler.get_log_prob(discard=burn, flat=True)
    par_test = sampler.get_chain(discard=burn, flat=True)
    
    #Store best fit model:
    index = np.argsort(chisq)
    p0best = par_test[index[0],:]

    #Choose number of models to draw from the chains:
    nsamples = 1000
    sample_choose = np.random.randint(len(chisq), \
                        size=nsamples)

    #Set up arrays to store med,68%,95%,99% confidence intervals:
    Rpmin = np.min(R)/50.0
    Rpmax = np.max(R)*50.0
    Rplot = np.logspace(np.log10(Rpmin),np.log10(Rpmax),1000)
    surf_int = np.zeros((7,len(Rplot)))
    surf_store = np.zeros((len(Rplot),nsamples))
    Rhalf_int = np.zeros(7)
    Rhalf_store = np.zeros(nsamples)

    for i in range(nsamples):
        theta = par_test[sample_choose[i],:]
        surf = threeplumsurf(Rplot,theta[0],theta[1],theta[2],\
                             theta[3],theta[4],theta[5])
        Rhalf = Rhalf_func(theta[0],theta[1],theta[2],\
                           theta[3],theta[4],theta[5])
        surf_store[:,i] = surf
        Rhalf_store[i] = Rhalf

    #Solve for confidence intervals:
    for j in range(len(Rplot)):
        surf_int[0,j], surf_int[1,j], surf_int[2,j], \
            surf_int[3,j], \
            surf_int[4,j], surf_int[5,j], surf_int[6,j] = \
            calcmedquartnine(surf_store[j,:])
    Rhalf_int[0], Rhalf_int[1], Rhalf_int[2], Rhalf_int[3], \
        Rhalf_int[4], Rhalf_int[5], Rhalf_int[6] = \
        calcmedquartnine(Rhalf_store)

    return Rplot,surf_int,Rhalf_int,p0best
