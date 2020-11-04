import numpy as np
from scipy.integrate.quadrature import simps as integrator
from scipy.special import gamma
from constants import *
from functions import * 
import emcee

#Functions for emcee (surface density):
def tracerfit(R,surfden,surfdenerr,Rfitmin,Rfitmax,p0in_min,p0in_max):
    #Code to fit the surface density profile with
    #a three-Plummer model:

    #Functions:
    def lnprob_surf(theta, x, y, yerr):
        lp = lnprior_surf(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike_surf(theta, x, y, yerr)

    def lnlike_surf(theta, x, y, yerr):
        model = threeplumsurf(x,theta[0],theta[1],theta[2],\
            theta[3],theta[4],theta[5])

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
        if all(minarr < thetau < maxarr for minarr,thetau,maxarr in \
               zip(minarr,theta,maxarr)):
            return 0.0
        return -np.inf
        
    lnprior_surf = lambda theta: \
        lnprior_set_surf(theta,p0in_min,p0in_max)

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
    
    #Emcee parameters:
    nwalkers = 250
    nmodels = 1000

    #Starting guess
    ndims = len(p0in_min)
    pos = np.zeros((nwalkers, ndims), dtype='float')
    p0in_startmin = p0in_min
    p0in_startmax = p0in_max 
    for i in range(ndims):
        pos[:,i] = np.random.uniform(p0in_startmin[i],\
            p0in_startmax[i],nwalkers)

    #Run chains:
    sampler = emcee.EnsembleSampler(nwalkers, ndims, lnprob_surf, \
                args=(Rfit, surfdenfit, surfdenerrfit))
    sampler.run_mcmc(pos, nmodels)

    #Extract results + errors:
    burn = np.int(0.5*nmodels)
    chisq = -2.0 * sampler.lnprobability[:, burn:].reshape(-1)
    par_test = np.zeros((len(chisq),ndims), dtype='float')
    for i in range(ndims):
        par_test[:,i] = sampler.chain[:, burn:, i].reshape(-1)

    #Store best fit model:
    index = np.argsort(chisq)
    p0best = par_test[index[0],:]

    #Choose number of models to draw from the chains:
    nsamples = 1000
    sample_choose = np.random.randint(len(chisq), \
                        size=nsamples)

    #Set up arrays to store med,68%,95%,99% confidence intervals:
    Rplot = np.logspace(-3,3,1000)
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
        surf_int[0,j], surf_int[1,j], surf_int[2,j], surf_int[3,j], \
            surf_int[4,j], surf_int[5,j], surf_int[6,j] = \
            calcmedquartnine(surf_store[j,:])
    Rhalf_int[0], Rhalf_int[1], Rhalf_int[2], Rhalf_int[3], \
        Rhalf_int[4], Rhalf_int[5], Rhalf_int[6] = \
        calcmedquartnine(Rhalf_store)

    return Rplot,surf_int,Rhalf_int,p0best
