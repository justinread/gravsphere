import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.integrate.quadrature import simps as integrator
from scipy.special import gamma
from constants import *
from functions import * 
from figures import *
import emcee
import sys

#Functions for emcee (velocity data):
def velfit(R,vz,vzerr,ms,Nbin,\
           vfitmin,vfitmax,\
           p0vin_min,p0vin_max,p0best,\
           alpmin,alpmax,nsamples,outfile):
    #Code to fit the velocity data with
    #the velpdf model.
    
    #Arrays to store binned:
    #radius (rbin),
    #<vlos> (vzmeanbin), 
    #<vlos^2>^1/2 (vztwobin),
    #<vlos^4> (vzfourbin),
    #and their errors; and
    #the mean and dispersion of the background
    #model.
    rbin = np.zeros(len(R))
    vzmeanbin = np.zeros(len(R))
    vzmeanbinerr = np.zeros(len(R))    
    vztwobin = np.zeros(len(R))
    vztwobinerr = np.zeros(len(R))
    vzfourbin = np.zeros(len(R))
    vzfourbinerr = np.zeros(len(R))
    backampbin = np.zeros(len(R))
    backampbinerr = np.zeros(len(R))
    backmeanbin = np.zeros(len(R))
    backmeanbinerr = np.zeros(len(R))
    backsigbin = np.zeros(len(R))
    backsigbinerr = np.zeros(len(R))
    
    #This for storing the vzfour pdf
    #for calculating the VSPs:
    vzfour_pdf = np.zeros((nsamples,len(R)))

    #Loop through the bins, assuming Nbin stars
    #(weighted by ms) per bin:
    index = np.argsort(R)
    vzstore = np.zeros(len(R))
    vzerrstore = np.zeros(len(R))
    msstore = np.zeros(len(R))
    cnt = 0
    jsum = 0
    js = 0
    for i in range(len(R)):
        #Find stars in bin:
        if (jsum < Nbin):
            vzstore[js] = vz[index[i]]
            vzerrstore[js] = vzerr[index[i]]
            msstore[js] = ms[index[i]]
            rbin[cnt] = R[index[i]]
            jsum = jsum + ms[index[i]]
            js = js + 1
        if (jsum >= Nbin):
            #Fit the velpdf model to these stars:
            vzuse = vzstore[:js]
            vzerruse = vzerrstore[:js]
            msuse = msstore[:js]
            
            #Cut back to fit range. If doing this,
            #perform an initial centre for this
            #bin first to ensure a symmetric
            #"haircut" on the data.
            if (vfitmin != 0 or vfitmax != 0):
                vzuse = vzuse - \
                    np.sum(vzuse*msuse)/np.sum(msuse)
            if (vfitmin != 0):
                vzuse_t = vzuse[vzuse > vfitmin]
                vzerruse_t = vzerruse[vzuse > vfitmin]
                msuse_t = msuse[vzuse > vfitmin]
            else:
                vzuse_t = vzuse
                vzerruse_t = vzerruse
                msuse_t = msuse
            if (vfitmax != 0):
                vzuse = vzuse_t[vzuse_t < vfitmax]
                vzerruse = vzerruse_t[vzuse_t < vfitmax]
                msuse = msuse_t[vzuse_t < vfitmax]
            else:
                vzuse = vzuse_t
                vzerruse = vzerruse_t
                msuse = msuse_t
           
            vzmeanbin[cnt],vzmeanbinerr[cnt],vztwobin[cnt],\
            vztwobinerr[cnt],\
            vzfourbin[cnt],vzfourbinerr[cnt],\
            backampbin[cnt],backampbinerr[cnt],\
            backmeanbin[cnt],backmeanbinerr[cnt],\
            backsigbin[cnt],backsigbinerr[cnt],\
            vzfour_store,p0vbest = \
                velfitbin(vzuse,vzerruse,msuse,\
                          p0vin_min,p0vin_max,nsamples)
            vzfour_pdf[:,cnt] = vzfour_store

            #Output the fit values:
            print('Bin %d | vzmean %f+/-%f | vztwo %f+/-%f | vzfour %f+/-%f | bamp %f+/-%f | bmean %f+/-%f | bsig %f+/-%f\n' \
                  % (cnt,vzmeanbin[cnt],vzmeanbinerr[cnt],\
                     vztwobin[cnt],vztwobinerr[cnt],\
                     vzfourbin[cnt],vzfourbinerr[cnt],\
                     backampbin[cnt],backampbinerr[cnt],\
                     backmeanbin[cnt],backmeanbinerr[cnt],\
                     backsigbin[cnt],backsigbinerr[cnt])
                 )
            
            #Make a plot of the fit:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(mylinewidth)

            ax.minorticks_on()
            ax.tick_params('both', length=10, width=2, which='major')
            ax.tick_params('both', length=5, width=1, which='minor')
            plt.xticks(fontsize=myfontsize)
            plt.yticks(fontsize=myfontsize)

            plt.xlabel(r'v$_z$ [km/s]',fontsize=myfontsize)
            plt.ylabel(r'frequency',\
                fontsize=myfontsize)
            
            n, bins, patches = plt.hist(vzuse,50,weights=msuse,\
                                        facecolor='g',\
                                        alpha=0.75)
            vplot = np.linspace(-50,50,np.int(500))
            vperr = np.zeros(len(vplot))+\
                np.sum(vzerruse*msuse)/np.sum(msuse)
            pdf = velpdffast(vplot,vperr,p0vbest)
            plt.plot(vplot,pdf/np.max(pdf)*np.max(n),\
                     linewidth=mylinewidth)
            plt.savefig(outfile+'vzhist_%d.pdf' % (cnt),\
                bbox_inches='tight')

            #Move on to the next bin:
            jsum = 0.0
            js = 0
            cnt = cnt + 1

    #Cut back the output arrays:
    rbin = rbin[:cnt]
    vzmeanbin = vzmeanbin[:cnt]
    vzmeanbinerr = vzmeanbinerr[:cnt]
    vztwobin = vztwobin[:cnt]
    vztwobinerr = vztwobinerr[:cnt]
    vzfourbin = vzfourbin[:cnt]
    vzfourbinerr = vzfourbinerr[:cnt]
    backampbin = backampbin[:cnt]
    backampbinerr = backampbinerr[:cnt]
    backmeanbin = backmeanbin[:cnt]
    backmeanbinerr = backmeanbinerr[:cnt]
    backsigbin = backsigbin[:cnt]
    backsigbinerr = backsigbinerr[:cnt]
    
    #Calculate the VSPs with uncertainties. This
    #assumes negligible error in the surface density 
    #profile as compared to the velocity uncertainties.
    #This is usually fine, but something to bear in mind.
    ranal = np.logspace(-3,3,np.int(5e3))
    surfden = threeplumsurf(ranal,p0best[0],p0best[1],p0best[2],\
                            p0best[3],p0best[4],p0best[5])

    #This assumes a flat or linearly falling relation
    #beyond the last data point:
    vsp1 = np.zeros(nsamples)
    vsp2 = np.zeros(nsamples)
    vsp1_int = np.zeros(7)
    vsp2_int = np.zeros(7)
    vzfourstore = np.zeros((nsamples,len(ranal)))
    for i in range(nsamples):
        vzfour_thissample = vzfour_pdf[i,:cnt]
        alp = np.random.random()*(alpmax-alpmin)+alpmin
        vzfour = vzfourfunc(ranal,rbin,vzfour_thissample,alp)
        vzfourstore[i,:] = vzfour
        vsp1[i] = integrator(surfden*vzfour*ranal,ranal)
        vsp2[i] = integrator(surfden*vzfour*ranal**3.0,ranal)
    vsp1_int[0], vsp1_int[1], vsp1_int[2], vsp1_int[3], \
        vsp1_int[4], vsp1_int[5], vsp1_int[6] = \
        calcmedquartnine(vsp1)
    vsp2_int[0], vsp2_int[1], vsp2_int[2], vsp2_int[3], \
        vsp2_int[4], vsp2_int[5], vsp2_int[6] = \
        calcmedquartnine(vsp2)
    vsp1out = vsp1_int[0]
    vsp1outerr = (vsp1_int[2]-vsp1_int[1])/2.0
    vsp2out = vsp2_int[0]
    vsp2outerr = (vsp2_int[2]-vsp2_int[1])/2.0
    
    return rbin,vzmeanbin,vzmeanbinerr,vztwobin,vztwobinerr,\
        vzfourbin,vzfourbinerr,\
        backampbin,backampbinerr,\
        backmeanbin,backmeanbinerr,backsigbin,backsigbinerr,\
        vsp1out,vsp1outerr,vsp2out,vsp2outerr,\
        ranal,vzfourstore
    
def velfitbin(vz,vzerr,ms,p0vin_min,p0vin_max,nsamples):
    #Fit the model velocity pdf to a single bin:

    #Functions:
    def lnprob_vel(theta, y, yerr, ms):
        lp = lnprior_vel(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike_vel(theta, y, yerr, ms)

    def lnlike_vel(theta, y, yerr, ms):
        modelpdf = velpdffast(y,yerr,theta)
        lnlike_out = np.sum(np.log(modelpdf)*ms)
    
        if (lnlike_out != lnlike_out):
            lnlike_out = -np.inf
    
        return lnlike_out
    
    def lnprior_set_vel(theta,p0vin_min,p0vin_max):
        ndims = len(theta)
        minarr = np.zeros(ndims)
        maxarr = np.zeros(ndims)
        for i in range(ndims):
            minarr[i] = p0vin_min[i]
            maxarr[i] = p0vin_max[i]
       
        if all(minarr < thetau < maxarr for minarr,thetau,maxarr in \
               zip(minarr,theta,maxarr)):
            return 0.0
        return -np.inf

    lnprior_vel = lambda theta: \
                lnprior_set_vel(theta,p0vin_min,p0vin_max)

    #Emcee parameters:
    nwalkers = 250
    nmodels = 5000

    #Starting guess
    ndims = len(p0vin_min)
    pos = np.zeros((nwalkers, ndims), dtype='float')
    for i in range(ndims):
        pos[:,i] = np.random.uniform(p0vin_min[i],p0vin_max[i],nwalkers)

    #Run chains:
    sampler = emcee.EnsembleSampler(nwalkers, ndims, lnprob_vel, \
                args=(vz, vzerr, ms))
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
    sample_choose = np.random.randint(len(chisq), \
                        size=nsamples)    

    #Set up arrays to store med,68%,95%,99% confidence intervals:
    vzmean_int = np.zeros(7)
    vzmean_store = np.zeros(nsamples)
    vztwo_int = np.zeros(7)
    vztwo_store = np.zeros(nsamples)
    vzfour_int = np.zeros(7)
    vzfour_store = np.zeros(nsamples)
    backamp_int = np.zeros(7)
    backamp_store = np.zeros(nsamples) 
    backmean_int = np.zeros(7)
    backmean_store = np.zeros(nsamples)
    backsig_int = np.zeros(7)
    backsig_store = np.zeros(nsamples)
 
    for i in range(nsamples):
        theta = par_test[sample_choose[i],:]
        vzmean_store[i] = theta[0]
        vztwo_store[i] = vztwo_calc(theta)
        vzfour_store[i] = vzfour_calc(theta)
        backamp_store[i] = theta[3]
        backmean_store[i] = theta[4]
        backsig_store[i] = theta[5]

    #Solve for confidence intervals:
    vzmean_int[0], vzmean_int[1], vzmean_int[2], vzmean_int[3], \
        vzmean_int[4], vzmean_int[5], vzmean_int[6] = \
        calcmedquartnine(vzmean_store)
    vztwo_int[0], vztwo_int[1], vztwo_int[2], vztwo_int[3], \
        vztwo_int[4], vztwo_int[5], vztwo_int[6] = \
        calcmedquartnine(vztwo_store)
    vzfour_int[0], vzfour_int[1], vzfour_int[2], vzfour_int[3], \
        vzfour_int[4], vzfour_int[5], vzfour_int[6] = \
        calcmedquartnine(vzfour_store)
    backamp_int[0], backamp_int[1], backamp_int[2], backamp_int[3], \
        backamp_int[4], backamp_int[5], backamp_int[6] = \
        calcmedquartnine(backamp_store)
    backmean_int[0], backmean_int[1], backmean_int[2], backmean_int[3], \
        backmean_int[4], backmean_int[5], backmean_int[6] = \
        calcmedquartnine(backmean_store)
    backsig_int[0], backsig_int[1], backsig_int[2], backsig_int[3], \
        backsig_int[4], backsig_int[5], backsig_int[6] = \
        calcmedquartnine(backsig_store)

    #Output median and symmetrised 68% confidence errors.
    #Pass back full vzfour distribution for VSP 
    #calculation.
    return vzmean_int[0],(vzmean_int[2]-vzmean_int[1])/2.0,\
           vztwo_int[0],(vztwo_int[2]-vztwo_int[1])/2.0,\
           vzfour_int[0],(vzfour_int[2]-vzfour_int[1])/2.0,\
           backamp_int[0],(backamp_int[2]-backamp_int[1])/2.0,\
           backmean_int[0],(backmean_int[2]-backmean_int[1])/2.0,\
           backsig_int[0],(backsig_int[2]-backsig_int[1])/2.0,\
           vzfour_store,p0best
