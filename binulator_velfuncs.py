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

#Select velocity distribution functin:
velpdfuse = velpdffast

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
    #and their confidence intervals; and
    #the mean and dispersion of the background
    #model.
    rbin = np.zeros(len(R))
    right_bin_edge = np.zeros(len(R))
    vzmeanbin = np.zeros(len(R))
    vzmeanbinlo = np.zeros(len(R))    
    vzmeanbinhi = np.zeros(len(R))    
    vztwobin = np.zeros(len(R))
    vztwobinlo = np.zeros(len(R))
    vztwobinhi = np.zeros(len(R))
    vzfourbin = np.zeros(len(R))
    vzfourbinlo = np.zeros(len(R))
    vzfourbinhi = np.zeros(len(R))
    backampbin = np.zeros(len(R))
    backampbinlo = np.zeros(len(R))
    backampbinhi = np.zeros(len(R))
    backmeanbin = np.zeros(len(R))
    backmeanbinlo = np.zeros(len(R))
    backmeanbinhi = np.zeros(len(R))
    backsigbin = np.zeros(len(R))
    backsigbinlo = np.zeros(len(R))
    backsigbinhi = np.zeros(len(R))
    
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
            right_bin_edge[cnt] = R[index[i]]
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
           
            vzmeanbin[cnt],vzmeanbinlo[cnt],vzmeanbinhi[cnt],\
            vztwobin[cnt],vztwobinlo[cnt],vztwobinhi[cnt],\
            vzfourbin[cnt],vzfourbinlo[cnt],vzfourbinhi[cnt],\
            backampbin[cnt],backampbinlo[cnt],backampbinhi[cnt],\
            backmeanbin[cnt],backmeanbinlo[cnt],backmeanbinhi[cnt],\
            backsigbin[cnt],backsigbinlo[cnt],backsigbinhi[cnt],\
            vzfour_store,p0vbest = \
                velfitbin(vzuse,vzerruse,msuse,\
                          p0vin_min,p0vin_max,nsamples)
            vzfour_pdf[:,cnt] = vzfour_store

            #Output the fit values:
            print('Bin %d | vzmean %f+%f-%f | vztwo %f+%f-%f | vzfour %f+%f-%f | bamp %f+%f-%f | bmean %f+%f-%f | bsig %f+%f-%f\n' \
                  % (cnt,vzmeanbin[cnt],vzmeanbin[cnt]-vzmeanbinlo[cnt],vzmeanbinhi[cnt]-vzmeanbin[cnt],\
                     vztwobin[cnt],vztwobinhi[cnt]-vztwobin[cnt],vztwobin[cnt]-vztwobinlo[cnt],\
                     vzfourbin[cnt],vzfourbinhi[cnt]-vzfourbin[cnt],vzfourbin[cnt]-vzfourbinlo[cnt],\
                     backampbin[cnt],backampbinhi[cnt]-backampbin[cnt],backampbin[cnt]-backampbinlo[cnt],\
                     backmeanbin[cnt],backmeanbinhi[cnt]-backmeanbin[cnt],backmeanbin[cnt]-backmeanbinlo[cnt],\
                     backsigbin[cnt],backsigbinhi[cnt]-backsigbin[cnt],backsigbin[cnt]-backsigbinlo[cnt])
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
            
            n, bins, patches = plt.hist(vzuse,10,weights=msuse,\
                                        facecolor='g',\
                                        alpha=0.75)
            vplot = np.linspace(-150,150,np.int(500))
            vperr = np.zeros(len(vplot))+\
                np.sum(vzerruse*msuse)/np.sum(msuse)
            pdf = velpdfuse(vplot,vperr,p0vbest)
            plt.plot(vplot,pdf/np.max(pdf)*np.max(n),\
                     linewidth=mylinewidth)
            plt.xlim([-150,150])
            plt.savefig(outfile+'vzhist_%d.pdf' % (cnt),\
                bbox_inches='tight')

            #Move on to the next bin:
            if (cnt == 0):
                rbin[cnt] = right_bin_edge[cnt]/2.0
            else:
                rbin[cnt] = \
                    (right_bin_edge[cnt] + right_bin_edge[cnt-1])/2.0
            jsum = 0.0
            js = 0
            cnt = cnt + 1

    #Cut back the output arrays:
    rbin = rbin[:cnt]
    vzmeanbin = vzmeanbin[:cnt]
    vzmeanbinlo = vzmeanbinlo[:cnt]
    vzmeanbinhi = vzmeanbinhi[:cnt]
    vztwobin = vztwobin[:cnt]
    vztwobinlo = vztwobinlo[:cnt]
    vztwobinhi = vztwobinhi[:cnt]
    vzfourbin = vzfourbin[:cnt]
    vzfourbinlo = vzfourbinlo[:cnt]
    vzfourbinhi = vzfourbinhi[:cnt]
    backampbin = backampbin[:cnt]
    backampbinlo = backampbinlo[:cnt]
    backampbinhi = backampbinhi[:cnt]
    backmeanbin = backmeanbin[:cnt]
    backmeanbinlo = backmeanbinlo[:cnt]
    backmeanbinhi = backmeanbinhi[:cnt]
    backsigbin = backsigbin[:cnt]
    backsigbinlo = backsigbinlo[:cnt]
    backsigbinhi = backsigbinhi[:cnt]
    
    #Calculate the VSPs with uncertainties. This
    #assumes negligible error in the surface density 
    #profile as compared to the velocity uncertainties.
    #This is usually fine, but something to bear in mind.
    ranal = np.logspace(-3,3,np.int(5e3))
    surfden = threeplumsurf(ranal,p0best[0],p0best[1],p0best[2],\
                            p0best[3],p0best[4],p0best[5])

    #Regularize:
    regularize = 'no'
    if (regularize == 'yes'):
        from scipy.signal import savgol_filter
        #Set filter size and regularize, if
        #enough data to warrant it:
        if (cnt > 4):
            filt_size = cnt / 5.0
            if (filt_size % 2 == 0): filt_size += 1
            if (filt_size < 3):
                filt_size = 3
            if (filt_size > 21):
                filt_size = 21
            
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

        #Regularize:
        vzfour_reg = vzfour_thissample
        if (regularize == 'yes'):
            if (cnt > 4):
                vzfour_reg = \
                    savgol_filter(vzfour_thissample,np.int(filt_size),2)
        vzfour = vzfourfunc(ranal,rbin,vzfour_reg,alp)

        vzfourstore[i,:] = vzfour
        vsp1[i] = integrator(surfden*vzfour*ranal,ranal)
        vsp2[i] = integrator(surfden*vzfour*ranal**3.0,ranal)
    vsp1_int[0], vsp1_int[1], vsp1_int[2], vsp1_int[3], \
        vsp1_int[4], vsp1_int[5], vsp1_int[6] = \
        calcmedquartnine(vsp1)
    vsp2_int[0], vsp2_int[1], vsp2_int[2], vsp2_int[3], \
        vsp2_int[4], vsp2_int[5], vsp2_int[6] = \
        calcmedquartnine(vsp2)
    
    return rbin,vzmeanbin,vzmeanbinlo,vzmeanbinhi,\
        vztwobin,vztwobinlo,vztwobinhi,\
        vzfourbin,vzfourbinlo,vzfourbinhi,\
        backampbin,backampbinlo,backampbinhi,\
        backmeanbin,backmeanbinlo,backmeanbinhi,\
        backsigbin,backsigbinlo,backsigbinhi,\
        vsp1_int[0],vsp1_int[1],vsp1_int[2],\
        vsp2_int[0],vsp2_int[1],vsp2_int[2],\
        ranal,vzfourstore,vsp1,vsp2
    
def velfitbin(vz,vzerr,ms,p0vin_min,p0vin_max,nsamples):
    #Fit the model velocity pdf to a single bin:

    #Functions:
    def lnprob_vel(theta, y, yerr, ms):
        lp = lnprior_vel(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike_vel(theta, y, yerr, ms)

    def lnlike_vel(theta, y, yerr, ms):
        modelpdf = velpdfuse(y,yerr,theta)
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
    nmodels = 10000

    #Starting guess
    ndims = len(p0vin_min)
    pos = np.zeros((nwalkers, ndims), dtype='float')
    p0vin_startmin = p0vin_min
    p0vin_startmax = p0vin_max    
    for i in range(ndims):
        pos[:,i] = np.random.uniform(p0vin_startmin[i],\
            p0vin_startmax[i],nwalkers)

    #Run chains:
    sampler = emcee.EnsembleSampler(nwalkers, ndims, lnprob_vel, \
                args=(vz, vzerr, ms))
    sampler.run_mcmc(pos, nmodels)

    #Extract results + errors:
    burn = np.int(0.75*nmodels)
    chisq = -2.0 * \
            sampler.get_log_prob(discard=burn, flat=True)
    par_test = sampler.get_chain(discard=burn, flat=True)

    #Store best fit model:
    index = np.argsort(chisq)
    p0best = par_test[index[0],:]

    #Choose number of models to draw from the chains:
    min_chisq = np.min(chisq)
    index = np.where(chisq < min_chisq*500.0)[0]
    sample_choose = index[np.random.randint(len(index), \
                          size=nsamples)]

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

    #Output median and 68% confidence intervals.
    #Pass back full vzfour distribution for VSP 
    #calculation.
    return vzmean_int[0],vzmean_int[1],vzmean_int[2],\
           vztwo_int[0],vztwo_int[1],vztwo_int[2],\
           vzfour_int[0],vzfour_int[1],vzfour_int[2],\
           backamp_int[0],backamp_int[1],backamp_int[2],\
           backmean_int[0],backmean_int[1],backmean_int[2],\
           backsig_int[0],backsig_int[1],backsig_int[2],\
           vzfour_store,p0best
