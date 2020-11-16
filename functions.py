import numpy as np
from scipy.integrate.quadrature import simps as integrator
from scipy.misc.common import derivative
from scipy.special import gamma
from scipy.integrate import quad, dblquad
from constants import *
multimode = 'normal'

###########################################################
#For setting cosmology priors on coreNFWtides parameters.
def cosmo_cfunc(M200,h):
    #From Dutton & Maccio 2014. Requires as input masses 
    #defined in 200c system in units of Msun:
    c = 10.**(0.905 - 0.101 * (np.log10(M200*h)-12.))
    return c
    
def cosmo_cfunc_WDM(M200,h,OmegaM,rhocrit,mWDM):
    #Use formula in https://arxiv.org/pdf/1112.0330.pdf
    #to modify CDM M200-c200 relation to the WDM 
    #one. Assumes mWDM in keV, dimensionless h
    #M200 in Msun and rhocrit in Msun kpc^-3.
    cCDM = cosmo_cfunc(M200,h)
    gamma1 = 15.0
    gamma2 = 0.3
    lamfseff = 0.049*(mWDM)**(-1.11)*\
        (OmegaM/0.25)**(0.11)*(h/0.7)**(1.22)*1000.0
    lamhm = 13.93*lamfseff
    Mhm = 4.0/3.0*np.pi*rhocrit*(lamhm/2.0)**3.0
    cWDM = cCDM * (1.0 + gamma1*Mhm / M200)**(-gamma2)
    return cWDM

###########################################################
#For constraining particle DM models:
def rhoNFW(r,rhos,rs):
    return rhos/((r/rs)*(1.+(r/rs))**2.)

def sidm_novel(rc,M200,c,oden,rhocrit,GammaX):
    #Calculate SIDM model parameters from the coreNFWtides
    #model fit. For this to be valid, the coreNFWtides fit
    #should assume a pure-core model, with n=1. See
    #Read et al. 2018 for further details.
    #Returns cross section/particle mass in cm^2 / g.
    GammaX = 0.005/(1e9*year)
    Guse = G*Msun/kpc
    rho_unit = Msun/kpc**3.0
    rc = np.abs(rc)*10.0
    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=oden*c**3.*gcon/3.
    rv=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    rs=rv/c
    rhos=rhocrit*deltachar

    rhorc = rhoNFW(rc,rhos,rs)
    r = np.logspace(np.log10(rc),np.log10(rs*5000.0),50000)
    rho = rhoNFW(r,rhos,rs)
    mass = M200*gcon*(np.log(1.0 + r/rs)-r/rs/(1.0+r/rs))
    sigvtworc = Guse/rhorc*integrator(mass*rho/r**2.0,r)
    sigvrc = np.sqrt(sigvtworc)

    sigm = np.sqrt(np.pi)*GammaX/(4.0*rhorc*rho_unit*sigvrc)
    return sigm*100.0**2.0/1000.0

def radius_dsph(s, b, distance):
    return np.sqrt((distance * np.sin(b))**2. + s*s)

def integrand(s, b, distance, Mpars):
    value = np.sin(b) * rho(np.array([radius_dsph(s, b, distance)]), Mpars)**2
    return value
                
def get_J(Mpars, distance, r_max):
    #Min/max integration angles in radians:
    b_min = 0.0
    b_max = np.arcsin(r_max/distance)
    
    #This is an appropriate choice for Dwarf galaxies but
    #should be reconsidered for large mass systems:
    Rmaximum = 250.0
    
    #Upper/lower limits:
    s_min_bound = lambda b :  -(Rmaximum**2 - (distance*np.sin(b))**2 )**0.5
    s_max_bound = lambda b : (Rmaximum**2 - (distance*np.sin(b))**2 )**0.5
    
    #Computation J_max:
    Rmax_arr = Rmaximum
    Acc_arr = 1.0e-8
    J_max = dblquad(integrand,b_min,b_max,s_min_bound,\
                    s_max_bound,args=(distance,Mpars),\
                    epsabs=Acc_arr,epsrel=Acc_arr)
    J_max = J_max[0]*kpccm*2.*np.pi*Msunkpc3toGeVcm3**2.0

    #Error checking:
    if (J_max == np.inf):
        print('Argh! Infinite J_max!! Bye bye...')
        sys.exit(0)
        
    if (J_max < 0):
        print('Argh! Negative J_max!! Bye bye...')
        sys.exit(0)

    return J_max  # in GeV2.cm-5


###########################################################
#For DM mass profile:
def corenfw_tides_den(r,M200,c,rc,n,rt,delta):
    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=oden*c**3.*gcon/3.
    rv=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    rs=rv/c

    rhos=rhocrit*deltachar
    rhoanal = rhos/((r/rs)*(1.+(r/rs))**2.)
    manal = M200 * gcon * (np.log(1.0 + r/rs)-r/rs/(1.0+r/rs))

    x = r/np.abs(rc)
    f = np.tanh(x)
    my_manal = manal*f**n
    my_rhoanal = rhoanal*f**n + \
        1.0/(4.*np.pi*r**2.*np.abs(rc))*manal*(1.0-f**2.)*n*f**(n-1.0)
    frt = np.tanh(rt/np.abs(rc))
    manal_rt = M200 * gcon * (np.log(1.0 + rt/rs)-rt/rs/(1.0+rt/rs))
    my_rhoanal_rt = rhos/((rt/rs)*(1.+(rt/rs))**2.)*frt**n + \
        1.0/(4.*np.pi*rt**2.*np.abs(rc))*manal_rt*(1.0-frt**2.)*n*frt**(n-1.0)

    my_rhoanal[r > rt] = my_rhoanal_rt * (r[r > rt]/rt)**(-delta)

    return my_rhoanal

def corenfw_tides_mass(r,M200,c,rc,n,rt,delta):
    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=oden*c**3.*gcon/3.
    rv=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    rs=rv/c

    rhos=rhocrit*deltachar
    rhoanal = rhos/((r/rs)*(1.+(r/rs))**2.)
    manal = M200 * gcon * (np.log(1.0 + r/rs)-r/rs/(1.0+r/rs))

    x = r/np.abs(rc)
    f = np.tanh(x)
    my_manal = manal*f**n

    frt = np.tanh(rt/np.abs(rc))
    manal_rt = M200 * gcon * (np.log(1.0 + rt/rs)-rt/rs/(1.0+rt/rs))
    my_rhoanal_rt = rhos/((rt/rs)*(1.+(rt/rs))**2.)*frt**n + \
        1.0/(4.*np.pi*rt**2.*np.abs(rc))*manal_rt*(1.0-frt**2.)*n*frt**(n-1.0)
    Mrt = manal_rt*frt**n

    my_manal[r > rt] = Mrt + \
        4.0*np.pi*my_rhoanal_rt*rt**3.0/(3.0-delta)*\
        ((r[r > rt]/rt)**(3.0-delta)-1.0)

    return my_manal

def corenfw_tides_dlnrhodlnr(r,M200,c,rc,n,rt,delta):
    dden = derivative(\
        lambda x: corenfw_tides_den(x,M200,c,rc,n,rt,delta),\
        r,dx=1e-6)
    dlnrhodlnr = dden / corenfw_tides_den(r,M200,c,rc,n,rt,delta) * r
    return dlnrhodlnr
    
def vmax_func(M200,c200,h):
    oden = 200.0
    Guse = G*Msun/kpc
    r200=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    
    #This from Sigad et al. 2000 (via Schneider et al. 2017):
    vmax = 0.465*np.sqrt(Guse*M200/r200)/\
           np.sqrt(1.0/c200*np.log(1.0+c200)-(1.0+c200)**(-1.0))
    return vmax/kms
    

###########################################################
#For DM mass profile and VSPs of GC mocks (for overlaying 
#true solution):
def alpbetgamden(r,rho0,r0,alp,bet,gam):
    return rho0*(r/r0)**(-gam)*(1.0+(r/r0)**alp)**((gam-bet)/alp)

def alpbetgamdlnrhodlnr(r,rho0,r0,alp,bet,gam):
    return -gam + (gam-bet)*(r/r0)**alp*(1.0+(r/r0)**alp)**(-1.0)

def alpbetgammass(r,rho0,r0,alp,bet,gam):
    den = rho0*(r/r0)**(-gam)*(1.0+(r/r0)**alp)**((gam-bet)/alp)
    mass = np.zeros(len(r))
    for i in range(3,len(r)):
        mass[i] = integrator(4.0*np.pi*r[:i]**2.*den[:i],r[:i])
    return mass
    
def alpbetgamsigr(r,rho0s,r0s,alps,bets,gams,rho0,r0,alp,bet,gam,ra):
    nu = alpbetgamden(r,rho0s,r0s,alps,bets,gams)
    mass = alpbetgammass(r,rho0,r0,alp,bet,gam)
    gf = gfunc_osipkov(r,ra)
    sigr = np.zeros(len(r))
    for i in range(len(r)-3):
        sigr[i] = 1.0/nu[i]/gf[i] * \
            integrator(Guse*mass[i:]*nu[i:]/r[i:]**2.0*\
                       gf[i:],r[i:])
    return sigr

def osipkov(r,r0):
    return r**2.0/(r**2.0+r0**2.0)
    
def gfunc_osipkov(r,r0):
    n0 = 2.0
    bet0 = 0.0
    betinf = 1.0
    gfunc = r**(2.0*betinf)*\
        ((r0/r)**n0+1.0)**(2.0/n0*(betinf-bet0))
    return gfunc

def alpbetgamvsp(rho0s,r0s,alps,bets,gams,rho0,r0,alp,bet,gam,ra):
    intpnts = 1e4
    r = np.logspace(np.log10(r0s/50.0),np.log10(500.0*r0s),\
                    np.int(intpnts))
    nu = alpbetgamden(r,rho0s,r0s,alps,bets,gams)
    massnu = alpbetgamden(r,rho0s,r0s,alps,bets,gams)
    mass = alpbetgammass(r,rho0,r0,alp,bet,gam)
    sigr = alpbetgamsigr(r,rho0s,r0s,alps,bets,gams,rho0,\
                         r0,alp,bet,gam,ra)
    bet = osipkov(r,ra)
    sigstar = np.zeros(len(r))
    for i in range(1,len(r)-3):
        sigstar[i] = 2.0*integrator(nu[i:]*r[i:]/\
                               np.sqrt(r[i:]**2.0-r[i-1]**2.0),\
                               r[i:])
 
    #Normalise similarly to the data:
    norm = integrator(sigstar*2.0*np.pi*r,r)
    nu = nu / norm
    sigstar = sigstar / norm

    #VSPs:
    vsp1 = \
        integrator(2.0/5.0*Guse*mass*nu*(5.0-2.0*bet)*\
            sigr*r,r)/1.0e12
    vsp2 = \
        integrator(4.0/35.0*Guse*mass*nu*(7.0-6.0*bet)*\
            sigr*r**3.0,r)/1.0e12
        
    #Richardson & Fairbairn zeta parameters:
    Ntotuse = integrator(sigstar*r,r)
    sigint = integrator(sigstar*r**3.0,r)
    zeta_A = 9.0/10.0*Ntotuse*integrator(Guse*mass*nu*(\
        5.0-2.0*bet)*sigr*r,r)/\
        (integrator(Guse*mass*nu*r,r))**2.0
    zeta_B = 9.0/35.0*Ntotuse**2.0*\
        integrator(Guse*mass*nu*(7.0-6.0*bet)*sigr*r**3.0,r)/\
        ((integrator(Guse*mass*nu*r,r))**2.0*sigint)
    return vsp1, vsp2, zeta_A, zeta_B

#Richardson-Fairbairn VSP estimators:
def richfair_vsp(vz,Rkin,mskin):
    vsp1_RF = 1.0/(np.pi*2.0)*\
        np.sum(vz**4.0*mskin)/np.sum(mskin)
    vsp2_RF = 1.0/(np.pi*2.0)*\
        np.sum(vz**4.0*mskin*Rkin**2.0)/np.sum(mskin*Rkin**2.0)
    return vsp1_RF, vsp2_RF

###########################################################
#For Jeans modelling:
def sigp(r1,r2,nu,Sigfunc,M,beta,betaf,nupars,Mpars,betpars,\
         Mstar_rad,Mstar_prof,Mstar,G,rmin,rmax):
    #Calculate projected velocity dispersion profiles
    #given input *functions* nu(r); M(r); beta(r); betaf(r).
    #Also input is an array Mstar_prof(Mstar_rad) describing the 3D 
    #cumulative stellar mass profile. This should be normalised 
    #so that it peaks at 1.0. The total stellar mass is passed in Mstar.

    #Set up theta integration array:
    intpnts = 100
    thmin = 0.
    bit = 1.e-5
    thmax = np.pi/2.-bit
    th = np.linspace(thmin,thmax,intpnts)
    sth = np.sin(th)
    cth = np.cos(th)
    cth2 = cth**2.

    #Set up rint interpolation array: 
    rint = np.logspace(np.log10(rmin),np.log10(rmax),intpnts)

    #First calc sigr2(rint):
    sigr2 = np.zeros(len(rint))
    nur = nu(rint,nupars)
    betafunc = betaf(rint,betpars)
    for i in range(len(rint)):
        rq = rint[i]/cth
        if (Mstar > 0):
            Mq = M(rq,Mpars)+Mstar*np.interp(rq,Mstar_rad,Mstar_prof)
        else:
            Mq = M(rq,Mpars)
        nuq = nu(rq,nupars)
        betafuncq = betaf(rq,betpars)
        sigr2[i] = 1./nur[i]/rint[i]/betafunc[i] * \
            integrator(G*Mq*nuq*betafuncq*sth,th)
    
    #And now the sig_LOS projection: 
    Sig = Sigfunc(rint,nupars)
    sigLOS2 = np.zeros(len(rint))
    for i in range(len(rint)):
        rq = rint[i]/cth
        nuq = nu(rq,nupars)
        sigr2q = np.interp(rq,rint,sigr2,left=0,right=0)
        betaq = beta(rq,betpars)
        sigLOS2[i] = 2.0*rint[i]/Sig[i]*\
            integrator((1.0-betaq*cth2)*nuq*sigr2q/cth2,th)

    sigr2out = np.interp(r2,rint,sigr2,left=0,right=0)
    sigLOS2out = np.interp(r2,rint,sigLOS2,left=0,right=0)
    Sigout = np.interp(r1,rint,Sig,left=0,right=0)

    return sigr2out, Sigout, sigLOS2out

def sigp_vs(r1,r2,nu,Sigfunc,M,beta,betaf,nupars,Mpars,betpars,\
            Mstar_rad,Mstar_prof,Mstar,G,rmin,rmax):
    #Calculate projected velocity dispersion profiles
    #given input *functions* nu(r); M(r); beta(r); betaf(r).
    #Also input is an array Mstar_prof(Mstar_rad) describing the 3D
    #cumulative stellar mass profile. This should be normalised
    #so that it peaks at 1.0. The total stellar mass is passed in Mstar.
    #Finally, the routine calculates a dimensional version of the
    #fourth order "virial shape" parmaeters in Richardson & Fairbairn 2014
    #described in their equations 8 and 9.

    #Set up theta integration array:
    intpnts = 100
    thmin = 0.
    bit = 1.e-5
    thmax = np.pi/2.-bit
    th = np.linspace(thmin,thmax,intpnts)
    sth = np.sin(th)
    cth = np.cos(th)
    cth2 = cth**2.

    #Set up rint interpolation array:
    rintpnts = 100
    rint = np.logspace(np.log10(rmin),\
                       np.log10(rmax),rintpnts)

    #First calc sigr2(rint):
    sigr2 = np.zeros(len(rint))
    nur = nu(rint,nupars)
    betafunc = betaf(rint,betpars)
    for i in range(len(rint)):
        rq = rint[i]/cth
        if (Mstar > 0):
            Mq = M(rq,Mpars)+Mstar*np.interp(rq,Mstar_rad,Mstar_prof)
        else:
            Mq = M(rq,Mpars)
        nuq = nu(rq,nupars)
        betafuncq = betaf(rq,betpars)
        sigr2[i] = 1./nur[i]/rint[i]/betafunc[i] * \
            integrator(G*Mq*nuq*betafuncq*sth,th)

    #And now the sig_LOS projection:
    Sig = Sigfunc(rint,nupars)
    sigLOS2 = np.zeros(len(rint))
    for i in range(len(rint)):
        rq = rint[i]/cth
        nuq = nu(rq,nupars)
        sigr2q = np.interp(rq,rint,sigr2,left=0,right=0)
        betaq = beta(rq,betpars)
        sigLOS2[i] = 2.0*rint[i]/Sig[i]*\
            integrator((1.0-betaq*cth2)*nuq*sigr2q/cth2,th)

    #And now the dimensional fourth order "virial shape"
    #parameters:
    betar = beta(rint,betpars)
    Mr = M(rint,Mpars)+Mstar*np.interp(rint,Mstar_rad,Mstar_prof)
    vs1 = 2.0/5.0*integrator(nur*(5.0-2.0*betar)*sigr2*\
                             G*Mr*rint,rint)
    vs2 = 4.0/35.0*integrator(nur*(7.0-6.0*betar)*sigr2*\
                              G*Mr*rint**3.0,rint)

    sigr2out = np.interp(r2,rint,sigr2,left=0,right=0)
    sigLOS2out = np.interp(r2,rint,sigLOS2,left=0,right=0)
    Sigout = np.interp(r1,rint,Sig,left=0,right=0)

    return sigr2out, Sigout, sigLOS2out, vs1, vs2

def sigp_prop(r1,r2,nu,Sigfunc,M,beta,betaf,nupars,Mpars,betpars,\
              Mstar_rad,Mstar_prof,Mstar,G,rmin,rmax):
    #Calculate projected velocity dispersion profiles
    #given input *functions* nu(r); M(r); beta(r); betaf(r).
    #Also input is an array Mstar_prof(Mstar_rad) describing the 3D
    #cumulative stellar mass profile. This should be normalised
    #so that it peaks at 1.0. The total stellar mass is passed in Mstar.

    #Set up theta integration array:
    intpnts = 100
    thmin = 0.
    bit = 1.e-5
    thmax = np.pi/2.-bit
    th = np.linspace(thmin,thmax,intpnts)
    sth = np.sin(th)
    cth = np.cos(th)
    cth2 = cth**2.

    rint = np.logspace(np.log10(rmin),np.log10(rmax),intpnts)

    sigr2 = np.zeros(len(rint))
    nur = nu(rint,nupars)
    betafunc = betaf(rint,betpars)
    for i in range(len(rint)):
        rq = rint[i]/cth
        if (Mstar > 0):
            Mq = M(rq,Mpars)+Mstar*np.interp(rq,Mstar_rad,Mstar_prof)
        else:
            Mq = M(rq,Mpars)
        nuq = nu(rq,nupars)
        betafuncq = betaf(rq,betpars)
        sigr2[i] = 1./nur[i]/rint[i]/betafunc[i] * \
            integrator(G*Mq*nuq*betafuncq*sth,th)
 
    Sig = Sigfunc(rint,nupars)
    sigLOS2 = np.zeros(len(rint))
    sigpmr2 = np.zeros(len(rint))
    sigpmt2 = np.zeros(len(rint))
    for i in range(len(rint)):
        rq = rint[i]/cth
        nuq = nu(rq,nupars)
        sigr2q = np.interp(rq,rint,sigr2,left=0,right=0)
        betaq = beta(rq,betpars)
        sigLOS2[i] = 2.0*rint[i]/Sig[i]*\
            integrator((1.0-betaq*cth2)*nuq*sigr2q/cth2,th)
        sigpmr2[i] = 2.0*rint[i]/Sig[i]*\
            integrator((1.0-betaq+betaq*cth2)*nuq*sigr2q/cth2,th)
        sigpmt2[i] = 2.0*rint[i]/Sig[i]*\
            integrator((1.0-betaq)*nuq*sigr2q/cth2,th)

    sigr2out = np.interp(r2,rint,sigr2,left=0,right=0)
    sigLOS2out = np.interp(r2,rint,sigLOS2,left=0,right=0)
    sigpmr2out = np.interp(r2,rint,sigpmr2,left=0,right=0)
    sigpmt2out = np.interp(r2,rint,sigpmt2,left=0,right=0)
    Sigout = np.interp(r1,rint,Sig,left=0,right=0)
    
    return sigr2out, Sigout, sigLOS2out, sigpmr2out, sigpmt2out

def beta(r,betpars):
    bet0star = betpars[0]
    betinfstar = betpars[1]
    r0 = 10.**betpars[2]
    n = betpars[3]

    #Ensure stability at beta extremities:
    if (bet0star > 0.98): 
        bet0star = 0.98
    if (bet0star < -0.95):
        bet0star = -0.95
    if (betinfstar > 0.98):
        betinfstar = 0.98
    if (betinfstar < -0.95):
        betinfstar = -0.95
    bet0 = 2.0*bet0star / (1.0 + bet0star)
    betinf = 2.0*betinfstar / (1.0 + betinfstar)

    beta = bet0 + (betinf-bet0)*(1.0/(1.0 + (r0/r)**n))
    return beta

def betaf(r,betpars):
    bet0star = betpars[0]
    betinfstar = betpars[1]
    r0 = 10.**betpars[2]
    n = betpars[3]

    #Ensure stability at beta extremities:
    if (bet0star > 0.98):
        bet0star = 0.98
    if (bet0star < -0.95):
        bet0star = -0.95
    if (betinfstar > 0.98):
        betinfstar = 0.98
    if (betinfstar < -0.95):
        betinfstar = -0.95
    bet0 = 2.0*bet0star / (1.0 + bet0star)
    betinf = 2.0*betinfstar / (1.0 + betinfstar)

    betafn = r**(2.0*betinf)*((r0/r)**n+1.0)**(2.0/n*(betinf-bet0))
    
    return betafn


###########################################################
#For data binning:
def binthedata(R,ms,Nbin):
    #Nbin is the number of particles / bin:
    index = np.argsort(R)
    right_bin_edge = np.zeros(len(R))
    norm = np.zeros(len(R))
    cnt = 0
    jsum = 0

    for i in range(len(R)):
        if (jsum < Nbin):
            norm[cnt] = norm[cnt] + ms[index[i]]
            right_bin_edge[cnt] = R[index[i]]
            jsum = jsum + ms[index[i]]
        if (jsum >= Nbin):
            jsum = 0.0
            cnt = cnt + 1
    
    right_bin_edge = right_bin_edge[:cnt]
    norm = norm[:cnt]
    surfden = np.zeros(cnt)
    rbin = np.zeros(cnt)
    
    for i in range(len(rbin)):
        if (i == 0):
            surfden[i] = norm[i] / \
                (np.pi*right_bin_edge[i]**2.0)
            rbin[i] = right_bin_edge[i]/2.0
        else:
            surfden[i] = norm[i] / \
                (np.pi*right_bin_edge[i]**2.0-\
                 np.pi*right_bin_edge[i-1]**2.0)
            rbin[i] = (right_bin_edge[i]+right_bin_edge[i-1])/2.0
    surfdenerr = surfden / np.sqrt(Nbin)
    
    #Calculate the projected half light radius &
    #surface density integral:
    Rhalf, Menc_tot = surf_renorm(rbin,surfden)

    #And normalise the profile:
    surfden = surfden / Menc_tot
    surfdenerr = surfdenerr / Menc_tot

    return rbin, surfden, surfdenerr, Rhalf

def surf_renorm(rbin,surfden):
    #Calculate the integral of the surface density
    #so that it can then be renormalised.
    #Calcualte also Rhalf.
    ranal = np.linspace(0,10,np.int(5000))
    surfden_ranal = np.interp(ranal,rbin,surfden,left=0,right=0)
    Menc_tot = 2.0*np.pi*integrator(surfden_ranal*ranal,ranal)
    Menc_half = 0.0
    i = 3
    while (Menc_half < Menc_tot/2.0):
        Menc_half = 2.0*np.pi*\
            integrator(surfden_ranal[:i]*ranal[:i],ranal[:i])
        i = i + 1
    Rhalf = ranal[i-1]
    return Rhalf, Menc_tot


###########################################################
#For calculating confidence intervals: 
def calcmedquartnine(array):
    index = np.argsort(array,axis=0)
    median = array[index[np.int(len(array)/2.)]]
    sixlowi = np.int(16./100. * len(array))
    sixhighi = np.int(84./100. * len(array))
    ninelowi = np.int(2.5/100. * len(array))
    ninehighi = np.int(97.5/100. * len(array))
    nineninelowi = np.int(0.15/100. * len(array))
    nineninehighi = np.int(99.85/100. * len(array))

    sixhigh = array[index[sixhighi]]
    sixlow = array[index[sixlowi]]
    ninehigh = array[index[ninehighi]]
    ninelow = array[index[ninelowi]]
    nineninehigh = array[index[nineninehighi]]
    nineninelow = array[index[nineninelowi]]

    return median, sixlow, sixhigh, ninelow, ninehigh,\
        nineninelow, nineninehigh


###########################################################
#For fitting the surface brightness:
def multiplumden(r,pars):
    Mpars = pars[0:np.int(len(pars)/2.0)]
    apars = pars[np.int(len(pars)/2.0):len(pars)]
    nplum = len(Mpars)
    multplum = np.zeros(len(r))
    for i in range(len(Mpars)):
        if (multimode == 'seq'):
            if (i == 0):
                aparsu = apars[0]
            else:
                aparsu = apars[i] + apars[i-1]
        else:
            aparsu = apars[i]
        multplum = multplum + \
            3.0*Mpars[i]/(4.*np.pi*aparsu**3.)*\
            (1.0+r**2./aparsu**2.)**(-5./2.)
    return multplum

def multiplumsurf(r,pars):
    Mpars = pars[0:np.int(len(pars)/2.0)]
    apars = pars[np.int(len(pars)/2.0):len(pars)]
    nplum = len(Mpars)
    multplum = np.zeros(len(r))
    for i in range(len(Mpars)):
        if (multimode == 'seq'):
            if (i == 0):
                aparsu = apars[0]
            else:
                aparsu = apars[i] + apars[i-1]
        else:
            aparsu = apars[i]
        multplum = multplum + \
            Mpars[i]*aparsu**2.0 / \
            (np.pi*(aparsu**2.0+r**2.0)**2.0)
    return multplum

def multiplumdlnrhodlnr(r,pars):
    Mpars = pars[0:np.int(len(pars)/2.0)]
    apars = pars[np.int(len(pars)/2.0):len(pars)]
    nplum = len(Mpars)
    multplumden = np.zeros(len(r))
    multplumdden = np.zeros(len(r))
    for i in range(len(Mpars)):
        if (multimode == 'seq'):
            if (i == 0):
                aparsu = apars[0]
            else:
                aparsu = apars[i] + apars[i-1]
        else:
            aparsu = apars[i]
        multplumden = multplumden + \
            3.0*Mpars[i]/(4.*np.pi*aparsu**3.)*\
            (1.0+r**2./aparsu**2.)**(-5./2.)
        multplumdden = multplumdden - \
            15.0*Mpars[i]/(4.*np.pi*aparsu**3.)*\
            r/aparsu**2.*(1.0+r**2./aparsu**2.)**(-7./2.)
    return multplumdden*r/multplumden

def multiplummass(r,pars):
    Mpars = pars[0:np.int(len(pars)/2.0)]
    apars = pars[np.int(len(pars)/2.0):len(pars)]
    nplum = len(Mpars)
    multplum = np.zeros(len(r))
    for i in range(len(Mpars)):
        if (multimode == 'seq'):
            if (i == 0):
                aparsu = apars[0]
            else:
                aparsu = apars[i] + apars[i-1]
        else:
            aparsu = apars[i]
        multplum = multplum + \
            Mpars[i]*r**3./(r**2.+aparsu**2.)**(3./2.)
    return multplum

def threeplumsurf(r,M1,M2,M3,a1,a2,a3):
    return multiplumsurf(r,[M1,M2,M3,\
                            a1,a2,a3])
def threeplumden(r,M1,M2,M3,a1,a2,a3):
    return multiplumden(r,[M1,M2,M3,\
                           a1,a2,a3])
def threeplummass(r,M1,M2,M3,a1,a2,a3):
    return multiplummass(r,[M1,M2,M3,\
                            a1,a2,a3])

def Rhalf_func(M1,M2,M3,a1,a2,a3):
    #Calculate projected half light radius for
    #the threeplum model:
    ranal = np.logspace(-3,1,np.int(500))
    Mstar_surf = threeplumsurf(ranal,M1,M2,M3,a1,a2,a3)

    Menc_half = 0.0
    i = 3
    while (Menc_half < (M1+M2+M3)/2.0):
        Menc_half = 2.0*np.pi*\
            integrator(Mstar_surf[:i]*ranal[:i],ranal[:i])
        i = i + 1
    Rhalf = ranal[i-1]
    return Rhalf


###########################################################
#For fitting the velocity distribution in each bin [no errors]:
def monte(func,a,b,n):
    #Function to perform fast 1D Monte-Carlo integration
    #for convolution integrals:
    xrand = np.random.uniform(a,b,n)
    integral = func(xrand).sum()
    return (b-a)/np.float(n)*integral

def velpdf_noerr(vz,theta):
    vzmean = theta[0]
    alp = theta[1]
    bet = theta[2]
    backamp = theta[3]
    backmean = theta[4]
    backsig = theta[5]

    pdf = (1.0-backamp)*bet/(2.0*alp*gamma(1.0/bet))*\
        np.exp(-(np.abs(vz-vzmean)/alp)**bet) + \
        backamp/(np.sqrt(2.0*np.pi)*backsig)*\
        np.exp(-0.5*(vz-backmean)**2.0/backsig**2.0)
    return pdf

#For fitting the velocity distribution in each bin [fast]
#Uses an approximation to the true convolution integral.
def velpdffast(vz,vzerr,theta):
    vzmean = theta[0]
    bet = theta[2]
    fgamma = gamma(1.0/bet)/gamma(3.0/bet)
    alp = np.sqrt(theta[1]**2.0+vzerr**2.0*fgamma)    
    backamp = theta[3]
    backmean = theta[4]
    backsig = np.sqrt(theta[5]**2.0 + vzerr**2.0)
    
    pdf = (1.0-backamp)*bet/(2.0*alp*gamma(1.0/bet))*\
        np.exp(-(np.abs(vz-vzmean)/alp)**bet) + \
        backamp/(np.sqrt(2.0*np.pi)*backsig)*\
        np.exp(-0.5*(vz-backmean)**2.0/backsig**2.0)
    return pdf
    
def velpdf_func(vz,vzerr,vzint,theta):
    #Inner integral function for convolving
    #velpdf with a Gaussian error PDF. Change
    #this function to implement non-Gaussian
    #errors.
    vzmean = theta[0]
    alp = theta[1]
    bet = theta[2]
    backamp = theta[3]
    backmean = theta[4]
    backsig = theta[5]
    pdf = (1.0-backamp)*bet/(2.0*alp*gamma(1.0/bet))*\
                np.exp(-(np.abs(vzint-vzmean)/alp)**bet)*\
                1.0/(np.sqrt(2.0*np.pi)*vzerr)*\
                np.exp(-0.5*(vz-vzint)**2.0/vzerr**2.0)+\
                backamp/(np.sqrt(2.0*np.pi)*backsig)*\
                np.exp(-0.5*(vzint-backmean)**2.0/backsig**2.0)*\
                1.0/(np.sqrt(2.0*np.pi)*vzerr)*\
                np.exp(-0.5*(vz-vzint)**2.0/vzerr**2.0)
    return pdf

#For fitting the velocity distribution in each bin with
#full (expensive) convolution integral:
def velpdf(vz,vzerr,theta):
    #Generalised Gaussian + Gaussian convolved with
    #vzerr, assuming Gaussian errors:
    vzmean = theta[0]
    sig = vztwo_calc(theta)
    vzlow = -sig*10+vzmean
    vzhigh = sig*10+vzmean
    if (type(vz) == np.ndarray):
        pdf = np.zeros(len(vz))
        for i in range(len(vz)):
            pdf_func = lambda vzint : velpdf_func(vz[i],\
                vzerr[i],vzint,theta)
            pdf[i] = quad(pdf_func,vzlow,vzhigh)[0]
    else:
        pdf_func = lambda vzint : velpdf_func(vz,\
            vzerr,vzint,theta)
        pdf = quad(pdf_func,vzlow,vzhigh)[0]
    return pdf

def velpdfmonte(vz,vzerr,theta):
    #Generalised Gaussian + Gaussian convolved with
    #vzerr, assuming Gaussian errors:
    npnts = np.int(500)
    vzmean = theta[0]
    sig = vztwo_calc(theta)
    vzlow = -sig*10+vzmean
    vzhigh = sig*10+vzmean
    if (type(vz) == np.ndarray):
        pdf = np.zeros(len(vz))
        for i in range(len(vz)):
            pdf_func = lambda vzint : velpdf_func(vz[i],\
                vzerr[i],vzint,theta)
            pdf[i] = monte(pdf_func,vzlow,vzhigh,npnts)
    else:
        pdf_func = lambda vzint : velpdf_func(vz,\
            vzerr,vzint,theta)
        pdf = monte(pdf_func,vzlow,vzhigh,npnts)
    return pdf

def vztwo_calc(theta):
    #Calculate <vlos^2>^(1/2) from
    #generalised Gaussian parameters:
    alp = theta[1]
    bet = theta[2]
    return np.sqrt(alp**2.0*gamma(3.0/bet)/gamma(1.0/bet))
    
def vzfour_calc(theta):
    #Calculate <vlos^4> from
    #generalised Gaussian parameters:
    alp = theta[1]
    bet = theta[2]
    sig = vztwo_calc(theta)
    kurt = gamma(5.0/bet)*gamma(1.0/bet)/(gamma(3.0/bet))**2.0
    return kurt*sig**4.0
    
def kurt_calc(theta):
    #Calculate kurtosis from generalised
    #Gaussian parameters:
    alp = theta[1]
    bet = theta[2]
    kurt = gamma(5.0/bet)*gamma(1.0/bet)/(gamma(3.0/bet))**2.0
    return kurt

def vzfourfunc(ranal,rbin,vzfourbin):
    #Interpolate and extrapolate
    #vzfour(R) over and beyond the data:
    vzfour = np.interp(ranal,rbin,vzfourbin)
    return vzfour
    
#For calculating the Likelihood from the vsp array:
def vsppdf_calc(vsp):
    #First bin the data:
    nbins = 50
    bins_plus_one = np.linspace(np.min(vsp),np.max(vsp),nbins+1)
    bins = np.linspace(np.min(vsp),np.max(vsp),nbins)
    vsp_pdf, bins_plus_one = np.histogram(vsp, bins=bins_plus_one)
    vsp_pdf = vsp_pdf / np.max(vsp_pdf)
    binsout = bins[vsp_pdf > 0]
    vsp_pdfout = vsp_pdf[vsp_pdf > 0]
    return binsout, vsp_pdfout

def vsp_pdf(vsp,bins,vsp_pdf):
    return np.interp(vsp,bins,vsp_pdf,left=0,right=0)

