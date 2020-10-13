### SUBROUTINES & FUNCTIONS HERE ###
def vmax_func(M200,c200,h):
    oden = 200.0
    Guse = G*Msun/kpc
    r200=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    
    #This from Sigad et al. 2000 (via Schneider et al. 2017):
    vmax = 0.465*np.sqrt(Guse*M200/r200)/\
           np.sqrt(1.0/c200*np.log(1.0+c200)-(1.0+c200)**(-1.0))
    return vmax/kms

def radius_dsph(s, b, distance):
    return np.sqrt((distance * np.sin(b))**2. + s*s)

def integrand(s, b, distance, Mpars):
    value = np.sin(b) * rho(np.array([radius_dsph(s, b, distance)]), Mpars)**2
    return value

def Xfunc(s):
    #Needed for Wyn Jfac analytic calculation:
    if (s < 1.0):
        Xs = 1.0/np.sqrt(1.0-s**2.0)*np.arccosh(1.0/s)
    elif (s == 1.0):
        Xs = 1.0
    else:
        Xs = 1.0/np.sqrt(s**2.0-1.0)*np.arccos(1.0/s)
    return Xs                

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
    J_max = si.dblquad(integrand,b_min,b_max,s_min_bound,\
                       s_max_bound,args=(distance,Mpars),\
                       epsabs=Acc_arr,epsrel=Acc_arr)
    J_max = J_max[0]*kpccm*2.*np.pi*Msunkpc3toGeVcm3**2.0

    #Error checking:
    if (J_max == np.inf):
        print 'Argh! Infinite J_max!! Bye bye...'
        sys.exit(0)
        
    if (J_max < 0):
        print 'Argh! Negative J_max!! Bye bye...'
        sys.exit(0)

    return J_max  # in GeV2.cm-5

def get_J_NFW(Mpars, distance, r_max):
    M200 = Mpars[0]
    c = Mpars[1]
    rc = Mpars[2]
    n = Mpars[3]
    
    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=oden*c**3.*gcon/3.
    rv=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    rs=rv/c
    rhos=rhocrit*deltachar
                                       
    D = distance
    theta = np.arcsin(r_max/D)
    y = D*theta/rs
    Del = np.sqrt(1.0-y**2.0)
    J_max = np.pi*rhos**2.0*rs**3.0/\
            (3.0*D**2.0*Del**4.0)*(2.0*y*\
            (7.0*y-4.0*y**3.0+\
             3.0*np.pi*Del**4.0)+\
            6.0*(2.0*Del**6.0-2.0*Del**2.0-\
                 y**4.0)*\
            Xfunc(y))*kpccm*Msunkpc3toGeVcm3**2.
    
    return J_max

def cosmo_cfunc_DM(M200,h):
    #From Dutton & Maccio 2014
    c = 10.**(0.905 - 0.101 * (np.log10(M200*h)-12.))
    return c
    
def mass_hm(mWDM):
    """
    Returns the WDM half-mode mass.
    """

    # WDM mass function definitions
    h0       = cosmo.Hz(0)/100
    rho_bar0 = cosmo.rho_m(0) * h0**2
    OMEGA_WDM = cosmo.Om(0)
    MU        = 1.12 
    alpha_hm  = 49.0 * mWDM**-1.11 * \
                (OMEGA_WDM/0.25)**0.11 * (h0/0.7)**1.22 /h0  # in KPC
    lambda_hm = 2*np.pi * alpha_hm * (2**(MU/5.) - 1)**(-1./2/MU)
    return  4*np.pi/3 * rho_bar0 * (lambda_hm/2.)**3


def cosmo_cfunc_DJ(m200c, h, mWDM):
    """
    WDM mass-concentration relation.
    Requires as input masses defined in 200c system in units of MSUN.
    """

    # WDM concentration from Schneider+ 2012
    # parameters to convert cCDM to cWDM via their eq. 39
    GAMMA1 = 15.
    GAMMA2 = 0.3

    h0       = cosmo.Hz(0)/100 
    rho_bar0 = cosmo.rho_m(0) * h0**2 
    c200c = colossus_cNFW(m200c, '200c', 0.0, model='diemer15')
    m1m_div_h, r1m_div_h, c1m = \
            changeMassDefinition(m200c/h0, c200c, 0.0, '200c', '1m')
    m1m = m1m_div_h * h0
    return c200c*(1 + GAMMA1*mass_hm(mWDM)/m1m)**(-GAMMA2)

def corenfw_den(r,M200,c,oden,tSF,Rhalf,rhocrit,eta,kappa):
    #If on, this defines tdyn at rc rather than rs:
    newswitch = 'no'

    #Assumes input arrays in Msun, kpc, Gyrs:
    G = 6.67e-11
    kpc = 3.086e19
    Msun = 1.989e30
    Gyr = 365.*24.*60.*60.*1e9

    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=oden*c**3.*gcon/3.
    rv=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    rs=rv/c
    rhos=rhocrit*deltachar
    rhoanal = rhos/((r/rs)*(1.+(r/rs))**2.)
    manal = M200 * gcon * (np.log(1.0 + r/rs)-r/rs/(1.0+r/rs))

    if (newswitch == 'no'):
        Mrs = M200 * gcon * (np.log(2.0)-0.5)
        tdyn_rs = \
            2.*np.pi*np.sqrt((rs*kpc)**3./G/(Mrs*Msun))/Gyr
    else:
        rc = eta * Rhalf
        Mrc = M200 * gcon * (np.log(1.0 + rc/rs)-rc/rs/(1.0+rc/rs))
        tdyn_rs = \
            2.*np.pi*np.sqrt((rc*kpc)**3./G/(Mrc*Msun))/Gyr

    rc = eta * Rhalf
    x = r/rc
    f = (np.exp(x) - np.exp(-x))/(np.exp(x)+np.exp(-x))
    xx = tSF/tdyn_rs * kappa
    n = (np.exp(xx) - np.exp(-xx))/(np.exp(xx)+np.exp(-xx))
    my_manal = manal*f**n

    my_rhoanal = rhoanal*f**n + \
        1.0/(4.*np.pi*r**2.*rc)*manal*(1.0-f**2.)*n*f**(n-1.0)

    return my_rhoanal

def corenfw_mass(r,M200,c,oden,tSF,Rhalf,rhocrit,eta,kappa):
    newswitch = 'no'

    #Assumes input arrays in Msun, kpc:
    G = 6.67e-11
    kpc = 3.086e19
    Msun = 1.989e30
    Gyr = 365.*24.*60.*60.*1e9

    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=oden*c**3.*gcon/3.
    rv=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)
    rs=rv/c
    rhos=rhocrit*deltachar
    rhoanal = rhos/((r/rs)*(1.+(r/rs))**2.)
    manal = M200 * gcon * (np.log(1.0 + r/rs)-r/rs/(1.0+r/rs))

    if (newswitch == 'no'):
        Mrs = M200 * gcon * (np.log(2.0)-0.5)
        tdyn_rs = \
            2.*np.pi*np.sqrt((rs*kpc)**3./G/(Mrs*Msun))/Gyr
    else:
        rc = eta * Rhalf
        Mrc = M200 * gcon * (np.log(1.0 + rc/rs)-rc/rs/(1.0+rc/rs))
        tdyn_rs = \
            2.*np.pi*np.sqrt((rc*kpc)**3./G/(Mrc*Msun))/Gyr

    rc = eta * Rhalf
    x = r/rc
    f = (np.exp(x) - np.exp(-x))/(np.exp(x)+np.exp(-x))
    xx = tSF/tdyn_rs * kappa
    n = (np.exp(xx) - np.exp(-xx))/(np.exp(xx)+np.exp(-xx))
    my_manal = manal*f**n

    my_rhoanal = rhoanal*f**n + \
        1.0/(4.*np.pi*r**2.*rc)*manal*(1.0-f**2.)*n*f**(n-1.0)

    return my_manal

def rhoNFW(r,rhos,rs):
    return rhos/((r/rs)*(1.+(r/rs))**2.)

def sidm_novel(rc,M200,c,oden,rhocrit,rt):
#    GammaX = 0.4/(1e9*year)
    GammaX = 2.486/(1e9*year)
    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=oden*c**3.*gcon/3.
    rv=(3./4.*M200/(np.pi*oden*rhocrit))**(1./3.)*kpc
    rs=rv/c
    rhos=rhocrit*deltachar*Msun/kpc**3.0

    rhorc = rhoNFW(rc*kpc,rhos,rs)
    rout = rs*500.0
#    if (rt > 0):
#        rout = rt*kpc
    r = np.logspace(np.log10(rc*kpc),np.log10(rout),5000)
    rho = rhoNFW(r,rhos,rs)
    mass = M200*Msun*gcon*(np.log(1.0 + r/rs)-r/rs/(1.0+r/rs))
    sigvtworc = G/rhorc*simps(mass*rho/r**2.0,r)
    sigvrc = np.sqrt(sigvtworc)

#    sigm = np.sqrt(np.pi)*GammaX/(4.0*rhorc*sigvrc)
    sigm = np.sqrt(np.pi)*GammaX/(4.0*rhorc*sigvrc)*(rc*kpc/rs)**1.3
    return sigm*100.0**2.0/1000.0

def sidm_novel3(rc,M200,c,oden,rhocrit,GammaX):
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
    sigvtworc = Guse/rhorc*simps(mass*rho/r**2.0,r)
    sigvrc = np.sqrt(sigvtworc)

    sigm = np.sqrt(np.pi)*GammaX/(4.0*rhorc*rho_unit*sigvrc)
    return sigm*100.0**2.0/1000.0

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

def sigp(r1,r2,nu,Sigfunc,M,beta,betaf,nupars,Mpars,betpars,\
         Mstar_rad,Mstar_prof,Mstar,G):
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
            Mstar_rad,Mstar_prof,Mstar,G):
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
              Mstar_rad,Mstar_prof,Mstar,G):
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
    if (linbetr0 == 'yes'):
        r0 = betpars[2]
    else:
        r0 = 10.**betpars[2]
    n = betpars[3]

    if (betmode == 'betstar'):
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
    else:
        bet0 = bet0star
        betinf = betinfstar

    beta = bet0 + (betinf-bet0)*(1.0/(1.0 + (r0/r)**n))
    return beta

def betaf(r,betpars):
    bet0star = betpars[0]
    betinfstar = betpars[1]
    if (linbetr0 == 'yes'):
        r0 = betpars[2]
    else:
        r0 = 10.**betpars[2]
    n = betpars[3]
    if (rotation == 'yes'):
        Arot = betpars[4]

    if (betmode == 'betstar'):
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
    else:
        bet0 = bet0star
        betinf = betinfstar

    betafn = r**(2.0*betinf)*((r0/r)**n+1.0)**(2.0/n*(betinf-bet0))
    
    if (rotation == 'yes'):
        betafn = betafn * np.exp(-2.0*Arot*r/Rhalf)
    
    return betafn

def multiplumden(r,pars):
    Mpars = pars[0:len(pars)/2]
    apars = pars[len(pars)/2:len(pars)]
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
    Mpars = pars[0:len(pars)/2]
    apars = pars[len(pars)/2:len(pars)]
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
    Mpars = pars[0:len(pars)/2]
    apars = pars[len(pars)/2:len(pars)]
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
    Mpars = pars[0:len(pars)/2]
    apars = pars[len(pars)/2:len(pars)]
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

def gambindlnrhodlnr(r,gampars,binpars):
    dden = np.zeros(len(r))
    j = 0
    for i in range(len(r)):
        if (r[i] < binpars[j]):
            gamuse = np.abs(gampars[j+1])
            dden[i] = gamuse
        else:
            if (j < len(binpars)-1):
                j = j + 1
            gamuse = np.abs(gampars[j+1])
            dden[i] = gamuse
    return -dden

def gambinden(r,gampars,binpars):
    den = np.zeros(len(r))
    j = 0
    rho_j = gampars[0]
    for i in range(len(r)):
        if (r[i] < binpars[j]):
            gamuse = np.abs(gampars[j+1])
            den[i] = rho_j * (r[i]/binpars[j])**(-gamuse)
        else:
            if (j < len(binpars)-1):
                j = j + 1
                gamuse = np.abs(gampars[j+1])
                rho_j = rho_j * (binpars[j]/binpars[j-1])**(-gamuse) 
            den[i] = rho_j * (r[i]/binpars[j])**(-gamuse)
    return den

def mpower(rleft,rright,rho0,rc,gam):
    #*** WARNING *** This will not work for gamma > 3.0
    mypow = 4.0*np.pi*rho0*rc**3.0/\
        (3.0-gam)*\
        ((rright/rc)**(3.0-gam)-\
         (rleft/rc)**(3.0-gam))
    return mypow

def gambinmass(r,gampars,binpars):
    mass = np.zeros(len(r))
    j = 0
    rho_j = gampars[0]
    bin_j_min_one = 0.0
    mass_j_min_one = 0.0
    for i in range(len(r)):
        if (r[i] < binpars[j]):
            gamuse = np.abs(gampars[j+1])
            mass[i] = mass_j_min_one + \
                mpower(bin_j_min_one,r[i],rho_j,binpars[j],gamuse)
        else:
            if (j < len(binpars)-1):
                gamuse = np.abs(gampars[j+1])
                mass_j_min_one = mass_j_min_one + \
                    mpower(bin_j_min_one,binpars[j],rho_j,\
                           binpars[j],gamuse)
                j = j + 1
                gamuse = np.abs(gampars[j+1])
                bin_j_min_one = binpars[j-1]
                rho_j = rho_j * (binpars[j]/binpars[j-1])**(-gamuse)
            mass[i] = mass_j_min_one + \
                mpower(bin_j_min_one,r[i],rho_j,binpars[j],gamuse)
    return mass

def lnprior_set_single(theta,n_betpars,bet0min,bet0max,betinfmin,betinfmax,\
                       betr0min,betr0max,betnmin,betnmax,\
                       nu_components,nupars_min,nupars_max,\
                       n_mpars,logM200low,logM200high,\
                       clow,chigh,logrclow,logrchigh,\
                       nlow,nhigh,logrtlow,logrthigh,\
                       dellow,delhigh,\
                       Mstar_min,Mstar_max):

    ndims = len(theta)
    minarr = np.zeros(ndims)
    maxarr = np.zeros(ndims)
    minarr[0] = bet0min
    maxarr[0] = bet0max
    minarr[1] = betinfmin
    maxarr[1] = betinfmax
    minarr[2] = betr0min
    maxarr[2] = betr0max
    minarr[3] = betnmin
    maxarr[3] = betnmax

    if (cosmo_cprior == 'yes'):
        M200u = 10.0**theta[n_betpars+nu_components*2]
        clowu = 10.0**(np.log10(cosmo_cfunc(M200u,h))-0.1)
        chighu = 10.0**(np.log10(cosmo_cfunc(M200u,h))+0.1)
    else:
        clowu = clow
        chighu = chigh

    minarr[n_betpars:n_betpars+nu_components*2] = nupars_min
    maxarr[n_betpars:n_betpars+nu_components*2] = nupars_max
    minarr[n_betpars+nu_components*2] = logM200low
    maxarr[n_betpars+nu_components*2] = logM200high
    minarr[n_betpars+nu_components*2+1] = clowu
    maxarr[n_betpars+nu_components*2+1] = chighu
    minarr[n_betpars+nu_components*2+2] = logrclow
    maxarr[n_betpars+nu_components*2+2] = logrchigh
    minarr[n_betpars+nu_components*2+3] = nlow
    maxarr[n_betpars+nu_components*2+3] = nhigh
    minarr[n_betpars+nu_components*2+4] = logrtlow
    maxarr[n_betpars+nu_components*2+4] = logrthigh
    minarr[n_betpars+nu_components*2+5] = dellow
    maxarr[n_betpars+nu_components*2+5] = delhigh

    minarr[ndims-1] = Mstar_min
    maxarr[ndims-1] = Mstar_max

    if all(minarr < theta < maxarr for minarr,theta,maxarr in \
               zip(minarr,theta,maxarr)):
        return 0.0
    return -np.inf

def lnprior_set_single_rot(theta,n_betpars,bet0min,bet0max,betinfmin,betinfmax,\
                           betr0min,betr0max,betnmin,betnmax,\
                           Arotmin,Arotmax,\
                           nu_components,nupars_min,nupars_max,\
                           n_mpars,logM200low,logM200high,\
                           clow,chigh,logrclow,logrchigh,\
                           nlow,nhigh,logrtlow,logrthigh,\
                           dellow,delhigh,\
                           Mstar_min,Mstar_max):
    
    ndims = len(theta)
    minarr = np.zeros(ndims)
    maxarr = np.zeros(ndims)
    minarr[0] = bet0min
    maxarr[0] = bet0max
    minarr[1] = betinfmin
    maxarr[1] = betinfmax
    minarr[2] = betr0min
    maxarr[2] = betr0max
    minarr[3] = betnmin
    maxarr[3] = betnmax
    minarr[4] = Arotmin
    maxarr[4] = Arotmax

    if (cosmo_cprior == 'yes'):
        M200u = 10.0**theta[n_betpars+nu_components*2]
        clowu = 10.0**(np.log10(cosmo_cfunc(M200u,h))-0.1)
        chighu = 10.0**(np.log10(cosmo_cfunc(M200u,h))+0.1)
    else:
        clowu = clow
        chighu = chigh

    minarr[n_betpars:n_betpars+nu_components*2] = nupars_min
    maxarr[n_betpars:n_betpars+nu_components*2] = nupars_max
    minarr[n_betpars+nu_components*2] = logM200low
    maxarr[n_betpars+nu_components*2] = logM200high
    minarr[n_betpars+nu_components*2+1] = clowu
    maxarr[n_betpars+nu_components*2+1] = chighu
    minarr[n_betpars+nu_components*2+2] = logrclow
    maxarr[n_betpars+nu_components*2+2] = logrchigh
    minarr[n_betpars+nu_components*2+3] = nlow
    maxarr[n_betpars+nu_components*2+3] = nhigh
    minarr[n_betpars+nu_components*2+4] = logrtlow
    maxarr[n_betpars+nu_components*2+4] = logrthigh
    minarr[n_betpars+nu_components*2+5] = dellow
    maxarr[n_betpars+nu_components*2+5] = delhigh
    
    minarr[ndims-1] = Mstar_min
    maxarr[ndims-1] = Mstar_max
    
    if all(minarr < theta < maxarr for minarr,theta,maxarr in \
           zip(minarr,theta,maxarr)):
        return 0.0
    return -np.inf

def lnprior_set_split(theta,n_betpars,bet0min,bet0max,betinfmin,betinfmax,\
                       betr0min,betr0max,betnmin,betnmax,\
                       nu_components,nupars1_min,nupars1_max,\
                       nupars2_min,nupars2_max,\
                       n_mpars,log_rho0_min,log_rho0_max,\
                       gam_min,gam_max,Mstar_min,Mstar_max):

    #Deal with negative gamma:
    thetau = np.array(theta)
    thetau[n_betpars*2+nu_components*2*2+1:len(thetau)-2] = \
        np.abs(theta[n_betpars*2+nu_components*2*2+1:len(thetau)-2])
    gam_minz = gam_min
    if (gam_minz < 0):
        gam_minz = 0.0

    #Regularize?
    if (regprior == 'no'):
        gam_minu = gam_minz
        gam_maxu = gam_max
    else:
        gam_minu = np.zeros(n_mpars-1)
        gam_maxu = np.zeros(n_mpars-1)
        gam_minu[0] = gam_minz
        gam_maxu[0] = gam_max
        gam_minu[1:n_mpars-1] = \
            thetau[n_betpars*2+nu_components*2*2+1:len(thetau)-2]-0.01
        gam_maxu[1:n_mpars-1] = \
            thetau[n_betpars*2+nu_components*2*2+1:len(thetau)-2]+1.0
        for i in range(1,n_mpars-1):
            if (gam_minu[i] < gam_minz):
                gam_minu[i] = gam_minz
            if (gam_maxu[i] > gam_max):
                gam_maxu[i] = gam_max

    ndims = len(thetau)
    minarr = np.zeros(ndims)
    maxarr = np.zeros(ndims)

    minarr[0] = bet0min
    maxarr[0] = bet0max
    minarr[1] = betinfmin
    maxarr[1] = betinfmax
    minarr[2] = betr0min
    maxarr[2] = betr0max
    minarr[3] = betnmin
    maxarr[3] = betnmax

    minarr[4] = bet0min
    maxarr[4] = bet0max
    minarr[5] = betinfmin
    maxarr[5] = betinfmax
    minarr[6] = betr0min
    maxarr[6] = betr0max
    minarr[7] = betnmin
    maxarr[7] = betnmax

    minarr[n_betpars*2:n_betpars*2+nu_components*2] = nupars1_min
    maxarr[n_betpars*2:n_betpars*2+nu_components*2] = nupars1_max

    minarr[n_betpars*2+nu_components*2:\
               n_betpars*2+nu_components*2+nu_components*2] = \
               nupars2_min
    maxarr[n_betpars*2+nu_components*2:\
               n_betpars*2+nu_components*2+nu_components*2] = \
               nupars2_max

    minarr[n_betpars*2+nu_components*2*2] = log_rho0_min
    maxarr[n_betpars*2+nu_components*2*2] = log_rho0_max

    minarr[n_betpars*2+nu_components*2*2+1:\
               n_betpars*2+nu_components*2*2+n_mpars] = gam_minu
    maxarr[n_betpars*2+nu_components*2*2+1:\
               n_betpars*2+nu_components*2*2+n_mpars] = gam_maxu

    minarr[ndims-1] = Mstar_min
    maxarr[ndims-1] = Mstar_max

    if all(minarr < thetau < maxarr for minarr,thetau,maxarr in \
               zip(minarr,thetau,maxarr)):
        return 0.0
    return -np.inf

def lnprob_single(theta, x1, x2, y1, y1err, y2, y2err):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x1, x2, y1, y1err, y2, y2err)

def lnprob_single_vs(theta, x1, x2, y1, y1err, \
                     y2, y2err, y3, y3err, y4, y4err):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x1, x2, y1, y1err, \
                       y2, y2err, y3, y3err, y4, y4err)

def lnprob_single_prop(theta, x1, x2, y1, y1err, \
                       y2, y2err, y3, y3err, y4, y4err):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x1, x2, y1, y1err, \
                       y2, y2err, y3, y3err, y4, y4err)

def lnprob_split(theta, x1, x2, x3, x4, y1, y1err, y2, y2err, \
                 y3, y3err, y4, y4err):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x1, x2, x3, x4, y1, y1err, \
                       y2, y2err, y3, y3err, \
                       y4, y4err)

def lnprob_split_vs(theta, x1, x2, x3, x4, y1, y1err, y2, y2err, \
                    y3, y3err, y4, y4err, y5, y5err, \
                    y6, y6err, y7, y7err, y8, y8err):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x1, x2, x3, x4, y1, y1err, \
                       y2, y2err, y3, y3err, \
                       y4, y4err, y5, y5err, \
                       y6, y6err, y7, y7err, \
                       y8, y8err)

def lnlike_single(theta, x1, x2, y1, y1err, y2, y2err):
    betpars = theta[0:n_betpars]
    nupars = theta[n_betpars:n_betpars+nu_components*2]
    Mpars = theta[n_betpars+nu_components*2:\
                  n_betpars+nu_components*2+n_mpars]
    Mstar = theta[n_betpars+nu_components*2+n_mpars]
    nuparsu = np.array(nupars)
    Mparsu = np.array(Mpars)
    Mparsu[0] = 10.**Mpars[0]
    if (use_rclinear == 'no'):
        Mparsu[2] = 10.**Mpars[2]
    if (use_rtlinear == 'no'):
        Mparsu[4] = 10.**Mpars[4]

    sigr2, Sig, sigLOS2 = \
           sigp_fit(x1,x2,nuparsu,Mparsu,betpars,Mstar)
    
    model1 = Sig
    model2 = np.sqrt(sigLOS2)/1000.
    
    inv_sigma2_1 = 1.0/y1err**2
    inv_sigma2_2 = 1.0/y2err**2

    lnlike_out = -0.5*(np.sum((y1-model1)**2*inv_sigma2_1)+\
                       np.sum((y2-model2)**2*inv_sigma2_2))

    if (SIDM_Draco_priors == 'cosmo_c'):
        #Add the conc. to the likelihood function
        #as a Gaussian in logspace:
        sigc = 0.1
        M200 = Mparsu[0]
        log_conc = np.log10(Mparsu[1])
        log_cmean = np.log10(cosmo_cfunc(M200,h))
        lnlike_out = lnlike_out * \
                     np.exp(-(log_conc-log_cmean)**2.0/(2.0*sigc**2.0))
    
    if (lnlike_out != lnlike_out):
        print sigLOS2
        print M(x1,Mpars)
        sys.exit(0)

    return lnlike_out

def lnlike_single_vs(theta, x1, x2, y1, y1err, y2, y2err, \
                     y3, y3err, y4, y4err):
    betpars = theta[0:n_betpars]
    nupars = theta[n_betpars:n_betpars+nu_components*2]
    Mpars = theta[n_betpars+nu_components*2:\
                  n_betpars+nu_components*2+n_mpars]
    Mstar = theta[n_betpars+nu_components*2+n_mpars]
    nuparsu = np.array(nupars)
    Mparsu = np.array(Mpars)
    Mparsu[0] = 10.**Mpars[0]
    if (use_rclinear == 'no'):
        Mparsu[2] = 10.**Mpars[2]
    if (use_rtlinear == 'no'):
        Mparsu[4] = 10.**Mpars[4]

    sigr2, Sig, sigLOS2, vs1, vs2 = \
        sigp_fit_vs(x1,x2,nuparsu,Mparsu,betpars,Mstar)

    model1 = Sig
    model2 = np.sqrt(sigLOS2)/1000.
    model3 = vs1/1.0e12
    model4 = vs2/1.0e12

    inv_sigma2_1 = 1.0/y1err**2
    inv_sigma2_2 = 1.0/y2err**2
    inv_sigma2_3 = 1.0/y3err**2
    if (usevs2 == 'yes'):
        inv_sigma2_4 = 1.0/y4err**2
    else:
        inv_sigma2_4 = 0.0

    lnlike_out = -0.5*(np.sum((y1-model1)**2*inv_sigma2_1)+\
                 np.sum((y2-model2)**2*inv_sigma2_2)+\
                 np.sum((y3-model3)**2*inv_sigma2_3)+\
                 np.sum((y4-model4)**2*inv_sigma2_4))

    if (SIDM_Draco_priors == 'cosmo_c'):
        #Add the conc. to the likelihood function
        #as a Gaussian in logspace:
        sigc = 0.1
        M200 = Mparsu[0]
        log_conc = np.log10(Mparsu[1])
        log_cmean = np.log10(cosmo_cfunc(M200,h))
        lnlike_out = lnlike_out * \
                     np.exp(-(log_conc-log_cmean)**2.0/(2.0*sigc**2.0))
    
    if (lnlike_out != lnlike_out):
        print sigLOS2
        print M(x1,Mpars)
        sys.exit(0)

    return lnlike_out

def lnlike_single_prop(theta, x1, x2, y1, y1err, y2, y2err, \
                       y3, y3err, y4, y4err):
    betpars = theta[0:n_betpars]
    nupars = theta[n_betpars:n_betpars+nu_components*2]
    Mpars = theta[n_betpars+nu_components*2:\
                  n_betpars+nu_components*2+n_mpars]
    Mstar = theta[n_betpars+nu_components*2+n_mpars]
    nuparsu = np.array(nupars)
    Mparsu = np.array(Mpars)
    Mparsu[0] = 10.**Mpars[0]

    sigr2, Sig, sigLOS2, sigpmr2, sigpmt2 = \
        sigp_fit_prop(x1,x2,nuparsu,Mparsu,betpars,Mstar)

    model1 = Sig
    model2 = np.sqrt(sigLOS2)/1000.
    model3 = np.sqrt(sigpmr2)/1000.
    model4 = np.sqrt(sigpmt2)/1000.

    inv_sigma2_1 = 1.0/y1err**2
    inv_sigma2_2 = 1.0/y2err**2
    inv_sigma2_3 = 1.0/y3err**2
    inv_sigma2_4 = 1.0/y4err**2

    lnlike_out = -0.5*(np.sum((y1-model1)**2*inv_sigma2_1)+\
                       np.sum((y2-model2)**2*inv_sigma2_2)+\
                       np.sum((y3-model3)**2*inv_sigma2_3)+\
                       np.sum((y4-model4)**2*inv_sigma2_4))

    if (SIDM_Draco_priors == 'cosmo_c'):
        #Add the conc. to the likelihood function
        #as a Gaussian in logspace:
        sigc = 0.1
        M200 = Mparsu[0]
        log_conc = np.log10(Mparsu[1])
        log_cmean = np.log10(cosmo_cfunc(M200,h))
        lnlike_out = lnlike_out * \
                     np.exp(-(log_conc-log_cmean)**2.0/(2.0*sigc**2.0))
        
    if (lnlike_out != lnlike_out):
        lnlike_out = -np.inf            
    
    return lnlike_out

def lnlike_split(theta, x1, x2, x3, x4, y1, y1err, \
                y2, y2err, y3, y3err, y4, y4err):
    #N.B. (x1,y1) should contain the surface density; 
    #(x2,y2) the projected dispersion;
    #and similarly for (x3,y3) and (x4,y4) but for the second component.
    betpars1 = theta[0:n_betpars]
    betpars2 = theta[n_betpars:2*n_betpars]
    nupars1 = theta[n_betpars*2:n_betpars*2+nu_components*2]
    nupars2 = theta[n_betpars*2+nu_components*2:\
                    n_betpars*2+nu_components*2*2]
    Mpars = theta[n_betpars*2+nu_components*2*2:\
                  n_betpars*2+nu_components*2*2+n_mpars]
    Mstar = theta[n_betpars*2+nu_components*2*2+n_mpars]

    nupars1u = np.array(nupars1)
    nupars2u = np.array(nupars2)
    Mparsu = np.array(Mpars)
    Mparsu[0] = 10.**Mpars[0]

    sigr2_1, Sig_1, sigLOS2_1 = \
        sigp_fit(x1,x2,nupars1u,Mparsu,betpars1,Mstar)
    sigr2_2, Sig_2, sigLOS2_2 = \
        sigp_fit(x3,x4,nupars2u,Mparsu,betpars2,Mstar)

    model1 = Sig_1
    model2 = np.sqrt(sigLOS2_1)/1000.
    model3 = Sig_2
    model4 = np.sqrt(sigLOS2_2)/1000.

    inv_sigma2_1 = 1.0/y1err**2
    inv_sigma2_2 = 1.0/y2err**2
    inv_sigma2_3 = 1.0/y3err**2
    inv_sigma2_4 = 1.0/y4err**2

    return -0.5*(np.sum((y1-model1)**2*inv_sigma2_1)+\
                 np.sum((y2-model2)**2*inv_sigma2_2)+\
                 np.sum((y3-model3)**2*inv_sigma2_3)+\
                 np.sum((y4-model4)**2*inv_sigma2_4))

def lnlike_split_vs(theta, x1, x2, x3, x4, y1, y1err, \
                    y2, y2err, y3, y3err, y4, y4err, \
                    y5, y5err, y6, y6err, y7, y7err, \
                    y8, y8err):
    #N.B. (x1,y1) should contain the surface density;
    #(x2,y2) the projected dispersion;
    #and similarly for (x3,y3) and (x4,y4) but for the second component.
    #y5 - y8 contain the "virial parameters" and their errors:
    betpars1 = theta[0:n_betpars]
    betpars2 = theta[n_betpars:2*n_betpars]
    nupars1 = theta[n_betpars*2:n_betpars*2+nu_components*2]
    nupars2 = theta[n_betpars*2+nu_components*2:\
                    n_betpars*2+nu_components*2*2]

    Mpars = theta[n_betpars*2+nu_components*2*2:\
                  n_betpars*2+nu_components*2*2+n_mpars]
    Mstar = theta[n_betpars*2+nu_components*2*2+n_mpars]

    nupars1u = np.array(nupars1)
    nupars2u = np.array(nupars2)
    Mparsu = np.array(Mpars)
    Mparsu[0] = 10.**Mpars[0]

    sigr2_1, Sig_1, sigLOS2_1, vs1_1, vs2_1 = \
        sigp_fit_vs(x1,x2,nupars1u,Mparsu,betpars1,Mstar)
    sigr2_2, Sig_2, sigLOS2_2, vs1_2, vs2_2 = \
        sigp_fit_vs(x3,x4,nupars2u,Mparsu,betpars2,Mstar)

    model1 = Sig_1
    model2 = np.sqrt(sigLOS2_1)/1000.
    model3 = Sig_2
    model4 = np.sqrt(sigLOS2_2)/1000.
    model5 = vs1_1/1.0e12
    model6 = vs2_1/1.0e12
    model7 = vs1_2/1.0e12
    model8 = vs2_2/1.0e12

    inv_sigma2_1 = 1.0/y1err**2
    inv_sigma2_2 = 1.0/y2err**2
    inv_sigma2_3 = 1.0/y3err**2
    inv_sigma2_4 = 1.0/y4err**2
    inv_sigma2_5 = 1.0/y5err**2
    inv_sigma2_6 = 1.0/y6err**2
    inv_sigma2_7 = 1.0/y7err**2
    inv_sigma2_8 = 1.0/y8err**2

    return -0.5*(np.sum((y1-model1)**2*inv_sigma2_1)+\
                 np.sum((y2-model2)**2*inv_sigma2_2)+\
                 np.sum((y3-model3)**2*inv_sigma2_3)+\
                 np.sum((y4-model4)**2*inv_sigma2_4)+\
                 np.sum((y5-model5)**2*inv_sigma2_5)+\
                 np.sum((y6-model6)**2*inv_sigma2_6)+\
                 np.sum((y7-model7)**2*inv_sigma2_7)+\
                 np.sum((y8-model8)**2*inv_sigma2_8))

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
        nineninehigh, nineninelow

def binthedata(R,vz,ms,Nbin,verr):
    #Information on error calculation for binned data:
    #http://esoads.eso.org/abs/1993ASPC...50..357P

    #Nbin is the number of particles / bin:
    index = np.argsort(R)
    rbin = np.zeros(len(R))
    norm = np.zeros(len(R))
    sigpmean2 = np.zeros(len(R))
    vmean = np.zeros(len(R))
    surfden = np.zeros(len(R))
    cnt = 0
    jsum = 0
    for i in range(len(R)):
        if (jsum < Nbin):
            norm[cnt] = norm[cnt] + ms[index[i]]
            sigpmean2[cnt] = sigpmean2[cnt] + \
                vz[index[i]]**2.*ms[index[i]]
            vmean[cnt] = vmean[cnt] + \
                vz[index[i]]*ms[index[i]]
            rbin[cnt] = R[index[i]]
            jsum = jsum + ms[index[i]]
        if (jsum >= Nbin):
            jsum = 0.0
            cnt = cnt + 1
    
    rbin = rbin[:cnt]
    norm = norm[:cnt]
    sigpmean2 = sigpmean2[:cnt]
    vmean = vmean[:cnt]
    surfden = surfden[:cnt]

    for i in range(len(rbin)):
        if (i == 0):
            surfden[i] = norm[i] / (np.pi * rbin[i]**2.0)
        else:
            surfden[i] = norm[i] / (2.0*np.pi*rbin[i]*(rbin[i]-rbin[i-1]))
    surfdenerr = surfden / np.sqrt(Nbin)
    sigpmean2 = sigpmean2 / norm
    vmean = vmean / norm

    #And deal with errors on sigpmean (see above ref.):
    sigpmean2 = sigpmean2 - verr**2.0
    sigperr2 = (sigpmean2 + verr**2.0)**2.0 / (2.0*Nbin*sigpmean2)
    sigpmean = np.sqrt(sigpmean2)
    sigperr = np.sqrt(sigperr2)

    #Calculate the projected half light radius: 
    ranal = np.linspace(0,10,5000)
    surfden_ranal = np.interp(ranal,rbin,surfden,left=0,right=0)
    Menc_tot = 2.0*np.pi*integrator(surfden_ranal*ranal,ranal)
    Menc_half = 0.0
    i = 3
    while (Menc_half < Menc_tot/2.0):
        Menc_half = 2.0*np.pi*\
            integrator(surfden_ranal[:i]*ranal[:i],ranal[:i])
        i = i + 1
    Rhalf = ranal[i-1]

    #And set up the light profile in case Mstar > 0:
    pfits = tracerfit(p0in,p0in_min,p0in_max,\
                      rbin,surfden,surfdenerr)
    Mstar_rad = rbin
    norm = np.max(threeplummass(\
            np.linspace(0,Mstar_rlim,100),pfits[0],pfits[1],\
                pfits[2],\
                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/norm,pfits[1]/norm,\
                                   pfits[2]/norm,\
                                   pfits[3],pfits[4],pfits[5])

    return rbin, surfden, surfdenerr,\
        sigpmean, sigperr, Rhalf, vmean, Mstar_rad, Mstar_prof

def binthedata_prop(x,y,z,vx,vy,vz,ms,Nbin,verr):
    #Set up coordinate system which is plane polar on the sky
    #a la Strigari 2007 and van der Marel & Anderson 2010: 
    R = np.sqrt(x**2. + y**2.)
    vR = (vx*x+vy*y)/R
    vphi = (vx*y-vy*x)/R

    index = np.argsort(R)

    rbin = np.zeros(len(R))
    norm = np.zeros(len(R))
    sigpmean2 = np.zeros(len(R))
    sigpmr2 = np.zeros(len(R))
    sigpmt2 = np.zeros(len(R))
    surfden = np.zeros(len(R))
    cnt = 0
    j = 0

    for i in range(len(R)):
        if (j < Nbin):
            norm[cnt] = norm[cnt] + ms[index[i]]
            sigpmean2[cnt] = sigpmean2[cnt] + \
                vz[index[i]]**2.*ms[index[i]]
            sigpmr2[cnt] = sigpmr2[cnt] + \
                vR[index[i]]**2.*ms[index[i]]
            sigpmt2[cnt] = sigpmt2[cnt] + \
                vphi[index[i]]**2.*ms[index[i]]
            j = j + 1
            rbin[cnt] = R[index[i]]
        else:
            j = 0
            cnt = cnt + 1

    rbin = rbin[:cnt]
    norm = norm[:cnt]
    sigpmean2 = sigpmean2[:cnt]
    sigpmr2 = sigpmr2[:cnt]
    sigpmt2 = sigpmt2[:cnt]
    surfden = surfden[:cnt]

    for i in range(len(rbin)):
        if (i == 0):
            surfden[i] = norm[i] / (np.pi * rbin[i]**2.0)
        else:
            surfden[i] = norm[i] / (2.0*np.pi*rbin[i]*(rbin[i]-rbin[i-1]))
    surfdenerr = surfden / np.sqrt(Nbin)
    sigpmean2 = sigpmean2 / norm
    sigpmr2 = sigpmr2 / norm
    sigpmt2 = sigpmt2 / norm

    sigpmean2 = sigpmean2 - verr**2.0
    sigperr2 = (sigpmean2 + verr**2.0)**2.0 / (2.0*Nbin*sigpmean2)
    sigpmean = np.sqrt(sigpmean2)
    sigperr = np.sqrt(sigperr2)
    sigpmt2 = sigpmt2 - verr**2.0
    sigpmterr2 = (sigpmt2 + verr**2.0)**2.0 / (2.0*Nbin*sigpmt2)
    sigpmt = np.sqrt(sigpmt2)
    sigpmterr = np.sqrt(sigpmterr2)
    sigpmr2 = sigpmr2 - verr**2.0
    sigpmrerr2 = (sigpmr2 + verr**2.0)**2.0 / (2.0*Nbin*sigpmr2)
    sigpmr = np.sqrt(sigpmr2)
    sigpmrerr = np.sqrt(sigpmrerr2)

    ranal = np.linspace(0,10,5000)
    surfden_ranal = np.interp(ranal,rbin,surfden,left=0,right=0)
    Menc_tot = 2.0*np.pi*integrator(surfden_ranal*ranal,ranal)
    Menc_half = 0.0
    i = 3
    while (Menc_half < Menc_tot/2.0):
        Menc_half = 2.0*np.pi*\
            integrator(surfden_ranal[:i]*ranal[:i],ranal[:i])
        i = i + 1
    Rhalf = ranal[i-1]

    return rbin, surfden, surfdenerr,\
        sigpmean, sigperr, \
        sigpmt, sigpmterr, sigpmr, sigpmrerr, Rhalf

def alpbetgamden(r,rho0,r0,alp,bet,gam):
    return rho0*(r/r0)**(-gam)*(1.0+(r/r0)**alp)**((gam-bet)/alp)

def alpbetgamdlnrhodlnr(r,rho0,r0,alp,bet,gam):
    return -gam + (gam-bet)*(r/r0)**alp*(1.0+(r/r0)**alp)**(-1.0)

def alpbetgammass(r,rho0,r0,alp,bet,gam):
    den = rho0*(r/r0)**(-gam)*(1.0+(r/r0)**alp)**((gam-bet)/alp)
    mass = np.zeros(len(r))
    for i in range(1,len(r)):
        if (i == 1):
            mass[i] = mass[i-1] + 4.0*np.pi*r[i]**2.*den[i]*(r[i]-r[i-1])/2.0
        else:
            mass[i] = mass[i-1] + 4.0*np.pi*r[i]**2.*den[i]*(r[i]-r[i-1])
    return mass

def threeplumsurf(r,M1,M2,M3,a1,a2,a3):
    return multiplumsurf(r,[M1,M2,M3,\
                                a1,a2,a3])
def threeplummass(r,M1,M2,M3,a1,a2,a3):
    return multiplummass(r,[M1,M2,M3,\
                                a1,a2,a3])
def residual(params, x, data, eps_data):
    M1 = params['M1'].value
    M2 = params['M2'].value
    M3 = params['M3'].value
    a1 = params['a1'].value
    a2 = params['a2'].value
    a3 = params['a3'].value

    model = threeplumsurf(x,M1,M2,M3,a1,a2,a3)

    return (data-model)/eps_data

def residual_powline(params, x, data, eps_data):
    A = params['A'].value
    B = params['B'].value
    C = params['C'].value

    model = A*x**B+C

    return (data-model)/eps_data

def tracerfit(p0in,p0in_min,p0in_max,\
              rbin_phot,surfden,surfdenerr):
    params = lm.Parameters()
    params.add('M1', value=p0in[0],min=p0in_min[0],max=p0in_max[0])
    params.add('M2', value=p0in[1],min=p0in_min[1],max=p0in_max[1])
    params.add('M3', value=p0in[2],min=p0in_min[2],max=p0in_max[2])
    params.add('a1', value=p0in[3],min=p0in_min[3],max=p0in_max[3])
    params.add('a2', value=p0in[4],min=p0in_min[4],max=p0in_max[4])
    params.add('a3', value=p0in[5],min=p0in_min[5],max=p0in_max[5])

    out = lm.minimize(residual, params,\
                      args=(rbin_phot, surfden, surfdenerr))
    pfits = np.zeros(6)
    pfits[0] = out.params['M1'].value
    pfits[1] = out.params['M2'].value
    pfits[2] = out.params['M3'].value
    pfits[3] = out.params['a1'].value
    pfits[4] = out.params['a2'].value
    pfits[5] = out.params['a3'].value

    return pfits

def v4func(x,A,B,C):
    return A*x**B + C

def v4residual(params, x, data, eps_data):
    A = params['A'].value
    B = params['B'].value
    C = params['C'].value

    model = v4func(x,A,B,C)

    return (data-model)/eps_data

def v4fit(p0in,p0in_min,p0in_max,\
          rbin,v4,v4err):
    params = lm.Parameters()
    params.add('A', value=p0in[0],min=p0in_min[0],max=p0in_max[0])
    params.add('B', value=p0in[1],min=p0in_min[1],max=p0in_max[1])
    params.add('C', value=p0in[2],min=p0in_min[2],max=p0in_max[2])

    out = lm.minimize(v4residual, params, \
                      args=(rbin, v4, v4err))

    pfits = np.zeros(3)
    pfits[0] = out.params['A'].value
    pfits[1] = out.params['B'].value
    pfits[2] = out.params['C'].value

    return pfits

def tvl4func(rint,rbin_tmp,vlos4med,A,B,C,rout,gamout,zero):
    tvl4 = np.interp(rint,rbin_tmp,vlos4med,left=0,right=0)
    if (zero == 'no'):
        for i in range(len(tvl4)):
            if (rint[i] > np.max(rbin_tmp)):
                tvl4[i] = A*(rint[i]/np.max(rbin_tmp))**B+C
            if (rint[i] > rout):
                tvl4[i] = (A*(rout/np.max(rbin_tmp))**B+C)*\
                    (rint[i]/rout)**(-gamout)
    return tvl4

def fit_powline(rbin_tmp,vlos4med,vlos4err):
    max_r = np.max(rbin_tmp)
    j=0L
    while(rbin_tmp[j] < Rhalf):
        j=j+1
    ruse = rbin_tmp[j:]/max_r
    vuse = vlos4med[j:]
    verr = vlos4err[j:]
        
    params_powline = lm.Parameters()
    params_powline.add('A', value=np.max(vlos4med),\
               min=np.min(vlos4med)/100.0,\
               max=np.max(vlos4med)*100.0)
    params_powline.add('B', value=0.0,\
               min=-2.0,\
               max=2.0)
    params_powline.add('C', value=0.0,vary=False)

    out = lm.minimize(residual_powline, params_powline, \
                      args=(ruse, vuse, verr))
    pfits_powline = np.zeros(3)
    pfits_powline[0] = out.params['A'].value
    pfits_powline[1] = out.params['B'].value
    pfits_powline[2] = out.params['C'].value

    return pfits_powline

def walker_surf_sort(data_file,Rhalf_lum):
    f = open(data_file,'r')
    data_phot = np.genfromtxt(f)
    f.close()
    rbin_phot_t = data_phot[:,0]*dgal/arcmin
    surfden_t = data_phot[:,1]
    norm = integrator(surfden_t*2.0*np.pi*rbin_phot_t,rbin_phot_t)
    surfden_t = surfden_t / norm
    nstars_per_bin = np.zeros(len(data_phot[:,0]))

    for ii in range(len(data_phot[:,0])):
        if (ii < len(data_phot[:,0])-1):
            dR = data_phot[ii+1,0]-data_phot[ii,0]
        nstars_per_bin[ii] = data_phot[ii,1]*2.0*np.pi*data_phot[ii,0]*dR
    surfdenerr_t = surfden_t / np.sqrt(nstars_per_bin)

    rbin_phot = rbin_phot_t[rbin_phot_t < maxdatrad]
    surfden = surfden_t[rbin_phot_t < maxdatrad]
    surfdenerr = surfdenerr_t[rbin_phot_t < maxdatrad]
    rbin_photfit = rbin_phot_t[rbin_phot_t < maxdatfitrad]
    surfdenfit = surfden_t[rbin_phot_t < maxdatfitrad]
    surfdenerrfit = surfdenerr_t[rbin_phot_t < maxdatfitrad]

    pfits = tracerfit(p0in,p0in_min,p0in_max,\
                          rbin_photfit,surfdenfit,surfdenerrfit)
    Mstar_rad = rbin_phot
    norm = np.max(threeplummass(\
            np.linspace(0,Mstar_rlim,100),pfits[0],pfits[1],\
                pfits[2],\
                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/norm,pfits[1]/norm,\
                                   pfits[2]/norm,\
                                   pfits[3],pfits[4],pfits[5])
    Mstar_surf = threeplumsurf(ranal,pfits[0]/norm,pfits[1]/norm,\
                                   pfits[2]/norm,\
                                   pfits[3],pfits[4],pfits[5])
    if (Rhalf_lum < 0):
        Mcum_surf = 0.0
        i = 1
        while (Mcum_surf < (pfits[0]+pfits[1]+pfits[2])/2.0/norm):
            Mcum_surf = \
                Mcum_surf + \
                2.0*np.pi*ranal[i]*Mstar_surf[i]*(ranal[i]-ranal[i-1])
            i = i + 1
        Rhalf = ranal[i-1]
        print 'Rhalf calculated: ', Rhalf
    else:
        Rhalf = Rhalf_lum

    return rbin_phot, surfden, surfdenerr, \
        rbin_photfit, surfdenfit, surfdenerrfit, \
        Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits

def tucana_surf_sort(data_file,Rhalf_lum):
    f = open(data_file,'r')
    data_phot = np.genfromtxt(f)
    f.close()
    
    R = data_phot[:,5]
    ms = data_phot[:,8]
    
    cnt = 0
    jsum = 0.0
    norm = np.zeros(len(R))
    rbin_phot_t = np.zeros(len(R))
    surfden_t = np.zeros(len(R))
    index = np.argsort(R)
    for i in range(len(R)):
        if (jsum < Nbin):
            surfden_t[cnt] = surfden_t[cnt] + ms[index[i]]
            jsum = jsum + ms[index[i]]
            rbin_phot_t[cnt] = R[index[i]]
        if (jsum >= Nbin):
            norm[cnt] = jsum
            if (cnt == 0):
                area = np.pi*rbin_phot_t[cnt]**2.0
            else:
                area = np.pi*(rbin_phot_t[cnt]**2.0-\
                              rbin_phot_t[cnt-1]**2.0)
            surfden_t[cnt] = surfden_t[cnt]/area
            jsum = 0.0
            cnt = cnt + 1
            
    surfdenerr_t = surfden_t / np.sqrt(norm)
    rbin_phot_t = rbin_phot_t[:cnt]
    surfden_t = surfden_t[:cnt]
    surfdenerr_t = surfdenerr_t[:cnt]
    rbin_phot = rbin_phot_t[rbin_phot_t < maxdatrad]
    surfden = surfden_t[rbin_phot_t < maxdatrad]
    surfdenerr = surfdenerr_t[rbin_phot_t < maxdatrad]
    rbin_photfit = rbin_phot_t[rbin_phot_t < maxdatfitrad]
    surfdenfit = surfden_t[rbin_phot_t < maxdatfitrad]
    surfdenerrfit = surfdenerr_t[rbin_phot_t < maxdatfitrad]
    
    pfits = tracerfit(p0in,p0in_min,p0in_max,\
                      rbin_photfit,surfdenfit,surfdenerrfit)

    Mstar_rad = rbin_phot
    norm = np.max(threeplummass(\
                                np.linspace(0,Mstar_rlim,100),pfits[0],pfits[1],\
                                pfits[2],\
                                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])
    Mstar_surf = threeplumsurf(ranal,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])

    if (Rhalf_lum < 0):
        Mcum_surf = 0.0
        i = 1
        while (Mcum_surf < (pfits[0]+pfits[1]+pfits[2])/2.0/norm):
            Mcum_surf = \
                        Mcum_surf + \
                        2.0*np.pi*ranal[i]*Mstar_surf[i]*(ranal[i]-ranal[i-1])
            i = i + 1
        Rhalf = ranal[i-1]
        print 'Rhalf calculated: ', Rhalf
    else:
        Rhalf = Rhalf_lum
        
    return rbin_phot, surfden, surfdenerr, \
        rbin_photfit, surfdenfit, surfdenerrfit, \
        Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits

def walker_surf_sort_cra(data_file,Rhalf_lum):
    f = open(data_file,'r')
    data_phot = np.genfromtxt(f)
    f.close()

    R = data_phot[:,4]*dgal/arcmin
    ms = data_phot[:,31]

    cnt = 0
    jsum = 0.0
    norm = np.zeros(len(R))
    rbin_phot_t = np.zeros(len(R))
    surfden_t = np.zeros(len(R))
    index = np.argsort(R)
    for i in range(len(R)):
        if (jsum < Nbin):
            surfden_t[cnt] = surfden_t[cnt] + ms[index[i]]
            jsum = jsum + ms[index[i]]
            rbin_phot_t[cnt] = R[index[i]]
        if (jsum >= Nbin):
            norm[cnt] = jsum
            if (cnt == 0):
                area = np.pi*rbin_phot_t[cnt]**2.0
            else:
                area = np.pi*(rbin_phot_t[cnt]**2.0-\
                              rbin_phot_t[cnt-1]**2.0)
            surfden_t[cnt] = surfden_t[cnt]/area
            jsum = 0.0
            cnt = cnt + 1
    surfdenerr_t = surfden_t / np.sqrt(norm)
    rbin_phot_t = rbin_phot_t[:cnt]
    surfden_t = surfden_t[:cnt]
    surfdenerr_t = surfdenerr_t[:cnt]

    rbin_phot = rbin_phot_t[rbin_phot_t < maxdatrad]
    surfden = surfden_t[rbin_phot_t < maxdatrad]
    surfdenerr = surfdenerr_t[rbin_phot_t < maxdatrad]
    rbin_photfit = rbin_phot_t[rbin_phot_t < maxdatfitrad]
    surfdenfit = surfden_t[rbin_phot_t < maxdatfitrad]
    surfdenerrfit = surfdenerr_t[rbin_phot_t < maxdatfitrad]

    pfits = tracerfit(p0in,p0in_min,p0in_max,\
                      rbin_photfit,surfdenfit,surfdenerrfit)

    Mstar_rad = rbin_phot
    norm = np.max(threeplummass(\
            np.linspace(0,Mstar_rlim,100),pfits[0],pfits[1],\
                pfits[2],\
                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/norm,pfits[1]/norm,\
                                   pfits[2]/norm,\
                                   pfits[3],pfits[4],pfits[5])
    Mstar_surf = threeplumsurf(ranal,pfits[0]/norm,pfits[1]/norm,\
                                   pfits[2]/norm,\
                                   pfits[3],pfits[4],pfits[5])

    if (Rhalf_lum < 0):
        Mcum_surf = 0.0
        i = 1
        while (Mcum_surf < (pfits[0]+pfits[1]+pfits[2])/2.0/norm):
            Mcum_surf = \
                Mcum_surf + \
                2.0*np.pi*ranal[i]*Mstar_surf[i]*(ranal[i]-ranal[i-1])
            i = i + 1
        Rhalf = ranal[i-1]
        print 'Rhalf calculated: ', Rhalf
    else:
        Rhalf = Rhalf_lum

    return rbin_phot, surfden, surfdenerr, \
        rbin_photfit, surfdenfit, surfdenerrfit, \
        Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits

def genina_surf_sort(R,ms):
    cnt = 0
    jsum = 0.0
    norm = np.zeros(len(R))
    rbin_phot_t = np.zeros(len(R))
    surfden_t = np.zeros(len(R))
    index = np.argsort(R)
    
    for i in range(len(R)):
        if (jsum < Nbin):
            surfden_t[cnt] = surfden_t[cnt] + ms[index[i]]
            jsum = jsum + ms[index[i]]
            rbin_phot_t[cnt] = R[index[i]]
        if (jsum >= Nbin):
            norm[cnt] = jsum
            if (cnt == 0):
                area = np.pi*rbin_phot_t[cnt]**2.0
            else:
                area = np.pi*(rbin_phot_t[cnt]**2.0-\
                              rbin_phot_t[cnt-1]**2.0)
            surfden_t[cnt] = surfden_t[cnt]/area
            jsum = 0.0
            cnt = cnt + 1

    surfdenerr_t = surfden_t / np.sqrt(norm)
    rbin_phot_t = rbin_phot_t[:cnt]
    surfden_t = surfden_t[:cnt]
    surfdenerr_t = surfdenerr_t[:cnt]
    rbin_phot = rbin_phot_t[rbin_phot_t < maxdatrad]
    surfden = surfden_t[rbin_phot_t < maxdatrad]
    surfdenerr = surfdenerr_t[rbin_phot_t < maxdatrad]
    rbin_photfit = rbin_phot_t[rbin_phot_t < maxdatfitrad]
    surfdenfit = surfden_t[rbin_phot_t < maxdatfitrad]
    surfdenerrfit = surfdenerr_t[rbin_phot_t < maxdatfitrad]
    
    pfits = tracerfit(p0in,p0in_min,p0in_max,\
                      rbin_photfit,surfdenfit,surfdenerrfit)

    Mstar_rad = rbin_phot
    norm = np.max(threeplummass(\
                                np.linspace(0,Mstar_rlim,100),\
                                pfits[0],pfits[1],\
                                pfits[2],\
                                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])
    Mstar_surf = threeplumsurf(ranal,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])

    Mcum_surf = 0.0
    i = 1
    while (Mcum_surf < (pfits[0]+pfits[1]+pfits[2])/2.0/norm):
        Mcum_surf = \
                    Mcum_surf + \
                    2.0*np.pi*ranal[i]*Mstar_surf[i]*(ranal[i]-ranal[i-1])
        i = i + 1
    Rhalf = ranal[i-1]
    print 'Rhalf calculated: ', Rhalf

    return rbin_phot, surfden, surfdenerr, \
        rbin_photfit, surfdenfit, surfdenerrfit, \
        Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits

def calc_virial_moments(Rhalf,nmonte,R,vz,vzerr,ms,pfits):
    #Improve estimator using fitted surfden:
    rint = np.logspace(np.log10(Rhalf/100.0),\
                       np.log10(Rhalf*1000.0),10000)
    index = np.argsort(R)

    #First calculate and subtract the mean vz:
    use_meanvz_by_bin = 'no'
    if (use_meanvz_by_bin == 'yes'):
        cnt = 0
        jsum = 0.0
        norm = np.zeros(len(R))
        vzmean = np.zeros(len(R))
        for i in range(len(R)):
            if (jsum < Nbin):
                vzmean[cnt] = vzmean[cnt] + \
                    vz[index[i]]*ms[index[i]]
                jsum = jsum + ms[index[i]]
            if (jsum >= Nbin):
                norm[cnt] = jsum
                jsum = 0.0
                cnt = cnt + 1
        vzmean = vzmean / norm
        vzmeanerr = vzmean / np.sqrt(Nbin)
    else:
        vzmean = np.full(len(vz),np.sum(vz*ms)/np.sum(ms))
        vzmeanerr = np.zeros(len(vzmean))
        
    print '|Vzmean| min/max', np.min(np.abs(vzmean[np.abs(vzmean) < 1e30])),\
        np.max(np.abs(vzmean[np.abs(vzmean) < 1e30]))
        
    #And now the 2nd and 4th moments:
    cnt = 0
    jsum = 0.0
    norm = np.zeros(len(R))
    vlos4med = np.zeros(len(R))
    vlos2med = np.zeros(len(R))
    rbin_tmp = np.zeros(len(R))
    for i in range(len(R)):
        if (jsum < Nbin):
            vlos4med[cnt] = vlos4med[cnt] + \
                (vz[index[i]]-vzmean[cnt])**4.*ms[index[i]]
            vlos2med[cnt] = vlos2med[cnt] + \
                (vz[index[i]]-vzmean[cnt])**2.*ms[index[i]]
            rbin_tmp[cnt] = R[index[i]]
            jsum = jsum + ms[index[i]]
        if (jsum >= Nbin):
            norm[cnt] = jsum
            jsum = 0.0
            cnt = cnt + 1
    vlos4med = vlos4med[:cnt]
    vlos2med = vlos2med[:cnt]
    norm = norm[:cnt]
    vlos4med = vlos4med / norm
    vlos2med = vlos2med / norm
    rbin_tmp = rbin_tmp[:cnt]

    #And Monte-Carlo the errors:
    vlos4 = np.zeros((nmonte,len(R)))
    vlos2 = np.zeros((nmonte,len(R)))
    vlos4_pureerr = np.zeros((nmonte,len(R)))
    vlos2_pureerr = np.zeros((nmonte,len(R)))
    norm = np.zeros(len(R))
    for k in range(nmonte):
        cnt = 0
        jsum = 0.0
        for i in range(len(R)):
            vz_err = (vz[index[i]]-vzmean[cnt])+\
                np.random.normal(0.0,vzerr[index[i]])
            vz_pure_err = np.random.normal(0.0,vzerr[index[i]])
            if (jsum < Nbin):
                vlos4[k,cnt] = vlos4[k,cnt] + \
                    vz_err**4.*ms[index[i]]
                vlos2[k,cnt] = vlos2[k,cnt] + \
                    vz_err**2.*ms[index[i]]
                vlos4_pureerr[k,cnt] = vlos4_pureerr[k,cnt] + \
                    vz_pure_err**4.*ms[index[i]]
                vlos2_pureerr[k,cnt] = vlos2_pureerr[k,cnt] + \
                    vz_pure_err**2.*ms[index[i]]
                jsum = jsum + ms[index[i]]
            if (jsum >= Nbin):
                norm[cnt] = jsum
                jsum = 0.0
                cnt = cnt + 1
    vlos4tmp = np.zeros((nmonte,cnt))
    vlos4tmp = vlos4[:,:cnt]
    vlos2tmp = np.zeros((nmonte,cnt))
    vlos2tmp = vlos2[:,:cnt]
    vlos4_pe_tmp = np.zeros((nmonte,cnt))
    vlos4_pe_tmp = vlos4_pureerr[:,:cnt]
    vlos2_pe_tmp = np.zeros((nmonte,cnt))
    vlos2_pe_tmp = vlos2_pureerr[:,:cnt]
    norm = norm[:cnt]

    vlos4 = vlos4tmp / norm
    vlos2 = vlos2tmp / norm
    vlos4_pe = vlos4_pe_tmp / norm
    vlos2_pe = vlos2_pe_tmp / norm

    #And now estimate the full measurement error:
    vlos4err_meas = np.zeros(cnt)
    vlos2err_meas = np.zeros(cnt)
    vlos4_pe_meas = np.zeros(cnt)
    vlos2_pe_meas = np.zeros(cnt)
    for k in range(cnt):
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = calcmedquartnine(vlos4[:,k])
        vlos4err_meas[k] = (sixhigh-sixlow)/2.0
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = calcmedquartnine(vlos2[:,k])
        vlos2err_meas[k] = (sixhigh-sixlow)/2.0
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = calcmedquartnine(vlos4_pe[:,k])
        vlos4_pe_meas[k] = (sixhigh-sixlow)/2.0
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = calcmedquartnine(vlos2_pe[:,k])
        vlos2_pe_meas[k] = (sixhigh-sixlow)/2.0

    #Combine with the Poisson error: 
    vlos4err = np.sqrt(vlos4err_meas**2.0 + vlos4med**2.0/Nbin)
    vlos2err = np.sqrt(vlos2err_meas**2.0 + vlos2med**2.0/Nbin)
    vlos4med = vlos4med - vlos4_pe_meas
    vlos2med = vlos2med - vlos2_pe_meas

    #Demand positive:
    vzmean = vzmean[:cnt]
    vzmeanerr = vzmeanerr[:cnt]

    if (virialshape == 'yes'):
        vlos2med = vlos2med[vlos4med > 0]
        vlos2err = vlos2err[vlos4med > 0]
        vlos4err = vlos4err[vlos4med > 0]
        rbin_tmp = rbin_tmp[vlos4med > 0]
        vzmean = vzmean[vlos4med > 0]
        vzmeanerr = vzmeanerr[vlos4med > 0]
        vlos4med = vlos4med[vlos4med > 0]
    else:
        vlos2err = vlos2err[vlos2med > 0]
        vlos4err = vlos4err[vlos2med > 0]
        rbin_tmp = rbin_tmp[vlos2med > 0]
        vzmean = vzmean[vlos2med > 0]
        vzmeanerr = vzmeanerr[vlos2med > 0]
        vlos4med = vlos4med[vlos2med > 0]
        vlos2med = vlos2med[vlos2med > 0]

    #If zero not set, fit a smooth function for 
    #interpolation, else set "beyond-data" = 0:
    p0in = np.array([1.0,0.0,1000.0])
    p0in_min = np.zeros(3)
    p0in_max = np.array([1e5,1.0,1e5])
    if (zero == 'no'):
        print 'Using outer interpolation scheme for v4los'
        if (len(rbin_tmp) < 6):
            fitmin = 0
        else:
            fitmin = len(rbin_tmp)/2
        pfitmed = v4fit(p0in,p0in_min,p0in_max,\
                        rbin_tmp[fitmin:],vlos4med[fitmin:],\
                        vlos4err[fitmin:])
    else:
        pfitmed = p0in
        print 'Assigning zeros beyond data for v4los'

    #Cut back to maxdatrad:
    rbin_tmp_full = rbin_tmp
    vlos4err_full = vlos4err
    vlos2err_full = vlos2err
    vlos4med_full = vlos4med
    vlos2med_full = vlos2med
    vzmean_full = vzmean
    vzmeanerr_full = vzmeanerr
    vlos4err = vlos4err[rbin_tmp < maxdatrad]
    vlos2err = vlos2err[rbin_tmp < maxdatrad]
    vlos4med = vlos4med[rbin_tmp < maxdatrad]
    vlos2med = vlos2med[rbin_tmp < maxdatrad]
    vzmean = vzmean[rbin_tmp < maxdatrad]
    vzmeanerr = vzmeanerr[rbin_tmp < maxdatrad]
    rbin_tmp = rbin_tmp[rbin_tmp < maxdatrad]
    
    #Fit straight line in log-space to vlos4med:
    if (zero == 'no'):
        pfits_powline = fit_powline(rbin_tmp,vlos4med,vlos4err)
        router_min = np.max(rbin_tmp)
        router_max = 2.0*router_min

        if (router_max < router_min):
            router_max = router_min
        gamout_min = 1.0
        gamout_max = 3.0
        router = (router_min+router_max)/2.0
        gamout = (gamout_min+gamout_max)/2.0        
    else:
        pfits_powline = np.array([1.0,1.0,1.0])
        router_min = 1.0
        router_max = 1.0
        gamout_min = 1.0
        gamout_max = 1.0
        router = 1.0
        gamout = 1.0
    
    if (testmode == 'yes'):
        plt.figure()
        plt.loglog()
        plt.errorbar(rbin_tmp_full,vlos4med_full,vlos4err_full,color='k')
        plt.errorbar(rbin_tmp,vlos4med,vlos4err,color='b')
        tvl4 = tvl4func(rint,rbin_tmp,vlos4med,\
                        pfits_powline[0],pfits_powline[1],\
                        pfits_powline[2],router,gamout,zero)
        plt.plot(rint,tvl4,'r')
        plt.show()

        plt.figure()
        plt.loglog()
        plt.errorbar(rbin_tmp_full,np.sqrt(vlos2med_full),\
                     vlos2err_full/2.0/np.sqrt(vlos2med_full),color='k')
        plt.errorbar(rbin_tmp,np.sqrt(vlos2med),\
                     vlos2err/2.0/np.sqrt(vlos2med),color='b')
        plt.show()

    #Plot the 1st, 2nd and 4th moments:
    if (len(rbin_tmp_full) > 1):
        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)

        plt.loglog()
        plt.errorbar(rbin_tmp_full,vlos4med_full,vlos4err_full,color='k')
        plt.errorbar(rbin_tmp,vlos4med,vlos4err,color='b')
        tvl4 = tvl4func(rint,rbin_tmp,vlos4med,\
                        pfits_powline[0],pfits_powline[1],\
                        pfits_powline[2],router,gamout,zero)
        plt.plot(rint,tvl4,'r')
        plt.plot(rint,pfits_powline[0]*(rint/np.max(rbin_tmp))**\
                 pfits_powline[1]+\
                 pfits_powline[2])
        plt.xlim([Rhalf/10.0,Rhalf*100.0])
        plt.ylim([1.0,y_sigLOSmax**4.0*100.0])
        plt.xlabel(r'$R\,[{\rm kpc}]$',\
                       fontsize=myfontsize)
        plt.ylabel(r'$\langle v_{\rm los}^4\rangle\,({\rm km}^4\,{\rm s}^{-4})$',\
                       fontsize=myfontsize)
        plt.savefig(outdir+'output_vlos4.pdf',bbox_inches='tight')

        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)
            
        plt.loglog()
        plt.errorbar(rbin_tmp_full,np.sqrt(vlos2med_full),\
                         vlos2err_full/2.0/np.sqrt(vlos2med_full),color='k')
        plt.errorbar(rbin_tmp,np.sqrt(vlos2med),\
                         vlos2err/2.0/np.sqrt(vlos2med),color='b')
        plt.xlim([Rhalf/10.0,Rhalf*100.0])
        plt.ylim([1.0,y_sigLOSmax])

        plt.xlabel(r'$R\,[{\rm kpc}]$',\
                       fontsize=myfontsize)
        plt.ylabel(r'$\sigma_{\rm LOS}\,({\rm km}\,{\rm s}^{-1})$',\
                       fontsize=myfontsize)
        plt.savefig(outdir+'output_vlos2.pdf',bbox_inches='tight')

    #And calculate vs1 and vs2:
    tvl4 = tvl4func(rint,rbin_tmp,vlos4med,\
                    pfits_powline[0],pfits_powline[1],pfits_powline[2],router,gamout,zero)
    test_surfden = threeplumsurf(rint,pfits[0],pfits[1],pfits[2],\
                                     pfits[3],pfits[4],pfits[5])
    vs1imp = integrator(tvl4*test_surfden*rint,rint)
    vs2imp = integrator(tvl4*test_surfden*rint**3.0,rint)

    #Monte-Carlo to calculate the ~1sigma errors:
    vs1_samp = np.zeros(nmonte)
    vs2_samp = np.zeros(nmonte)
    for i in range(nmonte):
        vlos4_samp = vlos4med + \
            np.random.normal(0.0,vlos4err,len(vlos4med))
        if (zero == 'no'):
            #Fit straight line in log-space for this draw:
            pfits_powline = fit_powline(rbin_tmp,vlos4_samp,vlos4err)
            
            #Draw router and gamouts:
            router = np.random.random()*(router_max-router_min)+router_min
            gamout = np.random.random()*(gamout_max-gamout_min)+gamout_min            

        tvl4 = tvl4func(rint,rbin_tmp,vlos4_samp,\
                        pfits_powline[0],pfits_powline[1],pfits_powline[2],router,gamout,zero)
        vs1_samp[i] = integrator(tvl4*test_surfden*rint,rint)
        vs2_samp[i] = integrator(tvl4*test_surfden*rint**3.0,rint)

    median, sixlow, sixhigh, ninelow, ninehigh,\
        nineninehigh, nineninelow = calcmedquartnine(vs1_samp)
    vs1imperr = (sixhigh-sixlow)/2.0
    median, sixlow, sixhigh, ninelow, ninehigh,\
        nineninehigh, nineninelow = calcmedquartnine(vs2_samp)
    vs2imperr = (sixhigh-sixlow)/2.0

    print 'VirialShape vs1:', vs1imp,vs1imperr
    print 'VirialShape vs2:', vs2imp,vs2imperr

    vs1bin = vs1imp
    vs2bin = vs2imp
    vs1err = vs1imperr
    vs2err = vs2imperr

    #Output also 2nd moment (pruning any would-be NaN values):
    rbin_kin = rbin_tmp[vlos2med > 0]
    sigpmean = np.sqrt(vlos2med[vlos2med > 0])
    sigperr = vlos2err[vlos2med > 0]/2.0/\
        np.sqrt(vlos2med[vlos2med > 0])

    #And finally calculate (for comparison only) the
    #Richardson & Fairbairn estimators:
    zeta_A = np.float(len(R))*np.sum(vz**4.0)/np.sum(vz**2.0)**2.0
    zeta_B = np.float(len(R))**2.0*np.sum(vz**4.0*R**2.0)/\
        (np.sum(vz**2.0)**2.0*np.sum(R**2.0))
    print 'Richardson+Fairbairn estimators:'
    print 'Nstars, zeta_A, zeta_B', len(R), zeta_A, zeta_B

    return rbin_kin, sigpmean, sigperr, \
        vs1bin, vs2bin, vs1err, vs2err

### MAIN CODE ###

#Test surface density fit and/or other code tests:
testthreeplum = 'no'
testmode = 'no'

#Import libraries: 
import numpy as np
if (testthreeplum == 'no' and testmode == 'no'):
    #Forbid plots to screen so GravSphere can run
    #remotely:
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pylab as plt
from scipy.integrate.quadrature import simps as integrator
from scipy.integrate.quadrature import simps
from scipy.misc.common import derivative
import emcee
import corner 
from matplotlib import rcParams
import lmfit as lm
import sys
import h5py
import scipy.integrate as si
from colossus.cosmology import cosmology 
from colossus.halo.concentration import concentration as colossus_cNFW
from colossus.halo.mass_defs import changeMassDefinition
cosmo = cosmology.setCosmology('planck18')
mX_keV = 500.0
if (mX_keV > 400.):
    cosmo_cfunc = lambda M200, h: cosmo_cfunc_DJ(M200,h,1.0e4)
else:
    cosmo_cfunc = lambda M200, h: cosmo_cfunc_DJ(M200,h,mX_keV)

#Welcome blurb: 
print '###### GRAVSPHERE VERSION 1.0 ######\n'

#Constants: 
G = 6.67e-11
Msun = 1.989e30
kpc = 3.086e19
kms = 1000.
Guse = G * Msun / kpc
deg = 360./2./np.pi
arcsec = 360./2./np.pi * 60. * 60.
arcmin = arcsec / 60.
Msunkpc3toGeVcm3 = 3.7973132915271756e-08
kpccm = 1e3*3.0856775*1e18
alpha_Jfac_deg = 0.5
rhocrit = 135.05
oden = 200
h = 0.7
year = 365.0*24.0*60.0*60.0

### CODE PARAMETERS ###

#Plot parameters: 
figsize = 8
figx = figsize
figy = figsize
myfontsize = 35
myfontsize2 = 50
mylegendfontsize = 25
mylinewidth = 5
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#Default parameters: 
datadir = '/vol/ph/astro_temp/jread/'
data_file_type = 'mock'
overtrue = 'no'
plotminmax = 'no'
plotlegend = 'no'
plotlegend_mass = 'no'
reduce_xticks = 'no'
walker_compare = 'no'
plot_evolve = 'no'
#nwalkers = 1000
#nmodels = 5000
#nwalkers = 750
#nmodels = 4000
nwalkers = 1000
nmodels = 5000
rmin = -1.0
rmax = -1.0
testmorehighbins = 'no'
singleinbin = 'no'
singlesplit = 'single'
multimode = 'normal'
betmode = 'betstar'
betprofile = 'osipkov'
blobstart = 'no'
osipkovstart = 'no'
linbetr0 = 'no'
ymin_Sigstar = 10
ymax_Sigstar = 1e5
propermotion = 'no'
virialshape = 'no'
usevs2 = 'yes'
zero = 'yes'
rotation = 'no'
sample_size = '1000'
subsample = 'no'
alp3sig = 0.0
p0in = np.array([1e3,1e3,1e3,0.1,0.1,0.1])
p0in_min = np.array([1e2,1e2,1e2,0.01,0.01,0.01])
p0in_max = np.array([1e6,1e6,1e6,1.0,1.0,1.0])
tracertol = 0.5
yMlow = 1e5
yMhigh = 5e9
rbin_inner = -1
rbin_outer = -1
Rhalf_use = -1
plotdenresults = 'yes'
loggcut = 'no'
craterphot = 'no'
mock_multimass = 'no'
Mstar = -1.0
Mstar_rlim = 5.0
maxdatfitrad = 5.0
maxdatrad = 5.0
data_file_kin = ''
colorpop1 = 'tomato'
colorpop2 = 'steelblue'
colorpop1data = 'red'
colorpop2data = 'blue'
usesimplemem = 'no'
usesimpcut = 'no'
Nbinin = -1
vsfac = 1.0
include_jardel = 'no'
walker_style_vs_old = 'no'
data_file_phot_vs = ''
cosmo_cprior = 'no'
rcmock = -1.0
sigmmock = -1.0
SIDM_Draco_priors = 'no'
addinfake = 'no'
use_rclinear = 'no'
use_rtlinear = 'no'
pfits = np.zeros(1)
calc_Jfac = 'no'
get_Juse = get_J
#get_Juse = get_J_NFW
prob_cut = ''
if (sample_size == '10000'):
    Nbin = 100
elif (sample_size == '100000'):
    Nbin = 316
elif (sample_size == '1000'):
    Nbin = 32
elif (sample_size == '100'):
    Nbin = 10
elif (sample_size == '1000-8'):
    Nbin = 8
elif (sample_size == '1000-16'):
    Nbin = 16
elif (sample_size == '1000-64'):
    Nbin = 64
elif (sample_size == '1000-128'):
    Nbin = 128

#Data parameters:
codemode = 'run'
makecorner = 'no'

#whichdata = 'PlumCuspIso'
#whichdata = 'PlumCoreIso'
#whichdata = 'NonPlumCuspOm'
#whichdata = 'NonPlumCoreOm'

#whichdata = 'PlumCoreOm'
#whichdata = 'PlumCuspOm'
#whichdata = 'NonPlumCoreIso'
#whichdata = 'NonPlumCuspIso'

#whichdata = 'PlumCoreTan'
#whichdata = 'PlumCuspTan'
#whichdata = 'NonPlumCoreTan'
#whichdata = 'NonPlumCuspTan'

#whichdata = 'SplitCompCusp_mgseplarge'
#whichdata = 'SplitCompCore'
#whichdata = 'TriaxCore'
#whichdata = 'TriaxCusp'
#myaxis = 'Z'

#whichdata = 'PlumTheiaSegCuspIso'
#whichdata = 'NonPlumTheiaSegCuspIso'
#whichdata = 'PlumTheiaSegCoreOm'
#whichdata = 'NonPlumTheiaSegCoreOm'
#whichdata = 'PlumTheiaUmaCuspIso'
#whichdata = 'NonPlumTheiaUmaCuspIso'
#whichdata = 'PlumTheiaUmaCoreOm'
#whichdata = 'NonPlumTheiaUmaCoreOm'

#whichdata = 'PlumTheiaDraCuspIso'
#whichdata = 'NonPlumTheiaDraCuspIso'
#whichdata = 'PlumTheiaDraCoreOm'
#whichdata = 'NonPlumTheiaDraCoreOm'

#whichdata = 'DracoMock'
#whichdata = 'DracoMockSplit'
#whichdata = 'CraMock'

#whichdata = 'Fornax'
#whichdata = 'Draco'
#whichdata = 'UMi'
#whichdata = 'Carina'
#whichdata = 'LeoI'
#whichdata = 'LeoII'
#whichdata = 'Sextans'
#whichdata = 'Sculptor'
#whichdata = 'SegI'
whichdata = 'CVnI'
#whichdata = 'LeoT'
#whichdata = 'Tucana'
#whichdata = 'Crater'
#whichdata = 'Aquarius'

#whichdata = 'NGC1407'
#whichdata = 'Ocen'

#whichdata = 'genina'

betprior = 'SIDM'

#This is the constant velocity error for the mock data. For non-constant
#or non-Gaussian errors, the binning routine needs updating!! 
verr = 2.0

if (whichdata == 'PlumCuspIso'):
    if (sample_size == '10000'):
        data_file = datadir+'PlumCuspIso/gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_10000_0_err.dat'
        outdirbase = datadir+'PlumCuspIso/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'PlumCuspIso/gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'PlumCuspIso/Output/1000/'
    elif (sample_size == '100'):
        sample_select = '0'
        data_file = datadir+'PlumCuspIso/gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_100_'+\
            sample_select+'_err.dat'
        outdirbase = datadir+'PlumCuspIso/Output/100-'+sample_select+'/'
    elif (sample_size == '1000-8'):
        data_file = datadir+'PlumCuspIso/gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'PlumCuspIso/Output/1000-8/'
    elif (sample_size == '1000-16'):
        data_file = datadir+'PlumCuspIso/gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'PlumCuspIso/Output/1000-16/'
    elif (sample_size == '1000-64'):
        data_file = datadir+'PlumCuspIso/gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'PlumCuspIso/Output/1000-64/'
    elif (sample_size == '1000-128'):
        data_file = datadir+'PlumCuspIso/gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'PlumCuspIso/Output/1000-128/'

    virialshape = 'yes'
    if (propermotion == 'yes'):
        virialshape = 'no'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 15
    SIDM_Draco_priors = 'nvary'                            
        
    #True solution (alp-bet-gam model): 
    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 25./100.*r0
    ra = 1e30

    if (sample_size != '10000'):
        rbin_inner = 0.0416106781007
        rbin_outer = 1.83198219733
        Rhalf_use = 0.27
    y_sigLOSmax = 25
    reduce_xticks = 'yes'

if (whichdata == 'CraMock'):
    whichmock = '2'
    if (whichmock == '1'):
        data_file = datadir+'Dwarfs/Mock//h378_544_0_err.dat'
        outdirbase = datadir+'Dwarfs/Output/CraMock1/'
        Mstar = 3092949.71972
    elif (whichmock == '2'):
        data_file = datadir+'Dwarfs/Mock/h636_145_0_err.dat'
        outdirbase = datadir+'Dwarfs/Output/CraMock2/'
        Mstar = 824389.8085990001
    data_file_kin = ''
    Mstar_rlim = 1.0
    Mstar_err = Mstar * 0.25
    mock_multimass = 'yes'

    p0in = np.array([50.0,1e-11,1e-11,0.1,0.1,0.1])
    p0in_min = np.array([10.0,1e-12,1e-12,0.001,0.001,0.001])
    p0in_max = np.array([1000.0,1e-9,1e-9,5.0,5.0,5.0])

    tracertol = 0.5
    maxdatrad = 5.0
    maxdatfitrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e-4
    ymax_Sigstar = 1000.0
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'PlumCoreIso'):
    if (sample_size == '10000'):
        data_file = datadir+'PlumCoreIso/gs010_bs050_rcrs100_rarcinf_core_0400mpc3_df_10000_0_err.dat'
        outdirbase = datadir+'PlumCoreIso/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'PlumCoreIso/gs010_bs050_rcrs100_rarcinf_core_0400mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'PlumCoreIso/Output/1000/'

    virialshape = 'yes'
    if (propermotion == 'yes'):
        virialshape = 'no'        
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 15
    SIDM_Draco_priors = 'nvary'                
        
    #True solution (alp-bet-gam model):
    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 100./100.*r0
    ra = 1e30

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    if (sample_size == '1000'):
        ymin_Sigstar = 0.1
        ymax_Sigstar = 1e3

if (whichdata == 'PlumCoreOm'):
    if (sample_size == '10000'):
        data_file = datadir+'PlumCoreOm/gs010_bs050_rcrs025_rarc100_core_0400mpc3_df_10000_0_err.dat'
        outdirbase = datadir+'PlumCoreOm/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'PlumCoreOm/gs010_bs050_rcrs025_rarc100_core_0400mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'PlumCoreOm/Output/1000/'

    virialshape = 'yes'
    if (propermotion == 'yes'):
        virialshape = 'no'        
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 15
    SIDM_Draco_priors = 'nvary'
        
    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 25./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    if (sample_size == '1000'):
        ymin_Sigstar = 0.1
        ymax_Sigstar = 1e3

if (whichdata == 'PlumCuspOm'):
    if (sample_size == '10000'):
        data_file = datadir+'PlumCuspOm/gs100_bs050_rcrs010_rarc100_cusp_0064mpc3_df_10000_0_err.dat'
        outdirbase = datadir+'PlumCuspOm/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'PlumCuspOm/gs100_bs050_rcrs010_rarc100_cusp_0064mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'PlumCuspOm/Output/1000/'
        
    virialshape = 'yes'
    if (propermotion == 'yes'):
        virialshape = 'no'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 15
    SIDM_Draco_priors = 'nvary'                        
        
    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 10./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    if (sample_size == '1000'):
        ymin_Sigstar = 0.1
        ymax_Sigstar = 1e3

if (whichdata == 'NonPlumCuspIso'):
    if (sample_size == '10000'):
        data_file = datadir+'NonPlumCuspIso/gs100_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_10000_0_err.dat'
        outdirbase = datadir+'NonPlumCuspIso/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'NonPlumCuspIso/gs100_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'NonPlumCuspIso/Output/1000/'

    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 25./100.*r0
    ra = 1e30

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    if (sample_size == '1000'):
        ymin_Sigstar = 0.1
        ymax_Sigstar = 1e3

if (whichdata == 'NonPlumCoreIso'):
    if (sample_size == '10000'):
        data_file = datadir+'NonPlumCoreIso/gs100_bs050_rcrs100_rarcinf_core_0400mpc3_df_10000_0_err.dat'
        outdirbase = datadir+'NonPlumCoreIso/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'NonPlumCoreIso/gs100_bs050_rcrs100_rarcinf_core_0400mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'NonPlumCoreIso/Output/1000/'

    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 100./100.*r0
    ra = 1e30

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    if (sample_size == '1000'):
        ymin_Sigstar = 0.1
        ymax_Sigstar = 1e3

if (whichdata == 'NonPlumCuspOm'):
    if (sample_size == '10000'):
        data_file = datadir+'NonPlumCuspOm/gs100_bs050_rcrs010_rarc100_cusp_0064mpc3_df_10000_0_err.dat'
        outdirbase = datadir+'NonPlumCuspOm/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'NonPlumCuspOm/gs100_bs050_rcrs010_rarc100_cusp_0064mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'NonPlumCuspOm/Output/1000/'
    elif (sample_size == '100'):
        data_file = datadir+'NonPlumCuspOm/gs100_bs050_rcrs010_rarc100_cusp_0064mpc3_df_100_0_err.dat'
        outdirbase = datadir+'NonPlumCuspOm/Output/100/'

    #Initial surface brightness guess:
    p0in = np.array([1e2,1e2,1e-11,0.1,0.1,0.1])
    p0in_min = np.array([1e1,1e1,1e-12,0.01,0.01,0.01])
    p0in_max = np.array([1e6,1e6,1e-10,3.0,3.0,3.0])

    #VSP pars:
    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'

    #True solution (alp-bet-gam model):
    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 10./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7

if (whichdata == 'NonPlumCoreOm'):
    if (sample_size == '10000'):
        data_file = datadir+'NonPlumCoreOm/gs100_bs050_rcrs025_rarc100_core_0400mpc3_df_10000_0_err.dat'
        outdirbase = datadir+'NonPlumCoreOm/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'NonPlumCoreOm/gs100_bs050_rcrs025_rarc100_core_0400mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'NonPlumCoreOm/Output/1000/'
    elif (sample_size == '100'):
        data_file = datadir+'NonPlumCoreOm/gs100_bs050_rcrs025_rarc100_core_0400mpc3_df_100_0_err.dat'
        outdirbase = datadir+'NonPlumCoreOm/Output/100/'
    elif (sample_size == '100000'):
        data_file = datadir+'NonPlumCoreOm/gs100_bs050_rcrs025_rarc100_core_0400mpc3_df_100000_0_err.dat'
        outdirbase = datadir+'NonPlumCoreOm/Output/100000/'
    elif (sample_size == '1000-128'):
        data_file = datadir+'NonPlumCoreOm/gs100_bs050_rcrs025_rarc100_core_0400mpc3_df_1000_0_err.dat'
        outdirbase = datadir+'NonPlumCoreOm/Output/1000-128/'

    #Initial surface brightness guess:
    p0in = np.array([100,100,1e-11,0.1,0.1,0.1])
    p0in_min = np.array([1e1,1e1,1e-12,0.01,0.01,0.01])
    p0in_max = np.array([1e6,1e6,1e-10,2.0,2.0,2.0])

    #VSP pars:
    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'

    #True solution (alp-bet-gam model):
    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 25./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 25
    ymin_Sigstar = 0.5
    ymax_Sigstar = 5e5

if (whichdata == 'PlumCoreTan'):
    if (sample_size == '10000'):
        data_file = datadir+'PlumCoreTan/data_c_rh4_rs175_gs01_ra0_b05n_10k_10000_0_err.dat'
        outdirbase = datadir+'PlumCoreTan/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'PlumCoreTan/data_c_rh4_rs175_gs01_ra0_b05n_10k_1000_0_err.dat'
        outdirbase = datadir+'PlumCoreTan/Output/1000/'

    p0in = np.array([1e3,1e3,1e3,0.1,0.1,0.1])
    p0in_min = np.array([1e2,1e2,1e2,0.01,0.01,0.01])
    p0in_max = np.array([1e6,1e6,1e6,2.0,2.0,2.0])

    overtrue = 'yes'
    rho0 = 3.021516E+07
    r0 = 4.0
    alp = 1.0
    bet = 4.0
    gam = 0.0
    rstar = 1.75
    betprofile = 'tangential'
    ra = -0.5

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    if (sample_size == '1000'):
        ymin_Sigstar = 0.1
        ymax_Sigstar = 1e3

if (whichdata == 'PlumCuspTan'):
    if (sample_size == '10000'):
        data_file = datadir+'PlumCuspTan/data_h_rh2_rs05_gs01_ra0_b05n_10k_10000_0_err.dat'
        outdirbase = datadir+'PlumCuspTan/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'PlumCuspTan/data_h_rh2_rs05_gs01_ra0_b05n_10k_1000_0_err.dat'
        outdirbase = datadir+'PlumCuspTan/Output/1000/'

    overtrue = 'yes'
    rho0 = 2.387329E+07
    r0 = 2.0
    alp = 1.0
    bet = 4.0
    gam = 1.0
    rstar = 0.5
    betprofile = 'tangential'
    ra = -0.5

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    if (sample_size == '1000'):
        ymin_Sigstar = 0.1
        ymax_Sigstar = 1e3
    if (sample_size == '10000'):
        rbin_inner = 0.0799391501894
        rbin_outer = 3.78593100838

if (whichdata == 'NonPlumCuspTan'):
    if (sample_size == '10000'):
        data_file = datadir+'NonPlumCuspTan/data_h_rh2_rs05_gs10_ra0_b05n_10k_10000_0_err.dat'
        outdirbase = datadir+'NonPlumCuspTan/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'NonPlumCuspTan/data_h_rh2_rs05_gs10_ra0_b05n_10k_1000_0_err.dat'
        outdirbase = datadir+'NonPlumCuspTan/Output/1000/'

    overtrue = 'yes'
    rho0 = 2.387329E+07
    r0 = 2.0
    alp = 1.0
    bet = 4.0
    gam = 1.0
    rstar = 0.5
    betprofile = 'tangential'
    ra = -0.5

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    if (sample_size == '1000'):
        ymin_Sigstar = 0.1
        ymax_Sigstar = 1e3

if (whichdata == 'NonPlumCoreTan'):
    if (sample_size == '10000'):
        data_file = datadir+'NonPlumCoreTan/data_c_rh4_rs175_gs10_ra0_b05n_10k_10000_0_err.dat'
        outdirbase = datadir+'NonPlumCoreTan/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'NonPlumCoreTan/data_c_rh4_rs175_gs10_ra0_b05n_10k_1000_0_err.dat'
        outdirbase = datadir+'NonPlumCoreTan/Output/1000/'

    p0in = np.array([1e3,1e3,1e3,0.1,0.1,0.1])
    p0in_min = np.array([1e2,1e2,1e2,0.01,0.01,0.01])
    p0in_max = np.array([1e6,1e6,1e6,2.0,2.0,2.0])

    overtrue = 'yes'
    rho0 = 3.021516E+07
    r0 = 4.0
    alp = 1.0
    bet = 4.0
    gam = 0.0
    rstar = 1.75
    betprofile = 'tangential'
    ra = -0.5

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    if (sample_size == '1000'):
        ymin_Sigstar = 0.1
        ymax_Sigstar = 1e3

if (whichdata == 'SplitCompCusp_mgseplarge'):
    if (sample_size == '10000'):
        data_file = datadir+'SplitCompCusp_mgseplarge/c1_100_050_050_100_cusp_c2_100_050_100_100_cusp_008_6d_mgseplarge.dat'
        outdirbase = datadir+'SplitCompCusp_mgseplarge/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'SplitCompCusp_mgseplarge/c1_100_050_050_100_cusp_c2_100_050_100_100_cusp_008_6d_mgseplarge.dat'
        outdirbase = datadir+'SplitCompCusp_mgseplarge/Output/1000/'

    #This chooses whether to draw perfect pops or use Mg:
    singlesplit = 'split'
    perfectdraw = 'yes'

    #True solution (alp-bet-gam model):
    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar1 = 50.*10./1000.
    ra1 = 100./100.*rstar1
    rstar2 = 100.*10./1000.
    ra2 = 100./100.*rstar2

    y_sigLOSmax = 25

if (whichdata == 'SplitCompCore'):
    if (sample_size == '10000'):
        data_file = datadir+'SplitCompCore/c1_100_050_050_100_core_c2_100_050_100_100_core_002_6d.dat'
        outdirbase = datadir+'SplitCompCore/Output/'
    elif (sample_size == '1000'):
        data_file = datadir+'SplitCompCore/c1_100_050_050_100_core_c2_100_050_100_100_core_002_6d.dat'
        outdirbase = datadir+'SplitCompCore/Output/1000/'

    singlesplit = 'split'
    perfectdraw = 'yes'

    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar1 = 50.*10./1000.
    ra1 = 100./100.*rstar1
    rstar2 = 100.*10./1000.
    ra2 = 100./100.*rstar2

    y_sigLOSmax = 35
    reduce_xticks = 'yes'

if (whichdata == 'DracoMockSplit'):
    use_unevol = 'no'
    data_file_kin = ''
    if (use_unevol == 'no'):
        Unevol = ''
        data_file_kin_vs1 = datadir+'Dwarfs/Mock/dracosplit_vel_sub1_fix.dat'
        data_file_phot1 = datadir+'Dwarfs/Mock/dracosplit_surfden_sub1_fix.dat'
        data_file_kin_vs2 = datadir+'Dwarfs/Mock/dracosplit_vel_sub2_fix.dat'
        data_file_phot2 = datadir+'Dwarfs/Mock/dracosplit_surfden_sub2_fix.dat'
        data_file_phot = datadir+'Dwarfs/Mock/dracosplit_surfden_comb_fix.dat'
    else:
        Unevol = 'Unevol/'
        data_file_kin_vs1 = datadir+'Dwarfs/Mock/dracosplit_vel_sub1_unevol.dat'
        data_file_phot1 = datadir+'Dwarfs/Mock/dracosplit_surfden_sub1_unevol.dat'
        data_file_kin_vs2 = datadir+'Dwarfs/Mock/dracosplit_vel_sub2_unevol.dat'
        data_file_phot2 = datadir+'Dwarfs/Mock/dracosplit_surfden_sub2_unevol.dat'
        data_file_phot = datadir+'Dwarfs/Mock/dracosplit_surfden_comb_unevol.dat'
    outdirbase = datadir+'Dwarfs/Output/DracoMockSplit/'+Unevol
    singlesplit = 'split'
    virialshape = 'yes'
    usevs2 = 'no'
    Nbin = np.sqrt(250.)
    data_file_type = 'walker_style'
    dgal = 76.0
    dgal_err = 6.0
    Rhalf_lum = -1
    Mstar = 0.29e6
    Mstar_rlim = 0.704
    Mstar_err = Mstar * 0.25

    if (use_unevol == 'yes'): 
        p0in = np.array([10.0,5.0,1.0,0.1,0.3,0.1])
        p0in_min = np.array([0.001,0.001,0.001,0.05,0.05,0.05])
        p0in_max = np.array([100.0,100.0,100.0,1.0,1.0,1.0])
        tracertol = 0.5
    else:
        p0in = np.array([10.0,10.0,10.0,0.1,0.3,0.5])
        p0in_min = np.array([-100.0,-100.0,-100.0,0.05,0.05,0.05])
        p0in_max = np.array([100.0,100.0,100.0,5.0,5.0,5.0])
        tracertol = 0.001
    maxdatrad = 5.0
    maxdatfitrad = 3.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e-4
    ymax_Sigstar = 1000.0
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'TriaxCore'):
    data_file = datadir+'Triaxial/StarsInCore3D.dat'
    if (sample_size == '10000'):
        outdirbase = datadir+'Triaxial/TriaxCore/Output/'
    elif (sample_size == '1000'):
        outdirbase = datadir+'Triaxial/TriaxCore/Output/1000/'
    data_file_type = 'triaxial'

    p0in = np.array([1e3,1e3,1e3,0.1,0.1,0.1])
    p0in_min = np.array([1e2,1e2,1e2,0.01,0.01,0.01])
    p0in_max = np.array([1e6,1e6,1e6,2.0,2.0,2.0])

    overtrue = 'yes'
    rho0 = 1.177e8
    r0 = 1.5
    alp = 1.0
    bet = 4.0
    gam = 0.23
    betprofile = 'triaxial'
    bet_eta = 0.5
    bet_bet0 = 0.0
    bet_betinf = 0.5
    bet_rs = 0.81

    y_sigLOSmax = 35
    reduce_xticks = 'yes'

if (whichdata == 'TriaxCusp'):
    data_file = datadir+'Triaxial/StarsInCusp3D.dat'
    if (sample_size == '10000'):
        outdirbase = datadir+'Triaxial/TriaxCusp/Output/'
    elif (sample_size == '1000'):
        outdirbase = datadir+'Triaxial/TriaxCusp/Output/1000/'
    data_file_type = 'triaxial'

    p0in = np.array([1e3,1e3,1e3,0.1,0.1,0.1])
    p0in_min = np.array([1e2,1e2,1e2,0.01,0.01,0.01])
    p0in_max = np.array([1e6,1e6,1e6,2.0,2.0,2.0])

    overtrue = 'yes'
    rho0 = 5.522e7
    r0 = 1.5
    alp = 1.0
    bet = 4.0
    gam = 1.0
    betprofile = 'triaxial'
    bet_eta = 0.5
    bet_bet0 = 0.0
    bet_betinf = 0.5
    bet_rs = 0.81

    y_sigLOSmax = 35
    reduce_xticks = 'yes'
    rbin_inner = 0.109620536488
    rbin_outer = 2.6063963536

if (whichdata == 'PlumTheiaSegCoreOm'):
    data_file = datadir+'Theia/seg1_gs010_bs050_rcrs025_rarc100_core_0400mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/PlumTheiaSegCoreOm/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 25./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 1.02020404081

if (whichdata == 'NonPlumTheiaSegCoreOm'):
    data_file = datadir+\
        'Theia/seg1_gs100_bs050_rcrs025_rarc100_core_0400mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/NonPlumTheiaSegCoreOm/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 25./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 0.232046409282

if (whichdata == 'PlumTheiaSegCuspIso'):
    data_file = datadir+'Theia/seg1_gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/PlumTheiaSegCuspIso/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 25./100.*r0
    ra = 1e30

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 0.270054010802

if (whichdata == 'NonPlumTheiaSegCuspIso'):
    data_file = datadir+'Theia/seg1_gs100_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/NonPlumTheiaSegCuspIso/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 25./100.*r0
    ra = 1e30

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 0.0900180036007

if (whichdata == 'PlumTheiaDraCoreOm'):
    data_file = datadir+'Theia/dra_gs010_bs050_rcrs025_rarc100_core_0400mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/PlumTheiaDraCoreOm/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 25./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 1.02020404081

if (whichdata == 'NonPlumTheiaDraCoreOm'):
    data_file = datadir+\
        'Theia/dra_gs100_bs050_rcrs025_rarc100_core_0400mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/NonPlumTheiaDraCoreOm/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 25./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 0.232046409282

if (whichdata == 'PlumTheiaDraCuspIso'):
    data_file = datadir+'Theia/dra_gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/PlumTheiaDraCuspIso/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 25./100.*r0
    ra = 1e30

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 0.270054010802

if (whichdata == 'NonPlumTheiaDraCuspIso'):
    data_file = datadir+'Theia/dra_gs100_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/NonPlumTheiaDraCuspIso/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 25./100.*r0
    ra = 1e30

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 0.0900180036007

if (whichdata == 'PlumTheiaUmaCoreOm'):
    data_file = datadir+'Theia/uma2_gs010_bs050_rcrs025_rarc100_core_0400mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/PlumTheiaUmaCoreOm/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 25./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 1.02020404081

if (whichdata == 'NonPlumTheiaUmaCoreOm'):
    data_file = datadir+'Theia/uma2_gs100_bs050_rcrs025_rarc100_core_0400mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/NonPlumTheiaUmaCoreOm/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 400./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 0.0
    rstar = 25./100.*r0
    ra = 100./100.*rstar

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 0.232046409282

if (whichdata == 'PlumTheiaUmaCuspIso'):
    data_file = datadir+'Theia/uma2_gs010_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/PlumTheiaUmaCuspIso/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 25./100.*r0
    ra = 1e30

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 0.270054010802

if (whichdata == 'NonPlumTheiaUmaCuspIso'):
    data_file = datadir+'Theia/uma2_gs100_bs050_rcrs025_rarcinf_cusp_0064mpc3_df_10000_0.bins'
    outdirbase = datadir+'Theia/NonPlumTheiaUmaCuspIso/'
    data_file_type = 'mock_binned'

    overtrue = 'yes'
    rho0 = 64./1000. * 1000.**3.
    r0 = 1.0
    alp = 1.0
    bet = 3.0
    gam = 1.0
    rstar = 25./100.*r0
    ra = 1e30

    y_sigLOSmax = 25
    ymin_Sigstar = 100
    ymax_Sigstar = 1e7
    Rhalf = 0.0900180036007

if (whichdata == 'Fornax'):
    data_file_kin = ''
    data_file_kin_vs = datadir+'Dwarfs/Final/for_justin1_spec.dat'
    data_file_phot_vs = datadir+'Dwarfs/Final/for_justin1_phot.dat'
    data_file_phot = ''
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/Fornax/'
    SIDM_Draco_priors = 'cosmo_c'
#    SIDM_Draco_priors = 'abundvary'
#    SIDM_Draco_priors = 'nvary'
#    SIDM_Draco_priors = 'NFW'
    
    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 25

    dgal = 138.0
    dgal_err = 8.0
    Rhalf_lum = 710./1000.
    Rhalf_lum_err = 77./1000.
    Mstar = 4.3e7
    Mstar_rlim = 1.6
    Mstar_err = Mstar * 0.25
    if (data_file_phot == ''):
        p0in = np.array([1000.0,1000.0,1000.0,0.1,0.1,0.1])
        p0in_min = np.array([-1e5,-1e5,-1e5,0.01,0.01,0.01])
        p0in_max = np.array([1e5,1e5,1e5,10.0,10.0,10.0])
        tracertol = 0.001
    else:
        p0in = np.array([-0.2,0.1,0.1,0.1,0.1,0.1])
        p0in_min = np.array([-8.0,0.001,0.001,0.01,0.01,0.01])
        p0in_max = np.array([0.0,4.0,4.0,4.0,4.0,4.0])
        tracertol = 0.001
    maxdatrad = 5.0
    maxdatfitrad = maxdatrad
    rotation = 'no'
    Arotmin = 0.0
    Arotmax = 2.0

    y_sigLOSmax = 25
    if (data_file_phot == ''):
        ymin_Sigstar = 1
        ymax_Sigstar = 1e6
    else:
        ymin_Sigstar = 1e-3
        ymax_Sigstar = 1.0

if (whichdata == 'DracoMock'):
    subsample = 'sidmbig'
    data_file_kin = ''
    data_file_phot = ''
    if (subsample == 'cusp'):
        data_file_kin_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_cusp_n1000_spec.dat'
        data_file_phot_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_cusp_n500_phot.dat'
        walker_style_vs_old = 'no'
    elif (subsample == 'cusp500'):
        data_file_kin_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_cusp_n500_spec.dat'
        data_file_phot_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_cusp_n500_phot.dat'
        walker_style_vs_old = 'no'
    elif (subsample == 'core'):
        data_file_kin_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_core_n1000_spec.dat'
        data_file_phot_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_core_n500_phot.dat'
        rcmock = 0.315
        sigmmock = 0.22
        walker_style_vs_old = 'no'
    elif (subsample == 'core500'):
        data_file_kin_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_core_n500_spec.dat'
        data_file_phot_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_core_n500_phot.dat'
        rcmock = 0.315
        sigmmock = 0.22
        walker_style_vs_old = 'no'
    elif (subsample == 'coreden500'):
        data_file_kin_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_coreden_n500_spec.dat'
        data_file_phot_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracosplit_matt_coreden_n500_phot.dat'
        walker_style_vs_old = 'no'
    elif (subsample == 'sidmbig'):
        data_file_kin_vs = datadir+'Dwarfs/Mock/RealDracoMock/dracoSIDMbig_vel.dat'
        data_file_phot = datadir+'Dwarfs/Mock/RealDracoMock/dracoSIDMbig_phot.dat'
        walker_style_vs_old = 'yes'
        rcmock = 1.0*1.75
        sigmmock = 130.0

#    SIDM_Draco_priors = 'yes'
    SIDM_Draco_priors = 'nvary'
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/DracoMock/'
    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    cosmo_cprior = 'no'
    Nbinin = 15
    if (subsample == 'core500' or subsample == 'sidmbig'):
        Nbinin = 10

    dgal = 82.0
    dgal_err = 6.0
    Rhalf_lum = -1
    Mstar = 0.29e6
    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25

    if (subsample == 'sidmbig'):
        p0in = np.array([1.0,1.0,1.0,0.1,0.1,0.1])
        p0in_min = np.array([1e-2,1e-2,1e-2,0.05,0.05,0.05])
        p0in_max = np.array([1e3,1e3,1e3,1.0,1.0,1.0])
        tracertol = 0.5
    else:
        p0in = np.array([1000.0,1000.0,1000.0,0.1,0.1,0.1])
        p0in_min = np.array([100,100,100,0.05,0.05,0.05])
        p0in_max = np.array([1e4,1e4,1e4,1.0,1.0,1.0])
        tracertol = 0.5
        
    maxdatrad = 5.0
    maxdatfitrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e-4
    ymax_Sigstar = 1000.0
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'Draco'):
    singlesplit = 'single'
    data_file_kin = ''

    if (singlesplit == 'single'):
        data_file_kin_vs = datadir+\
            'Dwarfs/Final/dra_justin1_spec.dat'
        data_file_phot_vs = datadir+\
            'Dwarfs/Final/dra_justin1_phot.dat'
        data_file_phot = ''
        outdirbase = datadir+'Dwarfs/Output/Draco/'
    else:
        data_file_kin_vs = datadir+\
            'Dwarfs/DraSplit/dra_rgbhb_wpvir_justin_specphot.dat'
        data_file_phot = datadir+'Dwarfs/dra_ih95.dat'
        data_file_phot1 = ''
        outdirbase = datadir+'Dwarfs/Output/DracoSplit/'

    data_file_type = 'walker_style_vs'
    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    cosmo_cprior = 'no'
    Nbinin = 15

    include_jardel = 'no'
    jardel_file = './Data/jardel_virus.txt'
#    SIDM_Draco_priors = 'abundvary'
#    SIDM_Draco_priors = 'nvary'
    SIDM_Draco_priors = 'cosmo_c'
    
    dgal = 76.0
    dgal_err = 6.0
    Rhalf_lum = 221./1000.
    Rhalf_lum_err = 19./1000.
    Mstar = 0.29e6

    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25
    p0in = np.array([1000.0,1000.0,1000.0,0.1,0.1,0.1])
    p0in_min = np.array([-1e5,-1e5,-1e5,0.01,0.01,0.01])
    p0in_max = np.array([1e5,1e5,1e5,15.0,15.0,15.0])
    tracertol = 0.001
    maxdatrad = 5.0
    maxdatfitrad = 5.0
    loggcut = 'no'

    y_sigLOSmax = 15
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e4
    yMhigh = 1e8
    vsfac = 1000.0
    plotdenresults = 'yes'

if (whichdata == 'Crater'):
    data_file_kin = ''
    data_file_kin_vs = datadir+'Dwarfs/cra2_justin_specphot.dat'
    data_file_phot = datadir+'Dwarfs/cra2counts.pop'
    craterphot = 'yes'
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/Crater/'

    dgal = 119.0
    dgal_err = 12.0
    Rhalf_lum = -1.0
    Rhalf_lum_err = 0.0
    Mstar = 3.3e5    #https://arxiv.org/pdf/1610.06189.pdf
    Mstar_rlim = 2.5
    Mstar_err = Mstar * 0.5
    p0in = np.array([250.0,\
                     250.0,\
                     250.0,\
                     1.0,\
                     1.0,\
                     1.0])
    p0in_min = np.array([100.0,100.0,100.0,0.1,0.1,0.1])
    p0in_max = np.array([1000.0,1000.0,1000.0,5.0,5.0,5.0])
    tracertol = 0.5
    maxdatrad = 5.0
    maxdatfitrad = 5.0
    y_sigLOSmax = 15
    ymin_Sigstar = 10.0
    ymax_Sigstar = 1000.0
    yMlow = 1e4
    yMhigh = 1e8
    plotdenresults = 'yes'

if (whichdata == 'UMi'):
    data_file_kin = ''
    data_file_phot = ''
    data_file_kin_vs = datadir+\
        'Dwarfs/Final/umi_justin1_spec.dat'
    data_file_phot_vs = datadir+\
        'Dwarfs/Final/umi_justin1_phot.dat'
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/UMi/'
#    SIDM_Draco_priors = 'abundvary'
#    SIDM_Draco_priors = 'nvary'
    SIDM_Draco_priors = 'cosmo_c'
        
    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 15

    dgal = 76.0
    dgal_err = 3.0
    Rhalf_lum = 181./1000.
    Rhalf_lum_err = 27./1000.
    Mstar = 0.29e6

    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25
    p0in = np.array([1000.0,1000.0,1000.0,0.1,0.1,0.1])
    p0in_min = np.array([1e-1,1e-1,1e-1,0.01,0.01,0.01])
    p0in_max = np.array([1e5,1e5,1e5,2.0,2.0,2.0])
    tracertol = 0.001
    maxdatrad = 5.0
    maxdatfitrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'Carina'):
    data_file_kin = ''
    data_file_phot = ''
    data_file_kin_vs = datadir+\
        'Dwarfs/Final/car_justin1_spec.dat'
    data_file_phot_vs = datadir+\
        'Dwarfs/Final/car_justin1_phot.dat'
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/Carina/'
#    SIDM_Draco_priors = 'abundvary'
#    SIDM_Draco_priors = 'nvary'
    SIDM_Draco_priors = 'cosmo_c'

    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 15

    dgal = 105.0
    dgal_err = 6.0
    Rhalf_lum = 250./1000.
    Rhalf_lum_err = 39./1000.
    Mstar = 0.38e6

    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25
    p0in = np.array([1000.0,1000.0,1000.0,0.1,0.1,0.1])
    p0in_min = np.array([-1e5,-1e5,-1e5,0.01,0.01,0.01])
    p0in_max = np.array([1e5,1e5,1e5,15.0,15.0,15.0])
    tracertol = 0.001
    maxdatrad = 5.0
    maxdatfitrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'LeoI'):
    data_file_kin = ''
    data_file_phot = ''
    data_file_kin_vs = datadir+\
        'Dwarfs/Final/leo1_justin1_spec.dat'
    data_file_phot_vs = datadir+\
        'Dwarfs/Final/leo1_justin1_phot.dat'
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/LeoI/'
#    SIDM_Draco_priors = 'abundvary'
#    SIDM_Draco_priors = 'nvary'
    SIDM_Draco_priors = 'cosmo_c'

    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 15

    dgal = 254.0
    dgal_err = 15.0
    Rhalf_lum = 251./1000.
    Rhalf_lum_err = 27./1000.
    Mstar = 5.5e6

    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25
    p0in = np.array([1000.0,1000.0,1000.0,0.1,0.1,0.1])
    p0in_min = np.array([-1e5,-1e5,-1e5,0.01,0.01,0.01])
    p0in_max = np.array([1e5,1e5,1e5,15.0,15.0,15.0])
    tracertol = 0.001
    maxdatrad = 5.0
    maxdatfitrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'LeoII'):
    data_file_kin = ''
    data_file_phot = ''
    data_file_kin_vs = datadir+\
        'Dwarfs/Final/leo2_justin1_spec.dat'
    data_file_phot_vs = datadir+\
        'Dwarfs/Final/leo2_justin1_phot.dat'
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/LeoII/'
#    SIDM_Draco_priors = 'abundvary'
#    SIDM_Draco_priors = 'nvary'
    SIDM_Draco_priors = 'cosmo_c'

    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 15

    dgal = 233.0
    dgal_err = 14.0
    Rhalf_lum = 176./1000.
    Rhalf_lum_err = 42./1000.
    Mstar = 0.74e6

    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25
    p0in = np.array([1000.0,1000.0,1000.0,0.1,0.1,0.1])
    p0in_min = np.array([-1e5,-1e5,-1e5,0.01,0.01,0.01])
    p0in_max = np.array([1e5,1e5,1e5,15.0,15.0,15.0])
    tracertol = 0.001
    maxdatrad = 5.0
    maxdatfitrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'Sextans'):
    data_file_kin = ''
    data_file_phot = ''
    data_file_kin_vs = datadir+\
        'Dwarfs/Final/sex_justin1_spec.dat'
    data_file_phot_vs = datadir+\
        'Dwarfs/Final/sex_justin1_phot.dat'
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/Sextans/'
#    SIDM_Draco_priors = 'abundvary'
#    SIDM_Draco_priors = 'nvary'
    SIDM_Draco_priors = 'cosmo_c'

    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 16

    dgal = 86.0
    dgal_err = 4.0
    Rhalf_lum = 695./1000.
    Rhalf_lum_err = 44./1000.
    Mstar = 0.44e6

    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25
    p0in = np.array([1000.0,1000.0,1000.0,0.1,0.1,0.1])
    p0in_min = np.array([-1e5,-1e5,-1e5,0.01,0.01,0.01])
    p0in_max = np.array([1e5,1e5,1e5,15.0,15.0,15.0])
    tracertol = 0.001
    maxdatrad = 5.0
    maxdatfitrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'Sculptor'):
    data_file_kin = ''
    data_file_phot = ''
    data_file_kin_vs = datadir+\
        'Dwarfs/Final/scl_justin1_spec.dat'
    data_file_phot_vs = datadir+\
        'Dwarfs/Final/scl_justin1_phot.dat'
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/Sculptor/'
#    SIDM_Draco_priors = 'abundvary'
#    SIDM_Draco_priors = 'nvary'
    SIDM_Draco_priors = 'cosmo_c'

    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 15

    dgal = 86.0
    dgal_err = 6.0
    Rhalf_lum = 283./1000.
    Rhalf_lum_err = 45./1000.
    Mstar = 2.3e6

    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25
    p0in = np.array([100.0,100.,100.,0.1,0.1,0.1])
    p0in_min = np.array([10.0,10.0,10.0,0.05,0.05,0.05])
    p0in_max = np.array([10000.0,10000.0,10000.0,2.0,2.0,2.0])
    tracertol = 0.5
    maxdatrad = 5.0
    maxdatfitrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'SegI'):
    data_file_kin = ''
    data_file_phot = datadir+\
                     'Dwarfs/SegIcut.txt'
    data_file_kin_vs = datadir+\
                       'Dwarfs/SegI_vel.txt'
    data_file_phot_vs = ''
    data_file_type = 'SegI'
    outdirbase = datadir+'Dwarfs/Output/SegI/'
#    SIDM_Draco_priors = 'abundvary'
    SIDM_Draco_priors = 'nvary'
    prob_cut = 'p09'
    
    virialshape = 'no'
    usevs2 = 'yes'
    zero = 'yes'
    if (prob_cut == 'p09'):
        Nbinin = 7
    else:
        Nbinin = 6

    dgal = 23.0
    Mstar = 1.0e3
    Mstar_rlim = 15.0
    Mstar_err = Mstar * 0.25
    Rhalf_lum = 29.0/1000.0
    p0in = np.array([0.1,0.1,0.1,1e-2,1e-2,1e-2])
    p0in_min = np.array([0.001,0.001,0.001,1e-3,1e-3,1e-3])
    p0in_max = np.array([10.0,10.0,10.0,0.1,0.1,0.1])
    tracertol = 0.5
    maxdatrad = 15.0
    maxdatfitrad = 0.04

    y_sigLOSmax = 25
    ymin_Sigstar = 1e-1
    ymax_Sigstar = 5e2
    yMlow = 1e4
    yMhigh = 1e8
    rbin_inner = 1e-3
    rbin_outer = 1.0

if (whichdata == 'CVnI'):
    data_file_kin = ''
    data_file_phot = ''
    data_file_kin_vs = datadir+\
                       'Dwarfs/Final/cvn1_justin1_spec.dat'
    data_file_phot_vs = datadir+\
                        'Dwarfs/Final/cvn1_justin1_phot.dat'
    data_file_type = 'walker_style_vs'
    outdirbase = datadir+'Dwarfs/Output/CVnI/'
    SIDM_Draco_priors = 'cosmo_c'
  
    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin = 6
    
    dgal = 216.0
    dgal_err = 8.0
    Rhalf_lum = 486./1000.
    Rhalf_lum_err = 14./1000.
    Mstar = 2.1e5
    
    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25
    p0in = np.array([1000.0,1000.0,1000.0,0.1,0.1,0.1])
    p0in_min = np.array([-1e5,-1e5,-1e5,0.01,0.01,0.01])
    p0in_max = np.array([1e5,1e5,1e5,15.0,15.0,15.0])
    tracertol = 0.001
    maxdatrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e4
    yMhigh = 1e8
    
if (whichdata == 'LeoT'):
    data_file_kin = datadir+'Dwarfs/leoT.txt'
    data_file_phot = datadir+'Dwarfs/leot_surface_bright.txt'
    data_file_type = 'LeoT'
    outdirbase = datadir+'Dwarfs/Output/LeoT/'

    Nbinin = 6
    singleinbin = 'no'
    virialshape = 'no'
    usevs2 = 'yes'
    zero = 'yes'

    dgal = 407.0
    dgal_err = 38.0
    Rhalf_lum = 0.106*1.68
    Rhalf_lum_err = 0.0
    Mstar = 0.0135*1e7

    Mstar_rlim = 0.704
    Mstar_err = Mstar * 0.25
    p0in = np.array([0.1,0.1,0.1,0.2,0.2,0.2])
    p0in_min = np.array([1e-5,1e-5,1e-5,0.001,0.001,0.001])
    p0in_max = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    tracertol = 0.5
    maxdatrad = 5.0
    maxdatfitrad = 0.8

    y_sigLOSmax = 25
    ymin_Sigstar = 1e-3
    ymax_Sigstar = 0.5
    yMlow = 1e4
    yMhigh = 1e8
    rbin_inner = 0.01
    rbin_outer = 1.0

if (whichdata == 'Tucana'):
    usemock = 'no'
    if (usemock == 'yes'):
        lowmedhi = 'high'
        mocknum = '0'
        samplenum = '36'
        usecore = ''
#        usecore = 'Core/'
        data_file_kin = datadir+'Dwarfs/Tucana_mock/'+usecore+'Tucana_'+\
                        lowmedhi+'_sample'+samplenum+'_foreground_'+\
                        mocknum+'_formatted.txt'
        data_file_phot = datadir+\
                         'Dwarfs/Tucana_mock/'+usecore+\
                         'Tucana_high_mockphot_2500sample_foreground50000.txt'
        outdirbase = datadir+'Dwarfs/Output/Tucana_mock_real_'+\
                     lowmedhi+'_'+samplenum+'_'+mocknum+'/'+usecore
        print 'Reading in kinematic data file:',\
            data_file_kin

        if (samplenum == '36'):
            Nbinin = 4
        else:
            Nbinin = 15
        Nbinin_phot = 15
        addinfake = 'no'
        fakehigh = 'high'
        singleinbin = 'no'
        virialshape = 'yes'
        usevs2 = 'yes'
        zero = 'no'

        Mstar = 0.56e6
        Mstar_rlim = 10.0
        Mstar_err = Mstar * 0.25
        p0in = np.array([1e2,1e2,1e2,0.2,0.2,0.2])
        p0in_min = np.array([1.0,1.0,1.0,1e-1,1e-1,1e-1])
        p0in_max = np.array([1e5,1e5,1e5,0.5,0.5,0.5])
        tracertol = 0.001
        maxdatrad = 10.0
        maxdatfitrad = 0.6
    else:
        data_file_kin = datadir+'Dwarfs/Tucana_vel.txt'
        data_file_phot = datadir+'Dwarfs/Tucana_phot.txt'
        outdirbase = datadir+'Dwarfs/Output/Tucana/'
        virialshape = 'yes'
        usevs2 = 'yes'
        zero = 'no'
        Nbinin = 4
        Nbinin_phot = 15
        addinfake = 'no'
        fakehigh = 'high'

        Mstar = 3.2e6
        Mstar_rlim = 10.0
        Mstar_err = Mstar * 0.25
        p0in = np.array([1e2,1e2,1e2,0.1,0.1,0.1])
        p0in_min = np.array([1.0,1.0,1.0,1e-2,1e-2,1e-2])
        p0in_max = np.array([5e4,5e4,5e4,0.7,0.7,0.7])
        tracertol = 0.001
        maxdatrad = 10.0
        maxdatfitrad = 0.7
                
    data_file_type = 'Tucana'       
    SIDM_Draco_priors = 'Tucana'    
    dgal = 887.0
    dgal_err = 49.0
    Rhalf_lum = 0.284
    Rhalf_lum_err = 0.054
    Mstar = 0.56e6
    Rhalf_lum = -1

    y_sigLOSmax = 35
    ymin_Sigstar = 1.0
    ymax_Sigstar = 1e4
    yMlow = 1e4
    yMhigh = 1e9
    rbin_inner = 0.1
    rbin_outer = 2.0
    
if (whichdata == 'Aquarius'):
    data_file_kin = datadir+'Dwarfs/Aquarius/aquarius_velocity_out.txt'
    data_file_phot = datadir+'Dwarfs/Aquarius/aquarius_phot.txt'
    data_file_type = 'Aquarius'
    outdirbase = datadir+'Dwarfs/Output/Aquarius/'
    singleinbin = 'yes'

    virialshape = 'no'
    usevs2 = 'yes'
    zero = 'yes'

    dgal = 977.0
    dgal_err = 45.0
    Rhalf_lum = 0.25*1.68
    Rhalf_lum_err = 0.0
    Mstar = (0.068+0.33)*1e7

    Mstar_rlim = 5.0
    Mstar_err = Mstar * 0.25
    p0in = np.array([1.0,1.0,1e-10,0.1,0.1,0.1])
    p0in_min = np.array([0.01,0.01,1e-11,0.01,0.01,0.01])
    p0in_max = np.array([10.0,10.0,1e-9,5.0,5.0,5.0])
    tracertol = 0.5
    maxdatrad = 5.0
    maxdatfitrad = 5.0

    y_sigLOSmax = 15
    ymin_Sigstar = 1e-3
    ymax_Sigstar = 0.5
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'Eridanus'):
    data_file_kin = datadir+'Dwarfs/Eridanus.txt'
    data_file_phot = datadir+'Dwarfs/eridanus_surface_bright.txt'
    data_file_type = 'LeoT'
    outdirbase = datadir+'Dwarfs/Output/Eridanus/'

    dgal = 366.0
    dgal_err = 17.0
    Rhalf_lum = 0.277
    Rhalf_lum_err = 0.014
    Mstar = 8.0e4

    Mstar_rlim = 0.704
    Mstar_err = Mstar * 0.25
    p0in = np.array([1.0,1.0,1.0,0.1,0.1,0.1])
    p0in_min = np.array([0.01,0.01,0.01,0.01,0.01,0.01])
    p0in_max = np.array([4.0,4.0,4.0,1.0,1.0,1.0])
    tracertol = 0.25
    maxdatrad = 5.0
    maxdatfitrad = 0.6

    y_sigLOSmax = 15
    ymin_Sigstar = 1e-3
    ymax_Sigstar = 0.5
    yMlow = 1e4
    yMhigh = 1e8

if (whichdata == 'NGC1407'):
    singlesplit = 'split'
    data_file_kin = ''
    data_file_kin_vs1 = datadir+'Sluggs/NGC1407_redGC_kin.dat'
    data_file_kin_vs2 = datadir+'Sluggs/NGC1407_blueGC_kin.dat'
    data_file_phot1 = datadir+'Sluggs/NGC1407_redGCs.txt'
    data_file_phot2 = datadir+'Sluggs/NGC1407_blueGCs.txt'
    data_file_phot = datadir+'Sluggs/NGC1407_stars.txt'
    data_file_type = 'sluggs'
    outdirbase = datadir+'Sluggs/Output/NGC1407/'
    virialshape = 'yes'
    usevs2 = 'yes'

    dgal = 28.05*1000.0
    dgal_err = dgal*0.1
    Rhalf_lum1 = -1.0
    Rhalf_lum2 = -1.0
    Rhalf_lum_err = 0.0
    Mstar = 10.0**11.8
    Mstar_rlim = 200.0
    Mstar_err = Mstar * 0.75
    p0in = np.array([500.0,500.0,500.0,10.0,10.0,10.0])
    p0in_min = np.array([100.0,100.0,100.0,1.0,1.0,1.0])
    p0in_max = np.array([100000.0,100000.0,100000.0,100.0,100.0,100.0])
    p0starin = np.array([100,100,100,10.0,10.0,10.0])
    p0starin_min = np.array([1e1,1e1,1e1,0.1,0.1,0.1])
    p0starin_max = np.array([1e4,1e4,1e4,100.0,100.0,100.0])
    tracertol = 0.5
    maxdatrad = 500.0
    maxdatfitrad = 500.0
    y_sigLOSmax = 350.0
    ymin_Sigstar = 1e-2
    ymax_Sigstar = 1e2
    yMlow = 1e9
    yMhigh = 1e13
    plotdenresults = 'yes'

if (whichdata == 'Ocen'):
    data_file_kin = datadir+'Omega_cen/Ocen_sigLOS.txt'
    data_file_phot = datadir+'Omega_cen/Ocen_surfden.txt'
    data_file_pmr = datadir+'Omega_cen/Ocen_sigPMR.txt'
    data_file_pmt = datadir+'Omega_cen/Ocen_sigPMT.txt'
    data_file_type = 'Ocen'
    outdirbase = datadir+'Omega_cen/Output/'

    SIDM_Draco_priors = 'Ocen'
    virialshape = 'no'
    propermotion = 'yes'
    
    Mstar = 4.0e6
    Mstar_err = Mstar*0.75
    Mstar_rlim = 50.0
    rbin_inner = 1e-5
    rbin_outer = 2.0
    
    dgal = 5.2
    Rhalf_lum = 4.8*dgal/arcmin
    if (dgal == 5.2):
        p0in = np.array([1e-3,1e-3,1e-3,1.5e-2,1.5e-2,1.5e-2])
        p0in_min = np.array([2e-7,2e-7,2e-7,7e-5,7e-5,7e-5])
        p0in_max = np.array([1,5e-7,5e-7,2e-1,2e-1,2e-1])
    elif (dgal == 4.8):
        p0in = np.array([1e-3,1e-3,1e-3,2e-2,2e-2,2e-2])
        p0in_min = np.array([2e-7,2e-7,2e-7,7e-5,7e-5,7e-5])
        p0in_max = np.array([1,1,1,3e-1,3e-1,3e-1])
    else:
        p0in = np.array([1e-2,1e-2,1e-2,1.5e-2*dgal/5.2,1.5e-2*dgal/5.2,1.5e-2*dgal/5.2])
        p0in_min = np.array([1e-7,1e-7,1e-7,5e-5*dgal/5.2,5e-5*dgal/5.2,5e-5*dgal/5.2])
        p0in_max = np.array([1e1,1e1,1e1,2e-1*dgal/5.2,2e-1*dgal/5.2,2e-1*dgal/5.2])
    tracertol = 0.5
    maxdatrad = 500.0
    maxdatfitrad = 500.0

    y_sigLOSmax = 25.0
    ymin_Sigstar = 1e-4
    ymax_Sigstar = 5
    yMlow = 1e-3
    yMhigh = 1e7
    plotdenresults = 'yes'
    
if (whichdata == 'genina'):
    genina_gal_num = 28
    if (genina_gal_num >= 0):
        data_file = datadir+'Dwarfs/Genina/Galaxy_%d.hdf5' % \
                    (genina_gal_num)
        outdirbase = datadir+'Dwarfs/Output/Genina/Galaxy_%d/' % \
                     (genina_gal_num)
    else:
        data_file = datadir+'Dwarfs/Genina/V1_9_2.hdf5'
        outdirbase = datadir+'Dwarfs/Output/Genina/'                    
    
    data_file_type = 'genina'
    SIDM_Draco_priors = 'Genina'
    virialshape = 'yes'
    usevs2 = 'yes'
    zero = 'no'
    Nbinin_phot = 15
    Nbinin_kin = 15

    Mstar_rlim = 50.0
    rbin_inner = 0.1
    rbin_outer = 3.0    
    if (genina_gal_num == 25 or \
        genina_gal_num == 26 or \
        genina_gal_num == 27 or \
        genina_gal_num == 28):
        p0in = np.array([100,100,100,0.25,0.5,0.75])
        p0in_min = np.array([10,10,10,0.1,0.1,0.1])
        p0in_max = np.array([1e4,1e4,1e4,20.0,20.0,20.0])
        SIDM_Draco_priors = 'Genina2'
        rbin_inner = 0.1
        rbin_outer = 15.0
    else:
        p0in = np.array([100,100,100,0.1,0.5,0.75])
        p0in_min = np.array([10,10,10,0.05,0.05,0.05])
        p0in_max = np.array([1e4,1e4,1e4,2.0,2.0,2.0])
    tracertol = 0.5
    maxdatrad = 50.0
    maxdatfitrad = 50.0
    
    y_sigLOSmax = 45
    ymin_Sigstar = 1e1
    ymax_Sigstar = 1e5
    yMlow = 1e6
    yMhigh = 1e10
    
#This for overlaying the true solutions (mock data): 
ranalmin = 0.001
ranalmax = 500.
ranalpnts = 10000
ranal = np.logspace(np.log10(ranalmin),np.log10(ranalmax),ranalpnts)

#Write the propermotion results to a separate folder:
if (usesimplemem == 'yes'):
    outdirbase = outdirbase + 'Simple/'
    if (usesimpcut == 'yes'):
        outdirbase = outdirbase + 'Cut/'
    elif (usesimpcut == 'oz'):
        outdirbase = outdirbase + 'Oz/'
if (propermotion == 'yes'):
    outdirbase = outdirbase + 'Propermotion/'
if (rotation == 'yes'):
    outdirbase = outdirbase + 'Rotation/'
if (virialshape == 'yes'):
    outdirbase = outdirbase + 'VirialShape/'
    if (usevs2 == 'no'):
        outdirbase = outdirbase + 'Novs2/'
    if (zero == 'no'):
        outdirbase = outdirbase + 'Vlos_interp/'
if (Nbinin > 0):
    outdirbase = outdirbase + 'Nbinin/'
if (include_jardel == 'yes'):
    outdirbase = outdirbase + 'Jar/'
if (addinfake == 'yes'):
    outdirbase = outdirbase + 'AddFake/'
    if (fakehigh == 'high'):
        outdirbase = outdirbase + 'High/'
    elif (fakehigh == 'mid'):
        outdirbase = outdirbase + 'Mid/'
    elif (fakehigh == 'low'):
        outdirbase = outdirbase + 'Low/'                
                
#Split pop with rotation not supported (yet): 
if (rotation == 'yes' and singlesplit == 'split'):
    print 'Split pop with rotation not supported. Sorry.'
    sys.exit(0)

#Priors and model parameters: 
if (linbetr0 == 'yes'):
    betr0min = 0.1
    betr0max = 1.0
else:
    betr0min = -2
    betr0max = 0.0
betnmin = 1.0
betnmax = 10.0
if (betprior == 'SIDM'):
    outdir = outdirbase + 'SIDM/'
    if (betmode == 'betstar'):
        bet0min = -1.0
        bet0max = 1.0
        betinfmin = -1.0
        betinfmax = 1.0
    else:
        bet0min = -5.0
        bet0max = 1.0
        betinfmin = -5.0
        betinfmax = 1.0
    betnmin = 1.0
    betnmax = 3.0
elif (betprior == 'SIDM2'):
    outdir = outdirbase + 'SIDM2/'
    if (betmode == 'betstar'):
        bet0min = -0.1
        bet0max = 0.1
        betinfmin = -0.2
        betinfmax = 0.4
    else:
        bet0min = -5.0
        bet0max = 1.0
        betinfmin = -5.0
        betinfmax = 1.0
    betnmin = 1.0
    betnmax = 3.0
        
if (singleinbin == 'yes'):
    outdir = outdir + 'Singleinbin/'
if (testmorehighbins == 'yes'):
    outdir = outdir + 'Testmorehighbins/'
if (cosmo_cprior == 'yes'):
    outdir = outdir + 'CosmoC/'
elif (subsample == 'cusp'):
    outdir = outdir + 'Cusp/'
elif (subsample == 'cusp500'):
    outdir = outdir + 'Cusp500/'
elif (subsample == 'core'):
    outdir = outdir + 'Core/'
elif (subsample == 'core500'):
    outdir = outdir + 'Core500/'
elif (subsample == 'coreden500'):
    outdir = outdir + 'Coreden500/'
elif (subsample == 'sidmbig'):
    outdir = outdir + 'SIDMbig/'
if (data_file_type == 'triaxial'):
    outdir = outdir + myaxis + '/'
if (SIDM_Draco_priors == 'no'):
    outdir = outdir + 'GenPrior/'
elif (SIDM_Draco_priors == 'Zavala'):
    outdir = outdir + 'Zavala/'
elif (SIDM_Draco_priors == 'nvary'):
    outdir = outdir + 'nvary/'    
elif (SIDM_Draco_priors == 'Ocen'):
    outdir = outdir + 'Ocen/'
elif (SIDM_Draco_priors == 'Tucana'):
    outdir = outdir + 'Tucana/'
elif (SIDM_Draco_priors == 'Genina'):
    outdir = outdir + 'Genina/'
elif (SIDM_Draco_priors == 'Genina2'):
    outdir = outdir + 'Genina2/'
elif (SIDM_Draco_priors == 'abund'):
    outdir = outdir + 'Abund/'
elif (SIDM_Draco_priors == 'abundvary'):
    outdir = outdir + 'Abund_nvary/'
elif (SIDM_Draco_priors == 'NFW'):
    outdir = outdir + 'NFW/'
elif (SIDM_Draco_priors == 'cosmo_c'):
    if (mX_keV > 400.0):
        outdir = outdir + 'Cosmo_c/'
    else:
        outdir = outdir + 'Cosmo_c_' + np.str(mX_keV) +'/'
if (prob_cut == 'p09'):
    outdir = outdir + 'p09/'

n_betpars = 4
if (rotation == 'yes'):
    n_betpars = 5

nu = multiplumden
nu_components = 3
Sigfunc = multiplumsurf
n_nupars = nu_components * 2

if (singlesplit == 'single'):
    if (propermotion == 'no'):
        if (virialshape == 'no'):
            lnlike = lnlike_single
            lnprob = lnprob_single
            lnprior_set = lnprior_set_single
        else:
            if (rotation == 'no'):
                lnlike = lnlike_single_vs
                lnprob = lnprob_single_vs
                lnprior_set = lnprior_set_single
            else:
                lnlike = lnlike_single_vs
                lnprob = lnprob_single_vs
                lnprior_set = lnprior_set_single_rot
    elif (propermotion == 'yes'):
        lnlike = lnlike_single_prop
        lnprob = lnprob_single_prop
        lnprior_set = lnprior_set_single
elif (singlesplit == 'split'):
    if (virialshape == 'no'):
        lnlike = lnlike_split
        lnprob = lnprob_split
        lnprior_set = lnprior_set_split
    else:
        lnlike = lnlike_split_vs
        lnprob = lnprob_split_vs
        lnprior_set = lnprior_set_split

if (rbin_inner > 0):
    print 'Setting plot range:', rbin_inner, rbin_outer

### READ IN AND BIN THE DATA ###

#Read in the data: 
if (data_file_type == 'mock'):
    f = open(data_file,'r')
    data = np.genfromtxt(f)
    f.close()

    #Bin the data:
    if (singlesplit == 'single'):
        R = np.sqrt(data[:,0]**2.0 + data[:,1]**2.0)
        vz = data[:,5]
        if (mock_multimass == 'no'):
            ms = np.zeros(len(R)) + 1.0
        else:
            ms = data[:,6]/np.sum(data[:,6])*np.float(len(data[:,6]))
            if (Nbinin < 0):
                Nbin = np.sqrt(np.sum(ms))
            else:
                Nbin = Nbinin
            print 'Number of stars per bin:', Nbin
        if (propermotion == 'no'):
            rbin, surfden, surfdenerr, sigpmean, \
                sigperr, Rhalf, vmean, Mstar_rad, Mstar_prof = \
                binthedata(R,vz,ms,Nbin,verr)
        elif (propermotion == 'yes'):
            x = data[:,0]
            y = data[:,1]
            z = data[:,2]
            vx = data[:,3]
            vy = data[:,4]
            vz = data[:,5]
            rbin, surfden, surfdenerr, sigpmean, sigperr, \
                sigpmt, sigpmterr, sigpmr, sigpmrerr, Rhalf = \
                binthedata_prop(x,y,z,vx,vy,vz,ms,Nbin,verr)
        rbin_phot = rbin
        rbin_kin = rbin
        print 'Half stellar mass radius:', Rhalf
    elif (singlesplit == 'split'):
        if (perfectdraw == 'yes'):
            R1 = np.sqrt((data[data[:,20] == 1,0]/1000.)**2.0 + \
                             (data[data[:,20] == 1,1]/1000.)**2.0)
            vz1 = data[data[:,20] == 1,11]
            ms1 = np.zeros(len(R1)) + 1.0
            R2 = np.sqrt((data[data[:,20] == 2,0]/1000.)**2.0 + \
                             (data[data[:,20] == 2,1]/1000.)**2.0)
            vz2 = data[data[:,20] == 2,11]
            ms2 = np.zeros(len(R2)) + 1.0
            if (sample_size == '1000'):
                #Down sample by a factor 10:
                select = np.random.randint(len(R1), \
                                           size=len(R1)/10)
                R1 = R1[select]
                vz1 = vz1[select]
                ms1 = ms1[select]
                select = np.random.randint(len(R2), \
                                           size=len(R2)/10)
                R2 = R2[select]
                vz2 = vz2[select]
                ms2 = ms2[select]
            else:
                Mg = data[:,13]
                #*** WARNING *** Code needs to be written here!
        rbin1, surfden1, surfdenerr1, sigpmean1, \
            sigperr1, Rhalf1, vmean1, Mstar_rad, Mstar_prof = \
            binthedata(R1,vz1,ms1,Nbin,verr)
        rbin2, surfden2, surfdenerr2, sigpmean2, \
            sigperr2, Rhalf2, vmean2, Mstar_rad, Mstar_prof = \
            binthedata(R2,vz2,ms2,Nbin,verr)
        rbin1_phot = rbin1
        rbin2_phot = rbin2
        rbin1_kin = rbin1
        rbin2_kin = rbin2
        rbin = rbin1
        rbin_phot = rbin1_phot
        rbin_kin = rbin1_kin
        print 'Half stellar mass radius 1:', Rhalf1
        print 'Half stellar mass radius 2:', Rhalf2

    #Set the stellar mass to zero for the mocks [if Mstar < 0].
    #Otherwise this takes the stellar distribution from the tracer
    #surface density profile:
    if (Mstar < 0):
        Mstar_rad = rbin
        Mstar_prof = np.zeros(len(Mstar_rad))
        Mstar = 0.0
        Mstar_err = 1.0

    #Overwrite Rhalf if req.:
    if (Rhalf_use > 0):
        Rhalf = Rhalf_use

    #Set up vzerr array (for VirialShape):
    vzerr = np.zeros(len(vz))+verr
elif (data_file_type == 'mock_binned'):
    f = open(data_file,'r')
    data = np.genfromtxt(f)
    f.close()
    
    rbin = data[:,2]
    surfden = data[:,6]
    surfdenerr = data[:,7]
    sigpmean = data[:,14]
    sigperr = data[:,15]
    sigpmt = data[:,18]
    sigpmterr = data[:,19]
    sigpmr = data[:,16]
    sigpmrerr = data[:,17]

    Mstar_rad = rbin
    Mstar_prof = np.zeros(len(Mstar_rad))
    Mstar = 0.0
    Mstar_err = 1.0

    rbin_phot = rbin
    rbin_kin = rbin

elif (data_file_type == 'triaxial'):
    f = open(data_file,'r')
    data = np.genfromtxt(f)
    f.close()

    #Down sample the data:
    if (sample_size == '1000'):
        mysize = 1000
    elif (sample_size == '10000'):
        mysize = 10000
    select = np.random.randint(len(data[:,0]), \
                                   size=mysize)

    #Bin the data:
    if (myaxis == 'Z'):
        R = np.sqrt(data[select,0]**2.0 + data[select,1]**2.0)
        vz = data[select,5]
    elif (myaxis == 'X'):
        R = np.sqrt(data[select,1]**2.0 + data[select,2]**2.0)
        vz = data[select,3]
    elif (myaxis == 'Y'):
        R = np.sqrt(data[select,0]**2.0 + data[select,2]**2.0)
        vz = data[select,4]
    ms = np.zeros(len(R)) + 1.0
    if (propermotion == 'no'):
        rbin, surfden, surfdenerr, sigpmean, \
            sigperr, Rhalf, vmean, Mstar_rad, Mstar_prof = \
            binthedata(R,vz,ms,Nbin,verr)
    elif (propermotion == 'yes'):
        if (myaxis == 'Z'):
            x = data[select,0]
            y = data[select,1]
            z = data[select,2]
            vx = data[select,3]
            vy = data[select,4]
            vz = data[select,5]
        elif (myaxis == 'X'):
            x = data[select,2]
            y = data[select,1]
            z = data[select,0]
            vx = data[select,5]
            vy = data[select,4]
            vz = data[select,3]
        elif (myaxis == 'Y'):
            x = data[select,0]
            y = data[select,2]
            z = data[select,1]
            vx = data[select,3]
            vy = data[select,5]
            vz = data[select,4]
        rbin, surfden, surfdenerr, sigpmean, sigperr, \
            sigpmt, sigpmterr, sigpmr, sigpmrerr, Rhalf = \
            binthedata_prop(x,y,z,vx,vy,vz,ms,Nbin,verr)
    rbin_phot = rbin
    rbin_kin = rbin
    print 'Half stellar mass radius:', Rhalf
    Mstar_rad = rbin
    Mstar_prof = np.zeros(len(Mstar_rad))
    Mstar = 0.0
    Mstar_err = 1.0
elif (data_file_type == 'walker_style' or data_file_type == 'walker_style_vs'):
    if (singlesplit == 'single'):
        if (data_file_kin != ''):
            f = open(data_file_kin,'r')
            data_kin = np.genfromtxt(f)
            f.close()
            rbin_kin = data_kin[:,0]*dgal/arcmin
            sigpmean = data_kin[:,4]
            sigperr = data_kin[:,5]

        if (data_file_phot != ''):
            if (craterphot == 'no'):
                rbin_phot, surfden, surfdenerr, \
                    rbin_photfit, surfdenfit, surfdenerrfit, \
                    Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits = \
                    walker_surf_sort(data_file_phot,Rhalf_lum)
            else:
                rbin_phot, surfden, surfdenerr, \
                    rbin_photfit, surfdenfit, surfdenerrfit, \
                    Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits = \
                    walker_surf_sort_cra(data_file_phot,Rhalf_lum)
        else:
            #Calculate surface density:
            f = open(data_file_phot_vs,'r')
            data_phot_vs = np.genfromtxt(f)
            f.close()
            R = data_phot_vs[:,4]*dgal/arcmin
            vz = np.zeros(len(R))
            vzerr = np.zeros(len(R))
            ms = data_phot_vs[:,10]
            if (Nbinin < 0):
                Nbin = np.sqrt(np.sum(ms))
            else:
                Nbin = Nbinin
            print 'Number of stars per photometric bin:', Nbin
            print 'Total effective no. of stars (photometric):', np.sum(ms)
            dummy_err = 2.0
            rbin_phot, surfden, surfdenerr, dummy_sigpmean, \
                dummy_sigperr, Rhalf, dummy_vmean, Mstar_rad, Mstar_prof = \
                binthedata(R,vz,ms,Nbin,dummy_err)
            rbin_photfit = rbin_phot[rbin_phot < maxdatrad]
            surfdenfit = surfden[rbin_phot < maxdatrad]
            surfdenerrfit = surfdenerr[rbin_phot < maxdatrad]
            pfits = tracerfit(p0in,p0in_min,p0in_max,\
                              rbin_photfit,surfdenfit,surfdenerrfit)
        
        #Half stellar mass radius:
        print 'Calculated projected half stellar mass radius:', Rhalf

        #For plotting:
        rbin = rbin_phot

        if (data_file_type == 'walker_style_vs'):
            if (walker_style_vs_old == 'yes'):
                f = open(data_file_kin_vs,'r')
                data_kin_vs = np.genfromtxt(f)
                f.close()
                gotvz = np.where(data_kin_vs[:,11])[0]
                R = data_kin_vs[gotvz,4]*dgal/arcmin
                vz = data_kin_vs[gotvz,10]
                vzerr = data_kin_vs[gotvz,11]
                ms = data_kin_vs[gotvz,19]
            else:   
                f = open(data_file_kin_vs,'r')
                data_kin_vs = np.genfromtxt(f)
                f.close()
                gotvz = np.where(data_kin_vs[:,6])[0]
                R = data_kin_vs[gotvz,4]*dgal/arcmin
                vz = data_kin_vs[gotvz,6]
                vzerr = data_kin_vs[gotvz,7]
                ms = data_kin_vs[gotvz,15]

            #Perform a hard cut on the membership probability:
            memcut = 0.0
            R = R[ms > memcut]
            vz = vz[ms > memcut]
            vzerr = vzerr[ms > memcut]
            ms = ms[ms > memcut]

            #Perform a cut on logg:
            if (loggcut == 'yes'):
                logg_cut = 3.7
                logg = data_kin_vs[gotvz,10]
                logg = logg[ms > memcut]
                R = R[logg < logg_cut]
                vz = vz[logg < logg_cut]
                vzerr = vzerr[logg < logg_cut]
                ms = ms[logg < logg_cut]

            #Set auto-binning after cuts:
            if (Nbinin < 0):
                Nbin = np.sqrt(np.sum(ms))
            else:
                Nbin = Nbinin
            print 'Number of stars per bin (kinematic):', Nbin
            print 'Effective no. of stars (kinematic):', np.sum(ms)

            if (data_file_kin == ''):
                nmonte = 1000
                rbin_kin, sigpmean, sigperr, \
                    vs1bin, vs2bin, vs1err, vs2err = \
                    calc_virial_moments(Rhalf,nmonte,R,vz,vzerr,ms,\
                                        pfits)
    else:
        if (data_file_phot1 == ''):
            f = open(data_file_kin_vs,'r')
            data_kin_vs = np.genfromtxt(f)
            f.close()
            R1 = data_kin_vs[:,4]*dgal/arcmin
            vz1 = data_kin_vs[:,10]
            vzerr1 = data_kin_vs[:,11]
            R2 = data_kin_vs[:,4]*dgal/arcmin
            vz2 = data_kin_vs[:,10]
            vzerr2 = data_kin_vs[:,11]
            ms1 = data_kin_vs[:,19]
            ms2 = data_kin_vs[:,24]
            Nbin1 = np.sqrt(np.sum(ms1))
            Nbin2 = np.sqrt(np.sum(ms2))
            print 'Nbin1 (photometric):', Nbin1
            print 'Nbin2 (photometric):', Nbin2
            dummy_err = 2.0
            rbin1_phot, surfden1, surfdenerr1, sigpmean, \
                sigperr, Rhalf1, vmean, Mstar_rad1, Mstar_prof1 = \
                binthedata(R1,vz1,ms1,Nbin1,dummy_err)
            rbin2_phot, surfden2, surfdenerr2, sigpmean, \
                sigperr, Rhalf2, vmean, Mstar_rad2, Mstar_prof2 = \
                binthedata(R2,vz2,ms2,Nbin2,dummy_err)

            rbin_photfit = rbin1_phot[rbin1_phot < maxdatrad]
            surfdenfit = surfden1[rbin1_phot < maxdatrad]
            surfdenerrfit = surfdenerr1[rbin1_phot < maxdatrad]
            pfits1 = tracerfit(p0in,p0in_min,p0in_max,\
                              rbin_photfit,surfdenfit,surfdenerrfit)
            rbin_photfit = rbin2_phot[rbin2_phot < maxdatrad]
            surfdenfit = surfden2[rbin2_phot < maxdatrad]
            surfdenerrfit = surfdenerr2[rbin2_phot < maxdatrad]
            pfits2 = tracerfit(p0in,p0in_min,p0in_max,\
                              rbin_photfit,surfdenfit,surfdenerrfit)
            print 'Half stellar mass radius 1:', Rhalf1
            print 'Half stellar mass radius 2:', Rhalf2
        else:
            rbin1_phot, surfden1, surfdenerr1, \
                rbin_photfit1, surfdenfit1, surfdenerrfit1, \
                Mstar_rad1, Mstar_prof1, Mstar_surf1, Rhalf1, pfits1 = \
                walker_surf_sort(data_file_phot1,Rhalf_lum)
            rbin2_phot, surfden2, surfdenerr2, \
                rbin_photfit2, surfdenfit2, surfdenerrfit2, \
                Mstar_rad2, Mstar_prof2, Mstar_surf2, Rhalf2, pfits2 = \
                walker_surf_sort(data_file_phot2,Rhalf_lum)
            
            print 'Half stellar mass radius 1:', Rhalf1
            print 'Half stellar mass radius 2:', Rhalf2

        if (Mstar > 0.0):
            #Read in the photometric light profile for the 
            #stellar mass:
            rbin_phot, surfden, surfdenerr, \
                rbin_photfit, surfdenfit, surfdenerrfit, \
                Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits = \
                walker_surf_sort(data_file_phot,Rhalf_lum)
        else:
            Mstar_rad = rbin1_phot
            Mstar_prof = np.zeros(len(Mstar_rad))
            Mstar = 0.0
            Mstar_err = 1.0

        if (data_file_phot1 != ''):
            f = open(data_file_kin_vs1,'r')
            data_kin_vs = np.genfromtxt(f)
            f.close()
            gotvz = np.where(data_kin_vs[:,10])[0]
            R1 = data_kin_vs[gotvz,4]*dgal/arcmin
            vz1 = data_kin_vs[gotvz,10]
            vzerr1 = data_kin_vs[gotvz,11]
            ms1 = data_kin_vs[gotvz,19]
            Nbin1 = np.sqrt(np.sum(ms1))
            f = open(data_file_kin_vs2,'r')
            data_kin_vs = np.genfromtxt(f)
            f.close()
            gotvz = np.where(data_kin_vs[:,10])[0]
            R2 = data_kin_vs[gotvz,4]*dgal/arcmin
            vz2 = data_kin_vs[gotvz,10]
            vzerr2 = data_kin_vs[gotvz,11]
            ms2 = data_kin_vs[gotvz,19]
            Nbin2 = np.sqrt(np.sum(ms2))
        else:
            f = open(data_file_kin_vs,'r')
            data_kin_vs = np.genfromtxt(f)
            f.close()
            gotvz = np.where(data_kin_vs[:,10])[0]
            R1 = data_kin_vs[gotvz,4]*dgal/arcmin
            vz1 = data_kin_vs[gotvz,10]
            vzerr1 = data_kin_vs[gotvz,11]
            R2 = data_kin_vs[gotvz,4]*dgal/arcmin
            vz2 = data_kin_vs[gotvz,10]
            vzerr2 = data_kin_vs[gotvz,11]
            ms1 = data_kin_vs[gotvz,19]
            ms2 = data_kin_vs[gotvz,24]
            Nbin1 = np.sqrt(np.sum(ms1))
            Nbin2 = np.sqrt(np.sum(ms2))
            print 'Nbin1 (kinematic):', Nbin1
            print 'Nbin2 (kinematic):', Nbin2

        nmonte = 1000
        Nbin = Nbin1
        rbin1_kin, sigpmean1, sigperr1, \
            vs1bin1, vs2bin1, vs1err1, vs2err1 = \
            calc_virial_moments(Rhalf1,nmonte,R1,vz1,vzerr1,ms1,\
                                pfits1)
        Nbin = Nbin2
        rbin2_kin, sigpmean2, sigperr2, \
            vs1bin2, vs2bin2, vs1err2, vs2err2 = \
            calc_virial_moments(Rhalf2,nmonte,R2,vz2,vzerr2,ms2,\
                                pfits2)

        #For plotting:
        rbin = rbin1_phot
        rbin1 = rbin1_phot
        rbin2 = rbin2_phot
        Rhalf = Rhalf1

    #Add in Jardel inner disperison data for Draco:
    if (include_jardel == 'yes'):
        print 'Including Jardel inner dispersion point ... '
        f = open(jardel_file,'r')
        data_jar = np.genfromtxt(f)
        f.close()
        RA_DRA_deg = 260.0516667  #McConnachie rev.
        DEC_DRA_deg = 57.9152778  #
        vsys_DRA_walker = -291.0  #

        print 'Systemtic vel of virus data (sanity check):',\
            np.sum(data_jar[:,2])/np.float(len(data_jar[:,2])),\
            vsys_DRA_walker

        rjar_arr = np.sqrt(((data_jar[:,0] - RA_DRA_deg)*\
                        np.cos(DEC_DRA_deg/360.0*2.0*np.pi))**2.0+\
                       (data_jar[:,1] - DEC_DRA_deg)**2.0)/\
                       360.0*2.0*np.pi*dgal

        rjar = np.sqrt(np.sum(rjar_arr**2.0)/np.float(len(rjar_arr))-\
                      (np.sum(rjar_arr)/np.float(len(rjar_arr)))**2.0)

        #And Monte-Carlo to get sigpjar and sigperrjar
        nmontejar = 1000
        vtemp = np.zeros(len(data_jar[:,2]))
        vpureerr = np.zeros(len(data_jar[:,2]))
        sigpjar_store = np.zeros(nmontejar)
        sigp_pureerr = np.zeros(nmontejar)
        vsys_store = np.zeros(nmontejar)
        for i in range(nmontejar):
            #Draw velocity from Gaussian error:
            vtemp = data_jar[:,2]+\
                np.random.normal(0.0,data_jar[:,3])
            vsys = np.sum(vtemp)/np.float(len(vtemp))
            vtemp = vtemp - vsys
            vpureerr = np.random.normal(0.0,data_jar[:,3])
            sigpjar_store[i] = \
                np.sum(vtemp**2.0)/np.float(len(vtemp))-\
                (np.sum(vtemp)/np.float(len(vtemp)))**2.0
            sigp_pureerr[i] = \
                np.sum(vpureerr**2.0)/np.float(len(vpureerr))-\
                (np.sum(vpureerr)/np.float(len(vpureerr)))**2.0
            vsys_store[i] = vsys

        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = calcmedquartnine(vsys_store)
        vsys = median
        vsys_err = (sixhigh-sixlow)/2.0
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = calcmedquartnine(sigpjar_store)
        sigpjar2 = median
        sigperr_meas = (sixhigh-sixlow)/2.0
        median, sixlow, sixhigh, ninelow, ninehigh,\
            nineninehigh, nineninelow = calcmedquartnine(sigp_pureerr)
        sigp_pe2 = median
        sigp_pe2err = (sixhigh-sixlow)/2.0

        sigpjar_err2 = np.sqrt(sigperr_meas**2.0 + sigpjar2**2.0/\
                               np.float(len(vtemp)))
        sigpjar2 = sigpjar2 - sigp_pe2
        
        sigpjar = np.sqrt(sigpjar2)
        sigperrjar = sigpjar_err2/2.0/sigpjar

        print 'Systematic velocity recovered:', vsys, vsys_err

        rbin_kin = np.concatenate((np.array([rjar]),rbin_kin))
        sigpmean = np.concatenate((np.array([sigpjar]),sigpmean))
        sigperr = np.concatenate((np.array([sigperrjar]),sigperr))
        
elif (data_file_type == 'LeoT'):
    f = open(data_file_kin,'r')
    data_kin = np.genfromtxt(f)
    f.close()
    f = open(data_file_phot,'r')
    data_phot = np.genfromtxt(f)
    f.close()
    
    R = data_kin[:,6]/1000.0
    vz = data_kin[:,1]
    vzerr = data_kin[:,2]
    ms = data_kin[:,7]

    ecut = 0.3
    Ruset = R[vzerr/vz < ecut]
    vzuset = vz[vzerr/vz < ecut]
    vzuset = vzuset - np.sum(vzuset)/np.float(len(vzuset))
    vzerruset = vzerr[vzerr/vz < ecut]
    msuset = ms[vzerr/vz < ecut]
    
    ecut = 0.3
    Ruset = R[vzerr/vz < ecut]
    vzuset = vz[vzerr/vz < ecut]
    vzuset = vzuset - np.sum(vzuset)/np.float(len(vzuset))
    vzerruset = vzerr[vzerr/vz < ecut]
    msuset = ms[vzerr/vz < ecut]
    
    print 'Initial numbers of stars:', np.sum(ms)
    print 'Stars left after ecut:',np.sum(msuset)
    
    pcut = 0.6
    Ruse = Ruset[msuset > pcut]
    vzuse = vzuset[msuset > pcut]
    vzuse = vzuse - np.sum(vzuse)/np.float(len(vzuse))
    vzerruse = vzerruset[msuset > pcut]
    msuse = msuset[msuset > pcut]

    print 'Stars left after pcut:',len(msuse)
    #Upsample the surface density profile for fitting and cut
    #on the maxdatfitrad:
    rbin_phot = data_phot[:,0]*dgal/arcmin
    surfden_t = data_phot[:,1]
    norm = integrator(surfden_t*2.0*np.pi*rbin_phot,rbin_phot)
    surfden = surfden_t / norm
    nstars_per_bin = data_phot[:,1]*2.0*np.pi*data_phot[:,0]
    surfdenerr = surfden / np.sqrt(nstars_per_bin)

    rbin_photfit = np.logspace(np.log10(np.min(rbin_phot)),\
                               np.log10(np.max(\
                               rbin_phot[rbin_phot < maxdatfitrad])),100)
    surfdenfit = np.interp(rbin_photfit,rbin_phot[rbin_phot < maxdatfitrad],\
                           surfden[rbin_phot < maxdatfitrad])
    surfdenerrfit = np.interp(rbin_photfit,\
                              rbin_phot[rbin_phot < maxdatfitrad],\
                              surfdenerr[rbin_phot < maxdatfitrad])
    pfits = tracerfit(p0in,p0in_min,p0in_max,\
                      rbin_phot,surfden,surfdenerr)
    Mstar_rad = rbin_phot
    norm = np.max(threeplummass(\
                  np.linspace(0,Mstar_rlim,100),pfits[0],pfits[1],\
                                pfits[2],\
                                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/\
                               norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])
    Mstar_surf = threeplumsurf(ranal,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])
    Rhalf = Rhalf_lum
                
    #Virial moments: 
    if (Nbinin < 0):
        Nbin = np.sqrt(np.sum(ms))
    else:
        Nbin = Nbinin
    print 'Nbin_kin = ', Nbin
    nmonte = 1000
    rbin_kin, sigpmean, sigperr, \
        vs1bin, vs2bin, vs1err, vs2err = \
                calc_virial_moments(Rhalf,nmonte,Ruse,vzuse,\
                                    vzerruse,msuse,\
                                    pfits)
    
    #For plotting:
    rbin = rbin_phot

elif (data_file_type == 'Ocen'):
    f = open(data_file_kin,'r')
    data_kin = np.genfromtxt(f)
    f.close()
    f = open(data_file_phot,'r')
    data_phot = np.genfromtxt(f)
    f.close()
    f = open(data_file_pmr,'r')
    data_pmr = np.genfromtxt(f)
    f.close()
    f = open(data_file_pmt,'r')
    data_pmt = np.genfromtxt(f)
    f.close()
    
    #Set up data:
    rbin_phot = data_phot[:,0]*dgal/arcsec
    surfden = data_phot[:,1]/np.max(data_phot[:,1])
    surfdenerr =  (data_phot[:,2]/np.max(data_phot[:,1])+\
                   data_phot[:,3]/np.max(data_phot[:,1]))/2.0
    rbin_kin = data_kin[:,0]*dgal/arcsec
    sigpmean = data_kin[:,1]
    sigperr = (-data_kin[:,2]+data_kin[:,3])/2.0
    sigpmr = np.interp(rbin_kin,data_pmr[:,0]*dgal/arcsec,\
                       data_pmr[:,1]*dgal/arcsec*1.0e-3/year*kpc/kms,left=0)
    sigpmrerr = np.interp(rbin_kin,data_pmr[:,0]*dgal/arcsec,\
                          (-data_pmr[:,2]+data_pmr[:,3])*\
                        dgal/arcsec*1.0e-3/year*kpc/kms,\
                          left=1e10)
    sigpmt = np.interp(rbin_kin,data_pmt[:,0]*dgal/arcsec,\
                       data_pmt[:,1]*dgal/arcsec*1.0e-3/year*kpc/kms,left=0)
    sigpmterr = np.interp(rbin_kin,data_pmt[:,0]*dgal/arcsec,\
                          (-data_pmt[:,2]+data_pmt[:,3])*\
                          dgal/arcsec*1.0e-3/year*kpc/kms,left=1e10)

    #Sort out the stellar mass profile:
    rbin_photfit = np.logspace(np.log10(np.min(rbin_phot)),\
                            np.log10(np.max(rbin_phot[rbin_phot < maxdatfitrad])),100)
    surfdenfit = np.interp(rbin_photfit,rbin_phot[rbin_phot < maxdatfitrad],\
                        surfden[rbin_phot < maxdatfitrad])
    surfdenerrfit = np.interp(rbin_photfit,\
                              rbin_phot[rbin_phot < maxdatfitrad],\
                            surfdenerr[rbin_phot < maxdatfitrad])
    pfits = tracerfit(p0in,p0in_min,p0in_max,\
                      rbin_photfit,surfdenfit,surfdenerrfit)
    
    Mstar_rad = np.logspace(np.log10(rbin_inner),\
                        np.log10(rbin_outer),100)                            
    norm = np.max(threeplummass(\
                                np.linspace(0,Mstar_rlim,100),\
                                pfits[0],pfits[1],\
                                pfits[2],\
                                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/\
                               norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])
    Mstar_surf = threeplumsurf(ranal,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])

    rbin = rbin_phot
    Rhalf = Rhalf_lum
        
elif (data_file_type == 'SegI'):
    f = open(data_file_kin_vs,'r')
    data_kin = np.genfromtxt(f)
    f.close()
    f = open(data_file_phot,'r')
    data_phot = np.genfromtxt(f)
    f.close()
                        
    R = data_kin[:,3]*dgal/arcmin
    vz = data_kin[:,1]
    vzerr = data_kin[:,2]
    ms = data_kin[:,6]
    
    ecut = 0.3
    Ruset = R[vzerr/vz < ecut]
    vzuset = vz[vzerr/vz < ecut]
    vzuset = vzuset - np.sum(vzuset)/np.float(len(vzuset))
    vzerruset = vzerr[vzerr/vz < ecut]
    msuset = ms[vzerr/vz < ecut]
    
    print 'Initial numbers of stars:', len(ms)
    print 'Stars left after ecut:',len(msuset)

    if (prob_cut == 'p09'):
        pcut = 0.9
    else:
        pcut = 0.95
    Ruse = Ruset[msuset > pcut]
    vzuse = vzuset[msuset > pcut]
    vzuse = vzuse - np.sum(vzuse)/np.float(len(vzuse))
    vzerruse = vzerruset[msuset > pcut]
    msuse = msuset[msuset > pcut]
    
    print 'Stars left after pcut:',len(msuse)

    rbin_phot = data_phot[:,0]*dgal/arcmin
    surfden = data_phot[:,1]
    surfdenerr = (-data_phot[:,2]+data_phot[:,3])/2.0
                
    #Upsample the surface density profile for fitting and cut
    #on the maxdatfitrad:
    Rhalf = Rhalf_lum
    rbin_photfit = np.logspace(np.log10(np.min(rbin_phot)),\
                np.log10(np.max(rbin_phot[rbin_phot < maxdatfitrad])),100)
    surfdenfit = np.interp(rbin_photfit,rbin_phot[rbin_phot < maxdatfitrad],\
                           surfden[rbin_phot < maxdatfitrad])
    surfdenerrfit = np.interp(rbin_photfit,\
                              rbin_phot[rbin_phot < maxdatfitrad],\
                              surfdenerr[rbin_phot < maxdatfitrad])
    pfits = tracerfit(p0in,p0in_min,p0in_max,\
                      rbin_photfit,surfdenfit,surfdenerrfit)

    Mstar_rad = rbin_phot
    norm = np.max(threeplummass(\
                                np.linspace(0,Mstar_rlim,100),\
                                pfits[0],pfits[1],\
                                pfits[2],\
                                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/\
                               norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])
    Mstar_surf = threeplumsurf(ranal,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])
    
    if (Nbinin < 0):
        Nbin = np.sqrt(np.sum(ms))
    else:
        Nbin = Nbinin
    print 'Nbin_kin = ', Nbin
    nmonte = 1000
    rbin_kin, sigpmean, sigperr, \
        vs1bin, vs2bin, vs1err, vs2err = \
            calc_virial_moments(Rhalf,nmonte,Ruse,vzuse,\
                                vzerruse,msuse,\
                                pfits)
    
    #For plotting:
    rbin = rbin_phot
    
elif (data_file_type == 'Tucana'):
    f = open(data_file_kin,'r')
    data_kin = np.genfromtxt(f)
    f.close()

    R = data_kin[:,8]
    vz = data_kin[:,1]
    vzerr = data_kin[:,2]
    ms = data_kin[:,17]
    vz = vz - np.sum(ms*vz)/np.sum(ms)

    if (Nbinin < 0):
        Nbin = np.sqrt(np.sum(ms))
    else:
        Nbin = Nbinin_phot
    print 'Number of kinematic data points (effective):', np.sum(ms)
    
    rbin_phot, surfden, surfdenerr, \
        rbin_photfit, surfdenfit, surfdenerrfit, \
        Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits = \
        tucana_surf_sort(data_file_phot,Rhalf_lum)

    if (addinfake == 'yes'):
        f = open('./fake'+fakehigh+'.txt','r')
        data = np.genfromtxt(f)
        R = data[:,0]
        vz = data[:,1]
        vzerr = data[:,2]
        ms = data[:,3]
    
    if (Nbinin > 0):
        Nbin = Nbinin
    nmonte = 1000
    rbin_kin, sigpmean, sigperr, \
        vs1bin, vs2bin, vs1err, vs2err = \
        calc_virial_moments(Rhalf,nmonte,R,vz,vzerr,ms,\
                            pfits)
    rbin = rbin_phot 
elif (data_file_type == 'Aquarius'):
    f = open(data_file_kin,'r')
    data_kin = np.genfromtxt(f)
    f.close()
    f = open(data_file_phot,'r')
    data_phot = np.genfromtxt(f)
    f.close()
    
    R = data_kin[:,0]
    vz = data_kin[:,1]
    vzerr = data_kin[:,2]
    ms = np.zeros(len(R))+1.0
    Nbin = 6
    print 'Aquarius Nbin = ', Nbin

    rbin_phot = data_phot[:,0]
    surfden = data_phot[:,1]
    surfdenerr = data_phot[:,2]

    pfits = tracerfit(p0in,p0in_min,p0in_max,\
                      rbin_phot,surfden,surfdenerr)
    Mstar_rad = rbin_phot
    norm = np.max(threeplummass(\
            np.linspace(0,Mstar_rlim,100),pfits[0],pfits[1],\
                pfits[2],\
                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])
    Mstar_surf = threeplumsurf(ranal,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])
    Rhalf = Rhalf_lum

    #Virial moments:
    nmonte = 1000
    rbin_kin, sigpmean, sigperr, \
        vs1bin, vs2bin, vs1err, vs2err = \
        calc_virial_moments(Rhalf,nmonte,R,vz,vzerr,ms,\
                                pfits)

    #For plotting:
    rbin = rbin_phot

elif (data_file_type == 'sluggs'):
    f = open(data_file_phot1,'r')
    data_phot1 = np.genfromtxt(f)
    f.close()
    f = open(data_file_phot2,'r')
    data_phot2 = np.genfromtxt(f)
    f.close()
    f = open(data_file_phot,'r')
    data_phot = np.genfromtxt(f)
    f.close()
    
    rbin_phot = data_phot[:,0]
    surfden = data_phot[:,1]
    surfdenerr = data_phot[:,2]
    rbin1_phot = data_phot1[:,0]
    surfden1 = data_phot1[:,1]
    surfdenerr1 = data_phot1[:,2]
    rbin2_phot = data_phot2[:,0]
    surfden2 = data_phot2[:,1]
    surfdenerr2 = data_phot2[:,2]

    pfits1 = tracerfit(p0in,p0in_min,p0in_max,\
                       rbin1_phot,surfden1,surfdenerr1)
    norm = np.max(threeplummass(\
            np.linspace(0,Mstar_rlim,100),pfits1[0],pfits1[1],\
                pfits1[2],\
                pfits1[3],pfits1[4],pfits1[5]))
    Mstar_surf = threeplumsurf(ranal,pfits1[0]/norm,pfits1[1]/norm,\
                                   pfits1[2]/norm,\
                                   pfits1[3],pfits1[4],pfits1[5])
    if (Rhalf_lum1 < 0):
        Mcum_surf = 0.0
        i = 1
        while (Mcum_surf < (pfits1[0]+pfits1[1]+pfits1[2])/2.0/norm):
            Mcum_surf = \
                Mcum_surf + \
                2.0*np.pi*ranal[i]*Mstar_surf[i]*(ranal[i]-ranal[i-1])
            i = i + 1
        Rhalf1 = ranal[i-1]
        print 'Rhalf calculated: ', Rhalf1
    else:
        Rhalf1 = Rhalf_lum1

    pfits2 = tracerfit(p0in,p0in_min,p0in_max,\
                       rbin2_phot,surfden2,surfdenerr2)
    norm = np.max(threeplummass(\
            np.linspace(0,Mstar_rlim,100),pfits2[0],pfits2[1],\
                pfits2[2],\
                pfits2[3],pfits2[4],pfits2[5]))
    Mstar_surf = threeplumsurf(ranal,pfits2[0]/norm,pfits2[1]/norm,\
                                   pfits2[2]/norm,\
                                   pfits2[3],pfits2[4],pfits2[5])
    if (Rhalf_lum2 < 0):
        Mcum_surf = 0.0
        i = 1
        while (Mcum_surf < (pfits2[0]+pfits2[1]+pfits2[2])/2.0/norm):
            Mcum_surf = \
                Mcum_surf + \
                2.0*np.pi*ranal[i]*Mstar_surf[i]*(ranal[i]-ranal[i-1])
            i = i + 1
        Rhalf2 = ranal[i-1]
        print 'Rhalf calculated: ', Rhalf2
    else:
        Rhalf2 = Rhalf_lum2

    pfits = tracerfit(p0starin,p0starin_min,p0starin_max,\
                       rbin_phot,surfden,surfdenerr)
    Mstar_rad = rbin_phot
    norm = np.max(threeplummass(\
            np.linspace(0,Mstar_rlim,100),pfits[0],pfits[1],\
                pfits[2],\
                pfits[3],pfits[4],pfits[5]))
    Mstar_prof = threeplummass(Mstar_rad,pfits[0]/norm,pfits[1]/norm,\
                               pfits[2]/norm,\
                               pfits[3],pfits[4],pfits[5])

    #TEST:
    if (testmode == 'yes'):
        Mstar_surf = threeplumsurf(ranal,pfits[0],pfits[1],\
                                   pfits[2],\
                                   pfits[3],pfits[4],pfits[5])
        plt.figure()
        plt.loglog()
        plt.errorbar(rbin_phot,surfden,yerr=surfdenerr,fmt='o')
        plt.plot(ranal,Mstar_surf)
        plt.show()

    f = open(data_file_kin_vs1,'r')
    data_kin1 = np.genfromtxt(f)
    f.close()
    f = open(data_file_kin_vs2,'r')
    data_kin2 = np.genfromtxt(f)
    f.close()
    R1 = data_kin1[:,0]
    vz1 = data_kin1[:,1]
    vzerr1 = data_kin1[:,2]
    ms1 = np.zeros(len(R1)) + 1.0
    R2 = data_kin2[:,0]
    vz2 = data_kin2[:,1]
    vzerr2 = data_kin2[:,2]
    ms2 = np.zeros(len(R2)) + 1.0
    Nbin1 = np.sqrt(np.sum(ms1))
    Nbin2 = np.sqrt(np.sum(ms2))
    print 'Nbin1:', Nbin1
    print 'Nbin2:', Nbin2

    nmonte = 1000
    Nbin = Nbin1
    rbin1_kin, sigpmean1, sigperr1, \
        vs1bin1, vs2bin1, vs1err1, vs2err1 = \
        calc_virial_moments(Rhalf1,nmonte,R1,vz1,vzerr1,ms1,\
                            pfits1)
    Nbin = Nbin2
    rbin2_kin, sigpmean2, sigperr2, \
        vs1bin2, vs2bin2, vs1err2, vs2err2 = \
        calc_virial_moments(Rhalf2,nmonte,R2,vz2,vzerr2,ms2,\
                            pfits2)
    
    #For plotting:
    rbin = rbin1_phot
    rbin1 = rbin1_phot
    rbin2 = rbin2_phot
    Rhalf = Rhalf1
elif (data_file_type == 'genina'):
    data = h5py.File(data_file, 'r')
    pos_kin = data['KinematicsPositions'].value
    R_kin = np.sqrt(pos_kin[:,0]**2.0+pos_kin[:,1]**2.0)/1000.0
    vz_kin = data['KinematicsVelocities'].value
    vzerr_kin = np.zeros(len(vz_kin))+2.0
    ms_kin = data['KinematicsMasses'].value
    ms_kin = ms_kin / np.sum(ms_kin)*len(ms_kin)

    #Subtract space velocity:
    vz_kin = vz_kin - np.sum(ms_kin*vz_kin)/np.sum(ms_kin)
    pos_phot = data['PhotometryPositions'].value
    R_phot = np.sqrt(pos_phot[:,0]**2.0+pos_phot[:,1]**2.0)/1000.0
    ms_phot = data['PhotometryMasses'].value
    ms_phot = ms_phot / np.sum(ms_phot)*len(ms_phot)
    Mstar = data['StellarMass3R'].value
    Mstar_err = Mstar*0.25
    
    if (Nbinin_phot < 0):
        Nbin = np.sqrt(np.sum(ms_phot))
    else:
        Nbin = Nbinin_phot
    print 'Number of photometric data points:', len(ms_phot)
    print 'Number of stars per photometric bin:', Nbin                                

    #Calculate the surface density:
    rbin_phot, surfden, surfdenerr, \
        rbin_photfit, surfdenfit, surfdenerrfit, \
        Mstar_rad, Mstar_prof, Mstar_surf, Rhalf, pfits = \
        genina_surf_sort(R_phot,ms_phot)

    #And calculate the velocity disperison profile, vs1 and vs2:
    if (Nbinin_kin < 0):
        Nbin = np.sqrt(np.sum(ms_kin))
    else:
        Nbin = Nbinin_kin
    print 'Number of kinematic data points:', len(ms_kin)
    print 'Number of stars per kinematic bin:', Nbin
    nmonte = 1000
    rbin_kin, sigpmean, sigperr, \
        vs1bin, vs2bin, vs1err, vs2err = \
        calc_virial_moments(Rhalf,nmonte,R_kin,vz_kin,vzerr_kin,ms_kin,\
                            pfits)

    #For calc_virial_moments which is run later for
    #calc. VSPs (not elegant, but needed to be compatible with
    #older data reading parts of the code)
    R = R_kin
    vz = vz_kin
    vzerr = vzerr_kin
    ms = ms_kin
        
    #For plotting:
    rbin = rbin_phot
    
#Calculate the model bins from the data:
if (singlesplit == 'single'):
    mbinpars = [0.25*Rhalf,0.5*Rhalf,Rhalf,2.0*Rhalf,4.0*Rhalf]
    if (testmorehighbins == 'yes'):
        mbinpars = [0.25*Rhalf,0.5*Rhalf,Rhalf,2.0*Rhalf,\
                    4.0*Rhalf,8.0*Rhalf]
    if (singleinbin == 'yes'):
        mbinpars = [Rhalf,2.0*Rhalf,4.0*Rhalf]
elif (singlesplit == 'split'):
    Rhalfarr = [Rhalf1,Rhalf2]
    Rhalfsmall = np.min(Rhalfarr)
    Rhalfbig = np.max(Rhalfarr)
    Rhalfmid = (Rhalfbig-Rhalfsmall)/2.0+Rhalfsmall
    mbinpars = [0.25*Rhalfsmall,0.5*Rhalfsmall,Rhalfsmall,\
                Rhalfbig,2.0*Rhalfbig,4.0*Rhalfbig]
    if (testmorehighbins == 'yes'):
        mbinpars = [0.25*Rhalfsmall,0.5*Rhalfsmall,Rhalfsmall,\
                    Rhalfbig,2.0*Rhalfbig,4.0*Rhalfbig,8.0*Rhalfbig]

#Set beta scale radius based also on Rhalf:
if (betprior == 'poszero' or betprior == 'default2' or betprior == 'zero'):
    if (singlesplit == 'single'):
        betr0min = np.log10(0.5*Rhalf)
        betr0max = np.log10(2.0*Rhalf)
    else:
        betr0min = np.log10(0.5*Rhalfsmall)
        betr0max = np.log10(2.0*Rhalfbig)
        
#Set r-grid based also on Rhalf:
if (rmax < 0):
    if (singlesplit == 'single'):
        rmin = Rhalf / 100.0
        rmax = Rhalf * 50.0
    else:
        rmin = Rhalfsmall / 100.0
        rmax = Rhalfbig * 50.0
print 'Inner/outer radial grid:', rmin, rmax
print 'Inner/outer radial bin:', np.min(rbin),np.max(rbin)

M = lambda r, Mpars: \
    corenfw_tides_mass(r,Mpars[0],Mpars[1],Mpars[2],\
                         Mpars[3],Mpars[4],Mpars[5])
rho = lambda r, Mpars: \
    corenfw_tides_den(r,Mpars[0],Mpars[1],Mpars[2],\
                        Mpars[3],Mpars[4],Mpars[5])
dlnrhodlnr = lambda r, Mpars: \
    corenfw_tides_dlnrhodlnr(r,Mpars[0],Mpars[1],Mpars[2],\
                               Mpars[3],Mpars[4],Mpars[5])
n_mpars = 6
if (SIDM_Draco_priors == 'yes'):
    print 'Using SIDM Draco priors:'
#    logM200low = 8.75
#    logM200high = 10.25

    logM200low = 9.0
    logM200high = 11.0
    
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)
    logchigh = np.log10(chigh)+0.1
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    use_rclinear = 'no'
    use_rtlinear = 'no'
    rclow = 1e-2
#    rchigh = 10.0**0.5
    rchigh = 500.0
    if (use_rclinear == 'yes'):
        logrclow = rclow
        logrchigh = rchigh
    else:
        logrclow = np.log10(rclow)
        logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 1.0
    nlow = 0.99
    nhigh = 1.0
    rtlow = 0.5*Rhalf
    rthigh = 10.0*Rhalf
    if (use_rtlinear == 'yes'):
        logrtlow = rtlow
        logrthigh = rthigh
    else:
        logrtlow = np.log10(rtlow)
        logrthigh = np.log10(rthigh)
    dellow = 2.5
    delhigh = 5.0
elif (SIDM_Draco_priors == 'Zavala'):
    print 'Using Zavala priors:'
    logM200low = 8.75
    logM200high = 10.25
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)*1.4
    logchigh = np.log10(chigh)+0.2
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    use_rclinear = 'no'
    use_rtlinear = 'no'
    rclow = 1e-2
    rchigh = 10.0**0.5
    if (use_rclinear == 'yes'):
        logrclow = rclow
        logrchigh = rchigh
    else:
        logrclow = np.log10(rclow)
        logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 1.0
    nlow = 0.99
    nhigh = 1.0
    rtlow = 0.5*Rhalf
    rthigh = 10.0*Rhalf
    if (use_rtlinear == 'yes'):
        logrtlow = rtlow
        logrthigh = rthigh
    else:
        logrtlow = np.log10(rtlow)
        logrthigh = np.log10(rthigh)
    dellow = 2.5
    delhigh = 5.0
elif (SIDM_Draco_priors == 'nvary'):
    print 'Using more generous SIDM priors + n-vary!'
    logM200low = 8.5
    logM200high = 10.5
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)*1.4
    logchigh = np.log10(chigh)+0.2
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    use_rclinear = 'no'
    use_rtlinear = 'no'
    rclow = 1e-2
    rchigh = 10.0**0.5
    if (use_rclinear == 'yes'):
        logrclow = rclow
        logrchigh = rchigh
    else:
        logrclow = np.log10(rclow)
        logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 5.0
    nlow = 0.0
    nhigh = 1.0
    rtlow = Rhalf
    rthigh = 10.0*Rhalf
    if (use_rtlinear == 'yes'):
        logrtlow = rtlow
        logrthigh = rthigh
    else:
        logrtlow = np.log10(rtlow)
        logrthigh = np.log10(rthigh)
    dellow = 3.01
    delhigh = 5.0
elif (SIDM_Draco_priors == 'cosmo_c'):
    print 'Using n-vary priors + cosmo_c:'
    use_rclinear = 'no'
    use_rtlinear = 'no'
    logM200low = 8.5
    logM200high = 10.5
    clow = 1.0
    chigh = 100.0
    print 'clow, chigh:', clow, chigh
    rclow = 1e-2
    rchigh = 10.0**0.5
    logrclow = np.log10(rclow)
    logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 5.0
    nlow = 0.0
    nhigh = 1.0
    rtlow = Rhalf
    rthigh = 10.0*Rhalf
    logrtlow = np.log10(rtlow)
    logrthigh = np.log10(rthigh)
    dellow = 3.01
    delhigh = 5.0
        
elif (SIDM_Draco_priors == 'Ocen'):
    print 'Using Ocen priors:'
    logM200low = 1.0
    logM200high = 10.0
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.2
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)
    logchigh = np.log10(chigh)+0.2
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    use_rclinear = 'no'
    use_rtlinear = 'no'
    rclow = 1e-2
    rchigh = 500.0
    if (use_rclinear == 'yes'):
        logrclow = rclow
        logrchigh = rchigh
    else:
        logrclow = np.log10(rclow)
        logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 50.0
    nlow = 0.0
    nhigh = 1.0
    rtlow = Rhalf
    rthigh = 10.0*Rhalf
    if (use_rtlinear == 'yes'):
        logrtlow = rtlow
        logrthigh = rthigh
    else:
        logrtlow = np.log10(rtlow)
        logrthigh = np.log10(rthigh)
    dellow = 3.5
    delhigh = 10.0    
elif (SIDM_Draco_priors == 'Tucana'):
    print 'Using Tucana priors:'
    logM200low = 9.0
    logM200high = 11.0
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)*1.4
    logchigh = np.log10(chigh)+0.2
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    rclow = 1e-2
    rchigh = 10.0**0.5
    logrclow = np.log10(rclow)
    logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 5.0
    nlow = 0.0
    nhigh = 1.0
    rtlow = Rhalf
    rthigh = 10.0*Rhalf
    logrtlow = np.log10(rtlow)
    logrthigh = np.log10(rthigh)
    dellow = 3.5
    delhigh = 5.0
elif (SIDM_Draco_priors == 'Genina'):
    print 'Using Genina priors:'
    logM200low = 7.0
    logM200high = 11.0
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)*1.4
    logchigh = np.log10(chigh)+0.2
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    rclow = 1e-2
    rchigh = 10.0**0.5
    logrclow = np.log10(rclow)
    logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 5.0
    nlow = 0.99
    nhigh = 1.0
    rtlow = Rhalf
    rthigh = 10.0*Rhalf
    logrtlow = np.log10(rtlow)
    logrthigh = np.log10(rthigh)
    dellow = 3.5
    delhigh = 5.0
elif (SIDM_Draco_priors == 'Genina2'):
    print 'Using Genina2 priors:'
    logM200low = 7.5
    logM200high = 11.5
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)*1.4
    logchigh = np.log10(chigh)+0.2
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    rclow = 1e-2
    rchigh = 500.0
    logrclow = np.log10(rclow)
    logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 5.0
    nlow = 0.99
    nhigh = 1.0
    rtlow = Rhalf
    rthigh = 10.0*Rhalf
    logrtlow = np.log10(rtlow)
    logrthigh = np.log10(rthigh)
    dellow = 3.5
    delhigh = 5.0
elif (SIDM_Draco_priors == 'abund'):
    print 'Using abundance matching priors:'
    f = open('./Data/mstar-mhalo-field.txt','r')
    data = np.genfromtxt(f)
    f.close()
    M200lowarr = data[:,0][np.where(data[:,1])]
    Mstarlow = data[:,1][np.where(data[:,1])]
    M200higharr = data[:,0][np.where(data[:,2])]
    Mstarhigh = data[:,2][np.where(data[:,2])]
    M200high = np.interp(Mstar+Mstar_err,Mstarlow,M200lowarr)
    M200low = np.interp(Mstar-Mstar_err,Mstarhigh,M200higharr) 
    if (whichdata == 'Fornax'):
        #Overwrite here to correct for quenching:
        M200low = (21.9-7.4*2.0)*1.0e9
        M200high = (21.9+7.4*2.0)*1.0e9
        if (Mstar == 2e7):
            M200low = (13.1-4.84)*1.0e9
            M200high = (13.1+4.84)*1.0e9
    elif (whichdata == 'Draco'):
        M200low = (1.8-0.7*2.0)*1.0e9
        M200high = (1.8+0.7*2.0)*1.0e9
    logM200low = np.log10(M200low)
    logM200high = np.log10(M200high)        
    print 'logM200low, logM200high:', logM200low, logM200high
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1    
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)*1.4
    logchigh = np.log10(chigh)+0.1
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    rclow = 1e-2
    rchigh = 500.0
    logrclow = np.log10(rclow)
    logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 5.0
    nlow = 0.99
    nhigh = 1.0
    rtlow = 3.0*Rhalf
    rthigh = 10.0*Rhalf
    logrtlow = np.log10(rtlow)
    logrthigh = np.log10(rthigh)
    dellow = 2.0
    delhigh = 5.0
elif (SIDM_Draco_priors == 'abundvary'):
    print 'Using abudnace matching priors + n-vary!'
    if (whichdata == 'UMi'):
        M200low = (2.8-1.1*2.0)*1.0e9
        M200high = (2.8+1.1*2.0)*1.0e9
    elif (whichdata == 'Draco'):
        M200low = (1.8-0.7*2.0)*1.0e9
        M200high = (1.8+0.7*2.0)*1.0e9
    elif (whichdata == 'Sculptor'):
        M200low = (5.7-2.3*2.0)*1.0e9
        M200high = (5.7+2.3*2.0)*1.0e9
    elif (whichdata == 'Sextans'):
        M200low = (2.0-0.8*2.0)*1.0e9
        M200high = (2.0+0.8*2.0)*1.0e9
    elif (whichdata == 'LeoI'):
        M200low = (5.6-2.2*2.0)*1.0e9
        M200high = (5.6+2.2*2.0)*1.0e9
    elif (whichdata == 'LeoII'):
        M200low = (1.6-0.7*2.0)*1.0e9
        M200high = (1.6+0.7*2.0)*1.0e9
    elif (whichdata == 'Carina'):
        M200low = (0.8-0.3*2.0)*1.0e9
        M200high = (0.8+0.3*2.0)*1.0e9
    elif (whichdata == 'Fornax'):
        M200low = (21.9-7.4*2.0)*1.0e9
        M200high = (21.9+7.4*2.0)*1.0e9
        if (Mstar == 2e7):
            M200low = (13.1-4.84)*1.0e9
            M200high = (13.1+4.84)*1.0e9
    elif (whichdata == 'SegI'):
        M200low = (1.0-0.4*2.0)*1.0e9
        M200high = (1.0+0.4*2.0)*1.0e9
        
    logM200low = np.log10(M200low)
    logM200high = np.log10(M200high)
        
    print 'logM200low, logM200high:', logM200low, logM200high
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)*1.4
    logchigh = np.log10(chigh)+0.1
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh

    rclow = 1e-2
    rchigh = 500.0
    logrclow = np.log10(rclow)
    logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 5.0
    nlow = 0.0
    nhigh = 1.0
    rtlow = 2.0*Rhalf
    rthigh = 25.0*Rhalf
    logrtlow = np.log10(rtlow)
    logrthigh = np.log10(rthigh)
    dellow = 3.0
    delhigh = 5.0
elif (SIDM_Draco_priors == 'NFW'):
    print 'Using NFW priors'
    logM200low = 8.0
    logM200high = 11.0
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)*1.4
    logchigh = np.log10(chigh)+0.2
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    rclow = 1e-2
    rchigh = 1e-1
    logrclow = np.log10(rclow)
    logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 5.0
    nlow = 1e-3
    nhigh = 1e-2
    rtlow = 5.0*Rhalf
    rthigh = 10.0*Rhalf
    logrtlow = np.log10(rtlow)
    logrthigh = np.log10(rthigh)
    dellow = 3.0
    delhigh = 4.0
else:
    print 'Using more generous SIDM priors:'
    logM200low = 8.5
    logM200high = 10.5
    clow = cosmo_cfunc(10.0**logM200high,h)
    logclow = np.log10(clow)-0.1
    clow = 10.0**logclow
    chigh = cosmo_cfunc(10.0**logM200low,h)
    logchigh = np.log10(chigh)+0.1
    chigh = 10.0**logchigh
    print 'clow, chigh:', clow, chigh
    use_rclinear = 'no'
    use_rtlinear = 'no'
    rclow = 1e-2
    rchigh = 10.0**0.5
    if (use_rclinear == 'yes'):
        logrclow = rclow
        logrchigh = rchigh
    else:
        logrclow = np.log10(rclow)
        logrchigh = np.log10(rchigh)
    sigmlow = 1e-3
    sigmhigh = 5.0
    nlow = 0.99
    nhigh = 1.0
    rtlow = Rhalf
    rthigh = 10.0*Rhalf
    if (use_rtlinear == 'yes'):
        logrtlow = rtlow
        logrthigh = rthigh
    else:
        logrtlow = np.log10(rtlow)
        logrthigh = np.log10(rthigh)
    dellow = 3.5
    delhigh = 5.0
                        
Mstar_min = Mstar - Mstar_err
Mstar_max = Mstar + Mstar_err

if (use_rclinear == 'yes'):
    outdir = outdir + 'rc_lin/'
if (use_rtlinear == 'yes'):
    outdir = outdir + 'rt_lin/'
if (codemode == 'run'):
    print 'Will write output to:', outdir

if (singlesplit == 'single'):
    ndim = n_betpars + n_nupars + n_mpars + 1
elif (singlesplit == 'split'):
    ndim = n_betpars*2 + n_nupars*2 + n_mpars + 1

if (propermotion == 'no'):
    if (virialshape == 'no'):
        sigp_fit = lambda r1, r2, nupars, Mpars, betpars, Mstar: \
            sigp(r1,r2,nu,Sigfunc,M,beta,betaf,nupars,Mpars,betpars,\
                 Mstar_rad,Mstar_prof,Mstar,Guse)
    else:
        sigp_fit_vs = lambda r1, r2, nupars, Mpars, betpars, Mstar: \
            sigp_vs(r1,r2,nu,Sigfunc,M,beta,betaf,nupars,Mpars,betpars,\
                    Mstar_rad,Mstar_prof,Mstar,Guse)
elif (propermotion == 'yes'):
    sigp_fit_prop = lambda r1, r2, nupars, Mpars, betpars, Mstar: \
        sigp_prop(r1,r2,nu,Sigfunc,M,beta,betaf,nupars,Mpars,betpars,\
                  Mstar_rad,Mstar_prof,Mstar,Guse)

#Calculate the tracer density profile:
if (singlesplit == 'single'):
    if (len(pfits) == 1):
        if (data_file_type == 'walker_style' or \
            data_file_type == 'walker_style_vs'):
            pfits = tracerfit(p0in,p0in_min,p0in_max,rbin_photfit,\
                              surfdenfit,surfdenerrfit)
        else:
            pfits = tracerfit(p0in,p0in_min,p0in_max,rbin_phot,\
                              surfden,surfdenerr)
    nupars_min = pfits*(1.0-tracertol)
    nupars_max = pfits*(1.0+tracertol)

    #Deal with negative parameters: 
    for i in range(len(nupars_min)):
        if (nupars_min[i] > nupars_max[i]):
            temp = nupars_max[i]
            nupars_max[i] = nupars_min[i]
            nupars_min[i] = temp
    
    if (testthreeplum == 'yes'):
        print 'Fit:', pfits
        test_surfden = threeplumsurf(ranal,pfits[0],pfits[1],pfits[2],\
                                     pfits[3],pfits[4],pfits[5])
        plt.figure()
        plt.loglog()
        plt.errorbar(rbin_phot,surfden,surfdenerr)
        plt.plot(ranal,test_surfden)
        plt.show()
        sys.exit()
elif (singlesplit == 'split'):
    pfits1 = tracerfit(p0in,p0in_min,p0in_max,\
                           rbin1_phot,surfden1,surfdenerr1)
    pfits2 = tracerfit(p0in,p0in_min,p0in_max,\
                           rbin2_phot,surfden2,surfdenerr2)
    nupars1_min = pfits1*(1.0-tracertol)
    nupars1_max = pfits1*(1.0+tracertol)
    nupars2_min = pfits2*(1.0-tracertol)
    nupars2_max = pfits2*(1.0+tracertol)

    #Deal with negative parameters:
    for i in range(len(nupars1_min)):
        if (nupars1_min[i] > nupars1_max[i]):
            temp = nupars1_max[i]
            nupars1_max[i] = nupars1_min[i]
            nupars1_min[i] = temp
    for i in range(len(nupars2_min)):
        if (nupars2_min[i] > nupars2_max[i]):
            temp = nupars2_max[i]
            nupars2_max[i] = nupars2_min[i]
            nupars2_min[i] = temp    

    if (testthreeplum == 'yes'):
        test_surfden1 = threeplumsurf(ranal,pfits1[0],pfits1[1],pfits1[2],\
                                     pfits1[3],pfits1[4],pfits1[5])
        test_surfden2 = threeplumsurf(ranal,pfits2[0],pfits2[1],pfits2[2],\
                                     pfits2[3],pfits2[4],pfits2[5])
        plt.figure()
        plt.loglog()
        plt.errorbar(rbin1_phot,surfden1,surfdenerr1)
        plt.plot(ranal,test_surfden1)
        plt.errorbar(rbin2_phot,surfden2,surfdenerr2)
        plt.plot(ranal,test_surfden2)
        plt.show()
        sys.exit()

if (virialshape == 'yes' and singlesplit == 'single'):
    nmonte = 1000
    rbin_kin_tmp, sigpmean_tmp, sigperr_tmp, \
        vs1bin, vs2bin, vs1err, vs2err = \
        calc_virial_moments(Rhalf,nmonte,R,vz,vzerr,ms,\
                            pfits)

    #If the second moment not loaded, set it:
    if (data_file_kin == ''):
        rbin_kin = rbin_kin_tmp
        sigpmean = sigpmean_tmp
        sigperr = sigperr_tmp
        
        if (include_jardel == 'yes'):
            #Add Virus-P data:
            print 'Adding in Virus-P data point:'
            rbin_kin = np.concatenate((np.array([rjar]),rbin_kin))
            sigpmean = np.concatenate((np.array([sigpjar]),sigpmean))
            sigperr = np.concatenate((np.array([sigperrjar]),sigperr))

### EMCEE FITTING CODE ### 

#Set up walkers:
if (codemode == 'run'):
    print 'Running in fitting mode ... '
    pos = np.zeros((nwalkers, ndim), dtype='float')

    if (blobstart == 'no'):
        if (singlesplit == 'single'):
            pos[:,0] = np.random.uniform(bet0min,bet0max,nwalkers)
            pos[:,1] = np.random.uniform(betinfmin,betinfmax,nwalkers)
            pos[:,2] = np.random.uniform(betr0min,betr0max,nwalkers)
            pos[:,3] = np.random.uniform(betnmin,betnmax,nwalkers)
            if (rotation == 'yes'):
                pos[:,4] = np.random.uniform(Arotmin,Arotmax,nwalkers)
            for i in range(len(nupars_min)):
                pos[:,n_betpars+i] = \
                    np.random.uniform(nupars_min[i],nupars_max[i],nwalkers)
            pos[:,n_betpars+nu_components*2] = \
                np.random.uniform(logM200low,logM200high,nwalkers)
            pos[:,n_betpars+nu_components*2+1] = \
                np.random.uniform(clow,chigh,nwalkers)
            pos[:,n_betpars+nu_components*2+2] = \
                np.random.uniform(logrclow,logrchigh,nwalkers)
            pos[:,n_betpars+nu_components*2+3] = \
                np.random.uniform(nlow,nhigh,nwalkers)
            pos[:,n_betpars+nu_components*2+4] = \
                np.random.uniform(logrtlow,logrthigh,nwalkers)
            pos[:,n_betpars+nu_components*2+5] = \
                np.random.uniform(dellow,delhigh,nwalkers)
            pos[:,ndim-1] = \
                np.random.uniform(Mstar_min,Mstar_max,nwalkers)
        elif (singlesplit == 'split'):
            print 'Not yet implmeneted. Bye.'
            sys.exit(0)
    else:
        #Start the blob in a reasonable part of parameter space:
        bet0min_start = 0.0
        bet0max_start = 0.2
        betr0min_start = np.log10(Rhalf)
        betr0max_start = np.log10(Rhalf*1.2)
        betnmin_start = 2.0
        betnmax_start = 2.2
        betinfmin_start = 0.0
        betinfmax_start = 0.2
        logM200low_start = np.log10(1e9)
        logM200high_start = np.log10(1e9*1.1)
        clow_start = 10.0
        chigh_start = 12.0
        nlow_start = nhigh-0.01
        nhigh_start = nhigh
        logrclow_start = logrclow
        logrchigh_start = logrchigh
        logrtlow_start = logrtlow
        logrthigh_start = logrthigh
        dellow_start = 3.5
        delhigh_start = 3.6

        if (singlesplit == 'single'):
            pos[:,0] = np.random.uniform(bet0min_start,bet0max_start,nwalkers)
            pos[:,1] = np.random.uniform(betinfmin_start,betinfmax_start,nwalkers)
            pos[:,2] = np.random.uniform(betr0min_start,betr0max_start,nwalkers)
            pos[:,3] = np.random.uniform(betnmin_start,betnmax_start,nwalkers)
            if (rotation == 'yes'):
                pos[:,4] = np.random.uniform(Arotmin,Arotmax,nwalkers)
            for i in range(len(nupars_min)):
                pos[:,n_betpars+i] = \
                    np.random.uniform(nupars_min[i],nupars_max[i],nwalkers)
            pos[:,n_betpars+nu_components*2] = \
                np.random.uniform(logM200low_start,logM200high_start,nwalkers)
            pos[:,n_betpars+nu_components*2+1] = \
                np.random.uniform(clow_start,chigh_start,nwalkers)
            pos[:,n_betpars+nu_components*2+2] = \
                np.random.uniform(logrclow_start,logrchigh_start,nwalkers)
            pos[:,n_betpars+nu_components*2+3] = \
                np.random.uniform(nlow_start,nhigh_start,nwalkers)
            pos[:,n_betpars+nu_components*2+4] = \
                np.random.uniform(logrtlow_start,logrthigh_start,nwalkers)
            pos[:,n_betpars+nu_components*2+5] = \
                np.random.uniform(dellow_start,delhigh_start,nwalkers)
            pos[:,ndim-1] = \
                np.random.uniform(Mstar_min,Mstar_max,nwalkers) 
        elif (singlesplit == 'split'):
            print 'Not yet implemented yet. Sorry. Bye.'
            sys.exit(0)

    #Set up fitting function and priors: 
    if (singlesplit == 'single'):
        if (propermotion == 'no'):
            if (virialshape == 'no'):
                x1 = rbin_phot
                x2 = rbin_kin
                y1 = surfden
                y1err = surfdenerr
                y2 = sigpmean
                y2err = sigperr
            else:
                x1 = rbin_phot
                x2 = rbin_kin
                y1 = surfden
                y1err = surfdenerr
                y2 = sigpmean
                y2err = sigperr
                y3 = vs1bin
                y3err = vs1err
                y4 = vs2bin
                y4err = vs2err
        elif (propermotion == 'yes'):
            x1 = rbin_phot
            x2 = rbin_kin
            y1 = surfden
            y1err = surfdenerr
            y2 = sigpmean
            y2err = sigperr
            y3 = sigpmr
            y3err = sigpmrerr
            y4 = sigpmt
            y4err = sigpmterr
    elif (singlesplit == 'split'):
        if (virialshape == 'no'):
            x1 = rbin1_phot
            x2 = rbin1_kin
            y1 = surfden1
            y1err = surfdenerr1
            y2 = sigpmean1
            y2err = sigperr1
            x3 = rbin2_phot
            x4 = rbin2_kin
            y3 = surfden2
            y3err = surfdenerr2
            y4 = sigpmean2
            y4err = sigperr2
        else:
            x1 = rbin1_phot
            x2 = rbin1_kin
            y1 = surfden1
            y1err = surfdenerr1
            y2 = sigpmean1
            y2err = sigperr1
            x3 = rbin2_phot
            x4 = rbin2_kin
            y3 = surfden2
            y3err = surfdenerr2
            y4 = sigpmean2
            y4err = sigperr2
            y5 = vs1bin1
            y5err = vs1err1
            y6 = vs2bin1
            y6err = vs2err1
            y7 = vs1bin2
            y7err = vs1err2
            y8 = vs2bin2
            y8err = vs2err2
    if (singlesplit == 'single'):
        if (rotation == 'no'):
            lnprior = lambda theta: \
                lnprior_set(theta,n_betpars,bet0min,bet0max,\
                                betinfmin,betinfmax,\
                                betr0min,betr0max,betnmin,betnmax,\
                                nu_components,nupars_min,nupars_max,\
                                n_mpars,logM200low,logM200high,\
                                clow,chigh,logrclow,logrchigh,\
                                nlow,nhigh,logrtlow,logrthigh,\
                                dellow,delhigh,\
                                Mstar_min,Mstar_max)
        else:
            lnprior = lambda theta: \
                lnprior_set(theta,n_betpars,bet0min,bet0max,\
                        betinfmin,betinfmax,\
                        betr0min,betr0max,betnmin,betnmax,\
                        Arotmin,Arotmax,\
                        nu_components,nupars_min,nupars_max,\
                        n_mpars,logM200low,logM200high,\
                        clow,chigh,logrclow,logrchigh,\
                        nlow,nhigh,logrtlow,logrthigh,\
                        dellow,delhigh,\
                        Mstar_min,Mstar_max)
            
    elif (singlesplit == 'split'):
        lnprior = lambda theta: \
            lnprior_set(theta,n_betpars,bet0min,bet0max,betinfmin,betinfmax,\
                            betr0min,betr0max,betnmin,betnmax,\
                            nu_components,nupars1_min,nupars1_max,\
                            nupars2_min,nupars2_max,\
                            n_mpars,log_rho0_min,log_rho0_max,\
                            gam_min,gam_max,Mstar_min,Mstar_max)

    print 'Running chains ... '
    if (singlesplit == 'single'):
        if (propermotion == 'no'):
            if (virialshape == 'no'):
                sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, \
                                     args=(x1, x2, y1, y1err, y2, y2err))
            else:
                sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, \
                                     args=(x1, x2, y1, y1err, y2, y2err, \
                                           y3, y3err, y4, y4err))

        elif (propermotion == 'yes'):
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, \
                                       args=(x1, x2, y1, y1err, y2, y2err, \
                                             y3, y3err, y4, y4err))
    elif (singlesplit == 'split'):
        if (virialshape == 'no'):
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, \
                                args=(x1, x2, x3, x4, y1, y1err, y2, y2err, \
                                      y3, y3err, y4, y4err))
        else:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, \
                                args=(x1, x2, x3, x4, y1, y1err, y2, y2err, \
                                      y3, y3err, y4, y4err, y5, y5err, \
                                      y6, y6err, y7, y7err, y8, y8err))

    sampler.run_mcmc(pos, nmodels)

    #Store the output (including the data):
    print 'Writing data to file ... '
    if (singlesplit == 'single'):
        if (propermotion == 'no'):
            f = open(outdir+'output_sigp.txt','w')
            for i in range(len(rbin_kin)):
                f.write('%f %f %f\n' % \
                        (rbin_kin[i], sigpmean[i], sigperr[i]))
            f.close()
            f = open(outdir+'output_surfden.txt','w')
            for i in range(len(rbin_phot)):
                f.write('%f %f %f\n' % \
                        (rbin_phot[i], surfden[i], surfdenerr[i]))
            f.close()
        elif (propermotion == 'yes'):
            f = open(outdir+'output_sigs.txt','w')
            for i in range(len(rbin_kin)):
                f.write('%f %f %f %f %f %f %f\n' % \
                        (rbin_kin[i],sigpmean[i],sigperr[i],\
                         sigpmr[i],sigpmrerr[i],sigpmt[i],sigpmterr[i]))
            f.close()
            f = open(outdir+'output_surfden.txt','w')
            for i in range(len(rbin_phot)):
                f.write('%f %f %f\n' % \
                        (rbin_phot[i], surfden[i], surfdenerr[i]))
            f.close()
    elif (singlesplit == 'split'):
        f = open(outdir+'output_sigp1.txt','w')
        for i in range(len(rbin1_kin)):
            f.write('%f %f %f\n' % \
                    (rbin1_kin[i], sigpmean1[i], sigperr1[i]))
        f.close()
        f = open(outdir+'output_surfden1.txt','w')
        for i in range(len(rbin1_phot)):
            f.write('%f %f %f\n' % \
                    (rbin1_phot[i], surfden1[i], surfdenerr1[i]))
        f.close()
        f = open(outdir+'output_sigp2.txt','w')
        for i in range(len(rbin2_kin)):
            f.write('%f %f %f\n' % \
                    (rbin2_kin[i], sigpmean2[i], sigperr2[i]))
        f.close()
        f = open(outdir+'output_surfden2.txt','w')
        for i in range(len(rbin2_phot)):
            f.write('%f %f %f\n' % \
                    (rbin2_phot[i], surfden2[i], surfdenerr2[i]))
        f.close()

    burn = np.int(0.5*nmodels)
    chisq = -2.0 * sampler.lnprobability[:, burn:].reshape(-1)
    par_test = np.zeros((len(chisq),ndim), dtype='float')
    for i in range(ndim):
        par_test[:,i] = sampler.chain[:, burn:, i].reshape(-1)

    f = open(outdir+'Boutput_chain.txt','w')
    for i in range(len(chisq)):
        outstr = str(chisq[i]) + ' '
        for j in range(ndim):
            outstr = outstr + str(par_test[i,j]) + ' '
        outstr = outstr + '\n'
        f.write(outstr)
    f.close()


### PLOTTING CODE ###
elif (codemode == 'plot'):
    print 'Running in plotting mode ... '
    print 'Loading data from:', outdir

    #Read in the data:
    if (singlesplit == 'single'):
        if (propermotion == 'no'):
            f = open(outdir+'output_sigp.txt','r')
            data_in = np.genfromtxt(f)
            f.close()
            rbin_kin = data_in[:,0]
            sigpmean = data_in[:,1]
            sigperr = data_in[:,2]
            f = open(outdir+'output_surfden.txt','r')
            data_in = np.genfromtxt(f)
            f.close()
            rbin_phot = data_in[:,0]
            surfden = data_in[:,1]
            surfdenerr = data_in[:,2]
        elif (propermotion == 'yes'):
            f = open(outdir+'output_sigs.txt','r')
            data_in = np.genfromtxt(f)
            f.close()
            rbin_kin = data_in[:,0]
            sigpmean = data_in[:,1]
            sigperr = data_in[:,2]
            sigpmr = data_in[:,3]
            sigpmrerr = data_in[:,4]
            sigpmt = data_in[:,5]
            sigpmterr = data_in[:,6]

            f = open(outdir+'output_surfden.txt','r')
            data_in = np.genfromtxt(f)
            f.close()
            rbin_phot = data_in[:,0]
            surfden = data_in[:,1]
            surfdenerr = data_in[:,2]

        #Define binning to use for the mass distribution
        #plots:
        rbin = rbin_phot
    elif (singlesplit == 'split'):
        f = open(outdir+'output_sigp1.txt','r')
        data_in = np.genfromtxt(f)
        f.close()
        rbin1_kin = data_in[:,0]
        sigpmean1 = data_in[:,1]
        sigperr1 = data_in[:,2]
        f = open(outdir+'output_surfden1.txt','r')
        data_in = np.genfromtxt(f)
        f.close()
        rbin1_phot = data_in[:,0]
        surfden1 = data_in[:,1]
        surfdenerr1 = data_in[:,2]
        f = open(outdir+'output_sigp2.txt','r')
        data_in = np.genfromtxt(f)
        f.close()
        rbin2_kin = data_in[:,0]
        sigpmean2 = data_in[:,1]
        sigperr2 = data_in[:,2]
        f = open(outdir+'output_surfden2.txt','r')
        data_in = np.genfromtxt(f)
        f.close()
        rbin2_phot = data_in[:,0]
        surfden2 = data_in[:,1]
        surfdenerr2 = data_in[:,2]

        #Define binning to use for the plots:
        rbin = rbin1_phot
        rbin1 = rbin1_phot
        rbin2 = rbin2_phot
        
    #Fix up if the plot range is wider than rbin!
    if (rbin_inner > 0):
        if (rbin_inner < np.min(rbin)):
            rleft = rbin_inner
        else:
            rleft = np.min(rbin)
        if (rbin_outer > np.max(rbin)):
            rright = rbin_outer
        else:
            rright = np.max(rbin)
        rbin = np.logspace(np.log10(rleft),np.log10(rright),\
                               len(rbin))

    f = open(outdir+'Boutput_chain.txt','r')
    data_in = np.genfromtxt(f)
    chisq = data_in[:,0]
    par_test = np.zeros((len(chisq),ndim), dtype='float')
    for i in range(1,ndim+1):
        par_test[:,i-1] = data_in[:,i]
    f.close()

    #Make sure no *really* bad models remain in the chains. In practice,
    #this cut makes no difference to the end result.
    min_chisq = np.min(chisq)
    index = np.where(chisq < min_chisq*500.0)[0]
    print 'Min/max chisq:', np.min(chisq), np.max(chisq)

    #Select models for the min-max lines. These show the range of 
    #"Good" models that lie at extremum beta*:
    chisqind = chisq[index]
    index_chisq = np.argsort(chisqind,axis=0)
    median_chisq = chisqind[index_chisq[np.int(len(chisqind)/2.)]]
    chisq_max = median_chisq*1.0
    index_good = np.where(chisq < chisq_max)[0]
    minmaxtype = 'allbet'

    #Cut the confidence intervals from the chains:
    nsamples = 5000
    sample_choose = index[np.random.randint(len(index), \
                                            size=nsamples)]

    #Make corner plot:
    if (makecorner == 'yes'):
        fig = corner.corner(par_test[sample_choose,:])
        plt.savefig(outdir+'triangle.png')                    

        #Plot the chisq distribution too:
        fig = plt.figure(figsize=(figx*5,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)

        plt.semilogy()

        plt.xlabel(r'$\#{\rm steps}$',\
                   fontsize=myfontsize)
        plt.ylabel(r'$\chi^2$',\
                   fontsize=myfontsize)
        
        plt.plot(chisq,color='black')
        plt.plot(chisqind,color='blue')
        plt.savefig(outdir+'output_chisq.png')
        
    if (singlesplit == 'single'):
        if (minmaxtype == 'allbet'):
            minmax_choose = np.array([\
                    index_good[\
                        np.where(\
                            par_test[index_good,0]+\
                                par_test[index_good,1] == \
                                np.max(par_test[index_good,0]+\
                                           par_test[index_good,1]))[0][0]],\
                        index_good[\
                        np.where(\
                            par_test[index_good,0]+\
                                par_test[index_good,1] == \
                                np.min(par_test[index_good,0]+\
                                           par_test[index_good,1]))[0][0]]])
        elif (minmaxtype == 'betinf'):
            minmax_choose = np.array([\
                    index_good[\
                        np.where(\
                            par_test[index_good,1] == \
                                np.max(par_test[index_good,1]))[0][0]],\
                        index_good[\
                        np.where(\
                            par_test[index_good,1] == \
                                np.min(par_test[index_good,1]))[0][0]]])
        elif (minmaxtype == 'gamcen'):
            minmax_choose = np.array([\
                    index_good[\
                        np.where(\
                            par_test[index_good,n_betpars+nu_components*2+1] == \
                                np.max(par_test[index_good,n_betpars+nu_components*2+1]))[0][0]],\
                        index_good[\
                        np.where(\
                            par_test[index_good,n_betpars+nu_components*2+1] == \
                                np.min(par_test[index_good,n_betpars+nu_components*2+1]))[0][0]]])
        elif (minmaxtype == 'allgam'):
            minmax_choose = np.array([\
                    index_good[\
                        np.where(\
                            par_test[index_good,n_betpars+nu_components*2+1]+\
                            par_test[index_good,n_betpars+nu_components*2+2]+\
                            par_test[index_good,n_betpars+nu_components*2+3]+\
                            par_test[index_good,n_betpars+nu_components*2+4]+\
                            par_test[index_good,n_betpars+nu_components*2+5] == \
                            np.max(\
                                par_test[index_good,n_betpars+nu_components*2+1]+\
                                par_test[index_good,n_betpars+nu_components*2+2]+\
                                par_test[index_good,n_betpars+nu_components*2+3]+\
                                par_test[index_good,n_betpars+nu_components*2+4]+\
                                par_test[index_good,n_betpars+nu_components*2+5]))[0][0]],\
                     index_good[\
                         np.where(\
                            par_test[index_good,n_betpars+nu_components*2+1]+\
                            par_test[index_good,n_betpars+nu_components*2+2]+\
                            par_test[index_good,n_betpars+nu_components*2+3]+\
                            par_test[index_good,n_betpars+nu_components*2+4]+\
                            par_test[index_good,n_betpars+nu_components*2+5] == \
                            np.min(\
                                par_test[index_good,n_betpars+nu_components*2+1]+\
                                par_test[index_good,n_betpars+nu_components*2+2]+\
                                par_test[index_good,n_betpars+nu_components*2+3]+\
                                par_test[index_good,n_betpars+nu_components*2+4]+\
                                par_test[index_good,n_betpars+nu_components*2+5]))[0][0]]])
    elif (singlesplit == 'split'):
        if (minmaxtype == 'allbet'):
            minmax_choose = np.array([\
                    index_good[\
                        np.where(\
                            par_test[index_good,0]+\
                                par_test[index_good,1]+\
                                par_test[index_good,4]+\
                                par_test[index_good,5] == \
                                np.max(par_test[index_good,0]+\
                                           par_test[index_good,1]+\
                                           par_test[index_good,4]+
                                       par_test[index_good,5]))[0][0]],\
                        index_good[\
                            np.where(\
                                 par_test[index_good,0]+\
                                     par_test[index_good,1]+\
                                     par_test[index_good,4]+\
                                     par_test[index_good,5] == \
                                     np.min(par_test[index_good,0]+\
                                                par_test[index_good,1]+\
                                                par_test[index_good,4]+\
                                                par_test[index_good,5]))[0][0]]])
        elif (minmaxtype == 'betinf'):
            minmax_choose = np.array([\
                    index_good[\
                        np.where(\
                                par_test[index_good,1]+\
                                par_test[index_good,5] == \
                                np.max(\
                                par_test[index_good,1]+\
                                    par_test[index_good,5]))[0][0]],\
                        index_good[\
                            np.where(\
                                par_test[index_good,1]+\
                                    par_test[index_good,5] == \
                                    np.min(\
                                par_test[index_good,1]+\
                                    par_test[index_good,5]))[0][0]]])            
    M_int = np.zeros((7,len(rbin)))
    rho_int = np.zeros((7,len(rbin)))
    dlnrhodlnr_int = np.zeros((7,len(rbin)))
    Mstar_int = np.zeros((7,len(Mstar_rad)))
    Mdynrat_int = np.zeros((7,len(Mstar_rad)))
    nu_int = np.zeros((7,len(Mstar_rad)))
    if (calc_Jfac == 'yes'):
        J_int = np.zeros(7)
    
    Mminmax = np.zeros((len(rbin),2))
    rhominmax = np.zeros((len(rbin),2))
    dlnrhodlnrminmax = np.zeros((len(rbin),2))

    Mstore = np.zeros((len(rbin),nsamples))
    rhostore = np.zeros((len(rbin),nsamples))
    dlnrhodlnrstore = np.zeros((len(rbin),nsamples))
    Mstarstore = np.zeros((len(Mstar_rad),nsamples))
    Mdynratstore = np.zeros((len(Mstar_rad),nsamples))
    nustore = np.zeros((len(Mstar_rad),nsamples))
    M200store = np.zeros(nsamples)
    vmaxstore = np.zeros(nsamples)
    cstore = np.zeros(nsamples)
    rcstore = np.zeros(nsamples)
    nstore = np.zeros(nsamples)
    rtstore = np.zeros(nsamples)
    delstore = np.zeros(nsamples)
    if (calc_Jfac == 'yes'):
        Jstore = np.zeros(nsamples)       
    
    if (singlesplit == 'single'):
        bet_int = np.zeros((7,len(rbin)))
        betstar_int = np.zeros((7,len(rbin)))
        Sig_int = np.zeros((7,len(rbin)))
        sigp_int = np.zeros((7,len(rbin)))
        if (virialshape == 'yes'):
            vs1_int = np.zeros((7,1))
            vs2_int = np.zeros((7,1))
        if (propermotion == 'yes'):
            sigpmr_int = np.zeros((7,len(rbin)))
            sigpmt_int = np.zeros((7,len(rbin)))
        betstore = np.zeros((len(rbin),nsamples))
        betstarstore = np.zeros((len(rbin),nsamples))
        Sigstore = np.zeros((len(rbin),nsamples))
        sigpstore = np.zeros((len(rbin),nsamples))
        betminmax = np.zeros((len(rbin),2))
        betstarminmax = np.zeros((len(rbin),2))
        Sigminmax = np.zeros((len(rbin),2))
        sigpminmax = np.zeros((len(rbin),2))
        if (virialshape == 'yes'):
            vs1store = np.zeros(nsamples)
            vs2store = np.zeros(nsamples)
        if (propermotion == 'yes'):
            sigpmrstore = np.zeros((len(rbin),nsamples))
            sigpmtstore = np.zeros((len(rbin),nsamples))
            sigpmrminmax = np.zeros((len(rbin),2))
            sigpmtminmax = np.zeros((len(rbin),2))
        if (rotation == 'yes'):
            Arotstore = np.zeros(nsamples)
    elif (singlesplit == 'split'):
        bet1_int = np.zeros((7,len(rbin1)))
        betstar1_int = np.zeros((7,len(rbin1)))
        Sig1_int = np.zeros((7,len(rbin1)))
        sigp1_int = np.zeros((7,len(rbin1)))
        bet2_int = np.zeros((7,len(rbin2)))
        betstar2_int = np.zeros((7,len(rbin2)))
        Sig2_int = np.zeros((7,len(rbin2)))
        sigp2_int = np.zeros((7,len(rbin2)))
        if (virialshape == 'yes'):
            vs1_1int = np.zeros((7,1))
            vs2_1int = np.zeros((7,1))
            vs1_2int = np.zeros((7,1))
            vs2_2int = np.zeros((7,1))
        bet1store = np.zeros((len(rbin1),nsamples))
        betstar1store = np.zeros((len(rbin1),nsamples))
        Sig1store = np.zeros((len(rbin1),nsamples))
        sigp1store = np.zeros((len(rbin1),nsamples))
        bet2store = np.zeros((len(rbin2),nsamples))
        betstar2store = np.zeros((len(rbin2),nsamples))
        Sig2store = np.zeros((len(rbin2),nsamples))
        sigp2store = np.zeros((len(rbin2),nsamples))
        if (virialshape == 'yes'):
            vs1_1store = np.zeros(nsamples)
            vs2_1store = np.zeros(nsamples)
            vs1_2store = np.zeros(nsamples)
            vs2_2store = np.zeros(nsamples)
        bet1minmax = np.zeros((len(rbin1),2))
        betstar1minmax = np.zeros((len(rbin1),2))
        Sig1minmax = np.zeros((len(rbin1),2))
        sigp1minmax = np.zeros((len(rbin1),2))
        bet2minmax = np.zeros((len(rbin2),2))
        betstar2minmax = np.zeros((len(rbin2),2))
        Sig2minmax = np.zeros((len(rbin2),2))
        sigp2minmax = np.zeros((len(rbin2),2))
    for i in range(nsamples):
        theta = par_test[sample_choose[i],:]
        if (singlesplit == 'single'):
            betpars = theta[0:n_betpars]
            nupars = theta[n_betpars:n_betpars+nu_components*2]
            Mpars = theta[n_betpars+nu_components*2:\
                          n_betpars+nu_components*2+n_mpars]
            Mstar = theta[ndim-1]
            nuparsu = np.array(nupars)
            Mparsu = np.array(Mpars)
            Mparsu[0] = 10.**Mpars[0]
            if (use_rclinear == 'no'):
                Mparsu[2] = 10.**Mpars[2]
            if (use_rtlinear == 'no'):
                Mparsu[4] = 10.**Mpars[4]
            if (rotation == 'yes'):
                Arotstore[i] = betpars[4]
        elif (singlesplit == 'split'):
            betpars1 = theta[0:n_betpars]
            betpars2 = theta[n_betpars:2*n_betpars]
            nupars1 = theta[n_betpars*2:n_betpars*2+nu_components*2]
            nupars2 = theta[n_betpars*2+nu_components*2:\
                            n_betpars*2+nu_components*2*2]
            Mpars = theta[n_betpars*2+nu_components*2*2:\
                          n_betpars*2+nu_components*2*2+n_mpars]
            Mstar = theta[ndim-1]
            nupars1u = np.array(nupars1)
            nupars2u = np.array(nupars2)
            Mparsu = np.array(Mpars)
            Mparsu[0] = 10.**Mpars[0]
            if (use_rclinear == 'no'):
                Mparsu[2] = 10.**Mpars[2]
            if (use_rtlinear == 'no'):
                Mparsu[4] = 10.**Mpars[4]

        #Calculate all profiles we want to plot:
        if (singlesplit == 'single'):
            if (propermotion == 'no'):
                if (virialshape == 'no'):
                    sigr2, Sig, sigLOS2 = sigp_fit(rbin,rbin,nuparsu,\
                                                   Mparsu,betpars,Mstar)
                else:
                    sigr2, Sig, sigLOS2, vs1, vs2 = sigp_fit_vs(rbin,rbin,nuparsu,\
                                                      Mparsu,betpars,Mstar)
            elif (propermotion == 'yes'):
                sigr2, Sig, sigLOS2, sigpmr2, sigpmt2 = \
                    sigp_fit_prop(rbin,rbin,nuparsu,Mparsu,betpars,Mstar)
            Mr = M(rbin,Mparsu)
            betar = beta(rbin,betpars)
            rhor = rho(rbin,Mparsu)
            dlnrhodlnrr = dlnrhodlnr(rbin,Mparsu)
            Mstarr = Mstar_prof*Mstar
            nu_mass_r = multiplummass(Mstar_rad,nuparsu)

            Mstore[:,i] = Mr
            betstore[:,i] = betar
            betstarstore[:,i] = betar/(2.0-betar)
            sigpstore[:,i] = np.sqrt(sigLOS2)/1000. 
            Sigstore[:,i] = Sig
            rhostore[:,i] = rhor
            dlnrhodlnrstore[:,i] = dlnrhodlnrr
            Mstarstore[:,i] = Mstarr
            Mdynratstore[:,i] = M(Mstar_rad,Mparsu)//Mstarr
            nustore[:,i] = nu_mass_r

            vmaxstore[i] = vmax_func(Mparsu[0],Mparsu[1],h)
            M200store[i] = Mparsu[0]
            cstore[i] = Mparsu[1]
            rcstore[i] = Mparsu[2]
            nstore[i] = Mparsu[3]
            rtstore[i] = Mparsu[4]
            delstore[i] = Mparsu[5]

            if (calc_Jfac == 'yes'):
                alpha_rmax = dgal*alpha_Jfac_deg/deg
                Jstore[i] = get_Juse(Mparsu,dgal,alpha_rmax)
            if (virialshape == 'yes'):
                vs1store[i] = vs1/1.0e12
                vs2store[i] = vs2/1.0e12
            if (propermotion == 'yes'):
                sigpmrstore[:,i] = np.sqrt(sigpmr2)/1000.
                sigpmtstore[:,i] = np.sqrt(sigpmt2)/1000.
        elif (singlesplit == 'split'):
            if (virialshape == 'no'):
                sigr2_1, Sig_1, sigLOS2_1 = \
                    sigp_fit(rbin1,rbin1,nupars1u,Mparsu,betpars1,Mstar)
                sigr2_2, Sig_2, sigLOS2_2 = \
                    sigp_fit(rbin2,rbin2,nupars2u,Mparsu,betpars2,Mstar)
            else:
                sigr2_1, Sig_1, sigLOS2_1, vs1_1, vs2_1 = \
                    sigp_fit_vs(rbin1,rbin1,nupars1u,Mparsu,betpars1,Mstar)
                sigr2_2, Sig_2, sigLOS2_2, vs1_2, vs2_2 = \
                    sigp_fit_vs(rbin2,rbin2,nupars2u,Mparsu,betpars2,Mstar)
            Mr = M(rbin,Mparsu)
            betar1 = beta(rbin1,betpars1)
            betar2 = beta(rbin2,betpars2)
            rhor = rho(rbin,Mparsu)
            dlnrhodlnrr = dlnrhodlnr(rbin,Mparsu)
            Mstarr = Mstar_prof*Mstar

            Mstore[:,i] = Mr
            bet1store[:,i] = betar1
            bet2store[:,i] = betar2
            betstar1store[:,i] = betar1/(2.0-betar1)
            betstar2store[:,i] = betar2/(2.0-betar2)
            sigp1store[:,i] = np.sqrt(sigLOS2_1)/1000.
            sigp2store[:,i] = np.sqrt(sigLOS2_2)/1000.
            Sig1store[:,i] = Sig_1
            Sig2store[:,i] = Sig_2
            rhostore[:,i] = rhor
            dlnrhodlnrstore[:,i] = dlnrhodlnrr
            Mstarstore[:,i] = Mstarr
            Mdynratstore[:,i] = Mr/Mstarr
            nustore[:,i] = nu_mass_r
            if (calc_Jfac == 'yes'):
                alpha_rmax = dgal*alpha_Jfac_deg/deg
                Jstore[i] = get_J(Mparsu,dgal,alpha_rmax)

            if (virialshape == 'yes'):
                vs1_1store[i] = vs1_1/1.0e12
                vs2_1store[i] = vs2_1/1.0e12
                vs1_2store[i] = vs1_2/1.0e12
                vs2_2store[i] = vs2_2/1.0e12

        #Solve for confidence intervals for each of these:
        for j in range(len(rbin)):
            M_int[0,j], M_int[1,j], M_int[2,j], M_int[3,j], \
                M_int[4,j], M_int[5,j], M_int[6,j] = \
                calcmedquartnine(Mstore[j,:])
            rho_int[0,j], rho_int[1,j], rho_int[2,j], rho_int[3,j], \
                rho_int[4,j], rho_int[5,j], rho_int[6,j] = \
                calcmedquartnine(rhostore[j,:])
            dlnrhodlnr_int[0,j], dlnrhodlnr_int[1,j], dlnrhodlnr_int[2,j], \
                dlnrhodlnr_int[3,j], \
                dlnrhodlnr_int[4,j], \
                dlnrhodlnr_int[5,j], \
                dlnrhodlnr_int[6,j] = \
                calcmedquartnine(dlnrhodlnrstore[j,:])
        for j in range(len(Mstar_rad)):
            Mstar_int[0,j], Mstar_int[1,j], Mstar_int[2,j], \
                Mstar_int[3,j], \
                Mstar_int[4,j], \
                Mstar_int[5,j], \
                Mstar_int[6,j] = \
                calcmedquartnine(Mstarstore[j,:])
        for j in range(len(Mstar_rad)):
            Mdynrat_int[0,j], Mdynrat_int[1,j], Mdynrat_int[2,j], \
                Mdynrat_int[3,j], \
                Mdynrat_int[4,j], \
                Mdynrat_int[5,j], \
                Mdynrat_int[6,j] = \
                calcmedquartnine(Mdynratstore[j,:])
        for j in range(len(Mstar_rad)):
            nu_int[0,j], nu_int[1,j], nu_int[2,j], \
                nu_int[3,j], \
                nu_int[4,j], \
                nu_int[5,j], \
                nu_int[6,j] = \
                calcmedquartnine(nustore[j,:])
        if (calc_Jfac == 'yes'):
            J_int = \
                    calcmedquartnine(Jstore[:])
            
        if (singlesplit == 'single'):
            for j in range(len(rbin)):
                bet_int[0,j], bet_int[1,j], bet_int[2,j], bet_int[3,j], \
                    bet_int[4,j], \
                    bet_int[5,j], \
                    bet_int[6,j] = \
                    calcmedquartnine(betstore[j,:])
                betstar_int[0,j], betstar_int[1,j], \
                    betstar_int[2,j], betstar_int[3,j], \
                    betstar_int[4,j], \
                    betstar_int[5,j], \
                    betstar_int[6,j] = \
                    calcmedquartnine(betstarstore[j,:])
                sigp_int[0,j], sigp_int[1,j], sigp_int[2,j], sigp_int[3,j], \
                    sigp_int[4,j], \
                    sigp_int[5,j], \
                    sigp_int[6,j] = \
                    calcmedquartnine(sigpstore[j,:])
                Sig_int[0,j], Sig_int[1,j], Sig_int[2,j], Sig_int[3,j], \
                    Sig_int[4,j], \
                    Sig_int[5,j], \
                    Sig_int[6,j] = \
                    calcmedquartnine(Sigstore[j,:])
                if (propermotion == 'yes'):
                    sigpmr_int[0,j], sigpmr_int[1,j], \
                        sigpmr_int[2,j], sigpmr_int[3,j], \
                        sigpmr_int[4,j], \
                        sigpmr_int[5,j], \
                        sigpmr_int[6,j] = \
                        calcmedquartnine(sigpmrstore[j,:])
                    sigpmt_int[0,j], sigpmt_int[1,j], \
                        sigpmt_int[2,j], sigpmt_int[3,j], \
                        sigpmt_int[4,j], \
                        sigpmt_int[5,j], \
                        sigpmt_int[6,j] = \
                        calcmedquartnine(sigpmtstore[j,:])
            if (virialshape == 'yes'):
                vs1_int[0], vs1_int[1], \
                    vs1_int[2], vs1_int[3], \
                    vs1_int[4], \
                    vs1_int[5], \
                    vs1_int[6] = \
                    calcmedquartnine(vs1store[:])
                vs2_int[0], vs2_int[1], \
                    vs2_int[2], vs2_int[3], \
                    vs2_int[4], \
                    vs2_int[5], \
                    vs2_int[6] = \
                    calcmedquartnine(vs2store[:])
        elif (singlesplit == 'split'):
            for j in range(len(rbin1)):
                bet1_int[0,j], bet1_int[1,j], bet1_int[2,j], \
                    bet1_int[3,j], \
                    bet1_int[4,j], \
                    bet1_int[5,j], \
                    bet1_int[6,j] = \
                    calcmedquartnine(bet1store[j,:])
                betstar1_int[0,j], betstar1_int[1,j], betstar1_int[2,j], \
                    betstar1_int[3,j], \
                    betstar1_int[4,j], \
                    betstar1_int[5,j], \
                    betstar1_int[6,j] = \
                    calcmedquartnine(betstar1store[j,:])
                sigp1_int[0,j], sigp1_int[1,j], sigp1_int[2,j], \
                    sigp1_int[3,j], \
                    sigp1_int[4,j], \
                    sigp1_int[5,j], \
                    sigp1_int[6,j] = \
                    calcmedquartnine(sigp1store[j,:])
                Sig1_int[0,j], Sig1_int[1,j], Sig1_int[2,j], \
                    Sig1_int[3,j], \
                    Sig1_int[4,j], \
                    Sig1_int[5,j], \
                    Sig1_int[6,j] = \
                    calcmedquartnine(Sig1store[j,:])
            for j in range(len(rbin2)):
                bet2_int[0,j], bet2_int[1,j], bet2_int[2,j], \
                    bet2_int[3,j], \
                    bet2_int[4,j], \
                    bet2_int[5,j], \
                    bet2_int[6,j] = \
                    calcmedquartnine(bet2store[j,:])
                betstar2_int[0,j], betstar2_int[1,j], betstar2_int[2,j], \
                    betstar2_int[3,j], \
                    betstar2_int[4,j], \
                    betstar2_int[5,j], \
                    betstar2_int[6,j] = \
                    calcmedquartnine(betstar2store[j,:])
                sigp2_int[0,j], sigp2_int[1,j], sigp2_int[2,j], \
                    sigp2_int[3,j], \
                    sigp2_int[4,j], \
                    sigp2_int[5,j], \
                    sigp2_int[6,j] = \
                    calcmedquartnine(sigp2store[j,:])
                Sig2_int[0,j], Sig2_int[1,j], Sig2_int[2,j], \
                    Sig2_int[3,j], \
                    Sig2_int[4,j], \
                    Sig2_int[5,j], \
                    Sig2_int[6,j] = \
                    calcmedquartnine(Sig2store[j,:])
            if (virialshape == 'yes'):
                vs1_1int[0], vs1_1int[1], \
                    vs1_1int[2], vs1_1int[3], \
                    vs1_1int[4], \
                    vs1_1int[5], \
                    vs1_1int[6] = \
                    calcmedquartnine(vs1_1store[:])
                vs2_1int[0], vs2_1int[1], \
                    vs2_1int[2], vs2_1int[3], \
                    vs2_1int[4], \
                    vs2_1int[5], \
                    vs2_1int[6] = \
                    calcmedquartnine(vs2_1store[:])
                vs1_2int[0], vs1_2int[1], \
                    vs1_2int[2], vs1_2int[3], \
                    vs1_2int[4], \
                    vs1_2int[5], \
                    vs1_2int[6] = \
                    calcmedquartnine(vs1_2store[:])
                vs2_2int[0], vs2_2int[1], \
                    vs2_2int[2], vs2_2int[3], \
                    vs2_2int[4], \
                    vs2_2int[5], \
                    vs2_2int[6] = \
                    calcmedquartnine(vs2_2store[:])

    #Now minmax calc:
    for i in range(len(minmax_choose)):
        theta = par_test[minmax_choose[i],:]
        if (singlesplit == 'single'):
            betpars = theta[0:n_betpars]
            nupars = theta[n_betpars:n_betpars+nu_components*2]
            Mpars = theta[n_betpars+nu_components*2:\
                          n_betpars+nu_components*2+n_mpars]
            nuparsu = np.array(nupars)
            Mparsu = np.array(Mpars)
            Mparsu[0] = 10.**Mpars[0]
            if (use_rclinear == 'no'):
                Mparsu[2] = 10.**Mpars[2]
            if (use_rtlinear == 'no'):
                Mparsu[4] = 10.**Mpars[4]
        elif (singlesplit == 'split'):
            betpars1 = theta[0:n_betpars]
            betpars2 = theta[n_betpars:2*n_betpars]
            nupars1 = theta[n_betpars*2:n_betpars*2+nu_components*2]
            nupars2 = theta[n_betpars*2+nu_components*2:\
                            n_betpars*2+nu_components*2*2]
            Mpars = theta[n_betpars*2+nu_components*2*2:\
                          n_betpars*2+nu_components*2*2+n_mpars]
            Mstar = theta[ndim-1]
            nupars1u = np.array(nupars1)
            nupars2u = np.array(nupars2)
            Mparsu = np.array(Mpars)
            Mparsu[0] = 10.**Mpars[0]
            if (use_rclinear == 'no'):
                Mparsu[2] = 10.**Mpars[2]
            if (use_rtlinear == 'no'):
                Mparsu[4] = 10.**Mpars[4]

        if (singlesplit == 'single'):
            if (propermotion == 'no'):
                if (virialshape == 'no'):
                    sigr2, Sig, sigLOS2 = \
                        sigp_fit(rbin,rbin,nuparsu,Mparsu,betpars,Mstar)
                else:
                    sigr2, Sig, sigLOS2, vs1, vs2 = \
                        sigp_fit_vs(rbin,rbin,nuparsu,Mparsu,betpars,Mstar)
            elif (propermotion == 'yes'):
                sigr2, Sig, sigLOS2, sigpmr2, sigpmt2 = \
                    sigp_fit_prop(rbin,rbin,nuparsu,Mparsu,betpars,Mstar)
            Mr = M(rbin,Mparsu)
            betar = beta(rbin,betpars)
            rhor = rho(rbin,Mparsu)
            dlnrhodlnrr = dlnrhodlnr(rbin,Mparsu)

            Mminmax[:,i] = Mr
            betminmax[:,i] = betar
            betstarminmax[:,i] = betar/(2.0-betar)
            sigpminmax[:,i] = np.sqrt(sigLOS2)/1000.
            Sigminmax[:,i] = Sig
            rhominmax[:,i] = rhor
            dlnrhodlnrminmax[:,i] = dlnrhodlnrr
        elif (singlesplit == 'split'):
            if (virialshape == 'no'):
                sigr2_1, Sig_1, sigLOS2_1 = \
                    sigp_fit(rbin1,rbin1,nupars1u,Mparsu,betpars1,Mstar)
                sigr2_2, Sig_2, sigLOS2_2 = \
                    sigp_fit(rbin2,rbin2,nupars2u,Mparsu,betpars2,Mstar)
            else:
                sigr2_1, Sig_1, sigLOS2_1, vs1_1, vs2_1 = \
                    sigp_fit_vs(rbin1,rbin1,nupars1u,Mparsu,betpars1,Mstar)
                sigr2_2, Sig_2, sigLOS2_2, vs1_2, vs2_2 = \
                    sigp_fit_vs(rbin2,rbin2,nupars2u,Mparsu,betpars2,Mstar)                
            Mr = M(rbin,Mparsu)
            betar1 = beta(rbin1,betpars1)
            betar2 = beta(rbin2,betpars2)
            rhor = rho(rbin,Mparsu)
            dlnrhodlnrr = dlnrhodlnr(rbin,Mparsu)
            Mstarr = Mstar_prof*Mstar

            Mminmax[:,i] = Mr
            bet1minmax[:,i] = betar1
            bet2minmax[:,i] = betar2
            betstar1minmax[:,i] = betar1/(2.0-betar1)
            betstar2minmax[:,i] = betar2/(2.0-betar2)
            sigp1minmax[:,i] = np.sqrt(sigLOS2_1)/1000.
            sigp2minmax[:,i] = np.sqrt(sigLOS2_2)/1000.
            Sig1minmax[:,i] = Sig_1
            Sig2minmax[:,i] = Sig_2
            rhominmax[:,i] = rhor
            dlnrhodlnrminmax[:,i] = dlnrhodlnrr

    #And now make the plots:

    #Surface density first:
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)

    plt.loglog()

    if (singlesplit == 'single'):
#        plt.errorbar(rbin_phot,surfden,surfdenerr,\
#                     color='b',ecolor='b',linewidth=mylinewidth)

        if (whichdata == 'Fornax'):
            plt.errorbar(rbin_phot,surfden,surfdenerr,\
                         color='b',ecolor='b',linewidth=2,alpha=0.25,\
                         fmt='o')
        else:
            plt.errorbar(rbin_phot,surfden,surfdenerr,\
                         color='b',ecolor='b',linewidth=2,alpha=0.75,\
                         fmt='o')

        plt.fill_between(rbin,Sig_int[5,:],Sig_int[6,:],\
                             facecolor='black',alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(rbin,Sig_int[3,:],Sig_int[4,:],\
                             facecolor='black',alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(rbin,Sig_int[1,:],Sig_int[2,:],\
                             facecolor='black',alpha=0.66,\
                             edgecolor='none')
        plt.plot(rbin,Sig_int[0,:],'k',linewidth=mylinewidth,\
                     label=r'Fit')
        plt.plot([Rhalf,Rhalf],[1e-30,1e30],'b',alpha=0.5,\
                     linewidth=mylinewidth)
#        if (plotminmax == 'yes'):
#            plt.plot(rbin,Sigminmax[:,0],'g',linewidth=4)
#            plt.plot(rbin,Sigminmax[:,1],'g',linewidth=2)
    elif (singlesplit == 'split'):
        plt.errorbar(rbin1_phot,surfden1,surfdenerr1,\
                     color=colorpop1data,ecolor=colorpop1data,\
                     linewidth=mylinewidth)
        plt.fill_between(rbin1,Sig1_int[5,:],Sig1_int[6,:],\
                             facecolor=colorpop1,alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(rbin1,Sig1_int[3,:],Sig1_int[4,:],\
                             facecolor=colorpop1,alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(rbin1,Sig1_int[1,:],Sig1_int[2,:],\
                             facecolor=colorpop1,alpha=0.66,\
                             edgecolor='none')
        plt.plot(rbin1,Sig1_int[0,:],color=colorpop1,linewidth=mylinewidth,\
                     label=r'Fit 1')
        plt.errorbar(rbin2_phot,surfden2,surfdenerr2,\
                     color=colorpop2data,ecolor=colorpop2data,\
                     linewidth=mylinewidth)
        plt.fill_between(rbin2,Sig2_int[5,:],Sig2_int[6,:],\
                             facecolor=colorpop2,alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(rbin2,Sig2_int[3,:],Sig2_int[4,:],\
                             facecolor=colorpop2,alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(rbin2,Sig2_int[1,:],Sig2_int[2,:],\
                             facecolor=colorpop2,alpha=0.66,\
                             edgecolor='none')
        plt.plot(rbin2,Sig2_int[0,:],color=colorpop2,linewidth=mylinewidth,\
                     label=r'Fit 2')

#        if (plotminmax == 'yes'):
#            plt.plot(rbin1,Sig1minmax[:,0],'g',linewidth=4)
#            plt.plot(rbin1,Sig1minmax[:,1],'g',linewidth=2)
#            plt.plot(rbin2,Sig2minmax[:,0],'g',linewidth=4)
#            plt.plot(rbin2,Sig2minmax[:,1],'g',linewidth=2)

        plt.plot([Rhalf1,Rhalf1],[1e-30,1e30],color=colorpop1,alpha=0.5,\
                     linewidth=mylinewidth)
        plt.plot([Rhalf2,Rhalf2],[1e-30,1e30],color=colorpop2,alpha=0.5,\
                     linewidth=mylinewidth)   

    plt.xlabel(r'$R\,[{\rm kpc}]$',\
                   fontsize=myfontsize)
    plt.ylabel(r'$\Sigma_*\,[N\,{\rm kpc}^{-2}]$',\
                   fontsize=myfontsize)
    if (rbin_inner < 0):
        plt.xlim([np.min(rbin),np.max(rbin)])
    else:
        plt.xlim([rbin_inner,rbin_outer])
    plt.ylim([ymin_Sigstar,ymax_Sigstar])

    if (plotlegend == 'yes'):
        plt.legend(fontsize=mylegendfontsize)

    plt.savefig(outdir+'output_Sigstar.pdf',bbox_inches='tight')

    #And the projected velocity dispersion: 
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    if (reduce_xticks == 'yes'):
        plt.locator_params(axis='x',nbins=4)
    
    if (singlesplit == 'single'):
        plt.errorbar(np.log10(rbin_kin),sigpmean,sigperr,\
                     linewidth=2,color='b',alpha=0.75,\
                     fmt='o')

        sel = sigp_int[0,:] > 0
        plt.fill_between(np.log10(rbin[sel]),sigp_int[5,:][sel],\
                             sigp_int[6,:][sel],\
                             facecolor='black',alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin[sel]),sigp_int[3,:][sel],\
                             sigp_int[4,:][sel],\
                             facecolor='black',alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin[sel]),sigp_int[1,:][sel],\
                             sigp_int[2,:][sel],\
                             facecolor='black',alpha=0.66,\
                             edgecolor='none')
        plt.plot(np.log10(rbin[sel]),sigp_int[0,:][sel],'k',linewidth=mylinewidth,\
                     label=r'Fit')
        if (plotminmax == 'yes'):
            plt.plot(np.log10(rbin),sigpminmax[:,0],'g',linewidth=4)
            plt.plot(np.log10(rbin),sigpminmax[:,1],'g',linewidth=2)
        plt.plot([np.log10(Rhalf),np.log10(Rhalf)],[0,100],'b',alpha=0.5,\
                     linewidth=mylinewidth)
    elif (singlesplit == 'split'):
        plt.errorbar(np.log10(rbin1_kin),sigpmean1,sigperr1,\
                     color=colorpop1data,ecolor=colorpop1data,\
                     linewidth=mylinewidth)
        plt.errorbar(np.log10(rbin2_kin),sigpmean2,sigperr2,\
                     color=colorpop2data,ecolor=colorpop2data,\
                     linewidth=mylinewidth)
        plt.fill_between(np.log10(rbin1),sigp1_int[5,:],sigp1_int[6,:],\
                             facecolor=colorpop1,alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin1),sigp1_int[3,:],sigp1_int[4,:],\
                             facecolor=colorpop1,alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin1),sigp1_int[1,:],sigp1_int[2,:],\
                             facecolor=colorpop1,alpha=0.66,\
                             edgecolor='none')
        plt.plot(np.log10(rbin1),sigp1_int[0,:],color=colorpop1,\
                     linewidth=mylinewidth,\
                     label=r'Fit 1')
        plt.fill_between(np.log10(rbin2),sigp2_int[5,:],sigp2_int[6,:],\
                             facecolor=colorpop2,alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin2),sigp2_int[3,:],sigp2_int[4,:],\
                             facecolor=colorpop2,alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin2),sigp2_int[1,:],sigp2_int[2,:],\
                             facecolor=colorpop2,alpha=0.66,\
                             edgecolor='none')
        plt.plot(np.log10(rbin2),sigp2_int[0,:],color=colorpop2,\
                     linewidth=mylinewidth,\
                     label=r'Fit 2')

        if (plotminmax == 'yes'):
            plt.plot(np.log10(rbin1),sigp1minmax[:,0],'g',\
                     linewidth=4)
            plt.plot(np.log10(rbin1),sigp1minmax[:,1],'g',\
                     linewidth=2)
            plt.plot(np.log10(rbin2),sigp2minmax[:,0],'g',\
                     linewidth=4)
            plt.plot(np.log10(rbin2),sigp2minmax[:,1],'g',\
                     linewidth=2)

        plt.plot([np.log10(Rhalf1),np.log10(Rhalf1)],[0,1e4],color=colorpop1,\
                     alpha=0.5,\
                     linewidth=mylinewidth)
        plt.plot([np.log10(Rhalf2),np.log10(Rhalf2)],[0,1e4],color=colorpop2,\
                     alpha=0.5,\
                     linewidth=mylinewidth)

    plt.xlabel(r'${\rm Log}_{10}[R/{\rm kpc}]$',\
                   fontsize=myfontsize)
    plt.ylabel(r'$\sigma_{\rm LOS}[{\rm km\,s}^{-1}]$',\
                   fontsize=myfontsize)
    plt.ylim([0,y_sigLOSmax])
    if (rbin_inner < 0):
        plt.xlim([np.log10(np.min(rbin)),np.log10(np.max(rbin))])
    else:
        plt.xlim([np.log10(rbin_inner),np.log10(rbin_outer)])

    if (plotlegend == 'yes'):
        plt.legend(fontsize=mylegendfontsize)
    
    plt.savefig(outdir+'output_sigLOS.pdf',bbox_inches='tight')

    #And the proper motion projected dispersions (if using proper motions): 
    if (propermotion == 'yes'):
        #First in the radial direction (on the sky):
        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)
        if (reduce_xticks == 'yes'):
            plt.locator_params(axis='x',nbins=4)

        if (singlesplit == 'single'):
            psel = sigpmr > 0
            plt.errorbar(np.log10(rbin_kin[psel]),sigpmr[psel],sigpmrerr[psel],\
                         linewidth=2,color='b',alpha=0.75,\
                         fmt='o')
            
            sel = sigpmr_int[0,:] > 0
            plt.fill_between(np.log10(rbin[sel]),sigpmr_int[5,:][sel],\
                                 sigpmr_int[6,:][sel],\
                                 facecolor='black',alpha=alp3sig,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin[sel]),sigpmr_int[3,:][sel],\
                                 sigpmr_int[4,:][sel],\
                                 facecolor='black',alpha=0.33,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin[sel]),sigpmr_int[1,:][sel],\
                                 sigpmr_int[2,:][sel],\
                                 facecolor='black',alpha=0.66,\
                                 edgecolor='none')
            plt.plot(np.log10(rbin[sel]),sigpmr_int[0,:][sel],'k',linewidth=mylinewidth,\
                         label=r'Fit')
            if (plotminmax == 'yes'):
                plt.plot(np.log10(rbin),sigpmrminmax[:,0],'g',\
                         linewidth=4)
                plt.plot(np.log10(rbin),sigpmrminmax[:,1],'g',\
                         linewidth=2)
            plt.plot([np.log10(Rhalf),np.log10(Rhalf)],[0,100],'b',\
                         alpha=0.5,\
                         linewidth=mylinewidth)
        elif (singlesplit == 'split'):
            plt.errorbar(np.log10(rbin1_kin),sigpmr1,sigpmr1,\
                             color=colorpop1data,ecolor=colorpop1data,\
                             linewidth=mylinewidth)
            plt.errorbar(np.log10(rbin2_kin),sigpmr2,sigpmr2,\
                             color=colorpop2data,ecolor=colorpop2data,\
                             linewidth=mylinewidth)
            plt.fill_between(np.log10(rbin1),sigpmr1_int[5,:],\
                                 sigpmr1_int[6,:],\
                                 facecolor=colorpop1,alpha=alp3sig,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin1),sigpmr1_int[3,:],\
                                 sigpmr1_int[4,:],\
                                 facecolor=colorpop1,alpha=0.33,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin1),sigpmr1_int[1,:],\
                                 sigpmr1_int[2,:],\
                                 facecolor=colorpop1,alpha=0.66,\
                                 edgecolor='none')
            plt.plot(np.log10(rbin1),sigpmr1_int[0,:],colorpop1,\
                     linewidth=mylinewidth,\
                     label=r'Fit 1')
            plt.fill_between(np.log10(rbin2),sigpmr2_int[5,:],\
                                 sigpmr2_int[6,:],\
                                 facecolor=colorpop2,alpha=alp3sig,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin2),sigpmr2_int[3,:],\
                                 sigpmr2_int[4,:],\
                                 facecolor=colorpop2,alpha=0.33,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin2),sigpmr2_int[1,:],\
                                 sigpmr2_int[2,:],\
                                 facecolor=colorpop2,alpha=0.66,\
                                 edgecolor='none')
            plt.plot(np.log10(rbin2),sigpmr2_int[0,:],color=colorpop2,\
                         linewidth=mylinewidth,\
                         label=r'Fit 2')
            plt.plot([np.log10(Rhalf1),np.log10(Rhalf1)],[0,100],'b',\
                         alpha=0.5,\
                         linewidth=mylinewidth)

            if (plotminmax == 'yes'):
                plt.plot(np.log10(rbin1),sigpmr1minmax[:,0],'g',\
                         linewidth=4)
                plt.plot(np.log10(rbin1),sigpmr1minmax[:,1],'g',\
                         linewidth=2)
                plt.plot(np.log10(rbin2),sigpmr2minmax[:,0],'g',\
                         linewidth=4)
                plt.plot(np.log10(rbin2),sigpmr2minmax[:,1],'g',\
                         linewidth=2)

            plt.plot([np.log10(Rhalf2),np.log10(Rhalf2)],[0,100],\
                         color=colorpop2,\
                         alpha=0.5,\
                         linewidth=mylinewidth)

        plt.xlabel(r'${\rm Log}_{10}[R/{\rm kpc}]$',\
                       fontsize=myfontsize)
        plt.ylabel(r'$\sigma_{\rm pmr}[{\rm km\,s}^{-1}]$',\
                       fontsize=myfontsize)
        plt.ylim([0,y_sigLOSmax])

        if (rbin_inner < 0):
            plt.xlim([np.log10(np.min(rbin)),np.log10(np.max(rbin))])
        else:
            plt.xlim([np.log10(rbin_inner),np.log10(rbin_outer)])

        if (plotlegend == 'yes'):
            plt.legend(fontsize=mylegendfontsize)

        plt.savefig(outdir+'output_sigpmr.pdf',bbox_inches='tight')

        #Then in the tangential direction (on the sky): 
        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)
        if (reduce_xticks == 'yes'):
            plt.locator_params(axis='x',nbins=4)

        if (singlesplit == 'single'):
            psel = sigpmt > 0
            plt.errorbar(np.log10(rbin_kin[psel]),sigpmt[psel],sigpmterr[psel],\
                         linewidth=2,color='b',alpha=0.75,\
                         fmt='o')
            
            sel = sigpmt_int[0,:] > 0
            plt.fill_between(np.log10(rbin[sel]),sigpmt_int[5,:][sel],\
                                 sigpmt_int[6,:][sel],\
                                 facecolor='black',alpha=alp3sig,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin[sel]),sigpmt_int[3,:][sel],\
                                 sigpmt_int[4,:][sel],\
                                 facecolor='black',alpha=0.33,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin[sel]),sigpmt_int[1,:][sel],\
                                 sigpmt_int[2,:][sel],\
                                 facecolor='black',alpha=0.66,\
                                 edgecolor='none')
            plt.plot(np.log10(rbin[sel]),sigpmt_int[0,:][sel],'k',linewidth=mylinewidth,\
                         label=r'Fit')
            if (plotminmax == 'yes'):
                plt.plot(np.log10(rbin),sigpmtminmax[:,0],'g',\
                         linewidth=4)
                plt.plot(np.log10(rbin),sigpmtminmax[:,1],'g',\
                         linewidth=2)
            plt.plot([np.log10(Rhalf),np.log10(Rhalf)],[0,100],'b',\
                         alpha=0.5,\
                         linewidth=mylinewidth)

        elif (singlesplit == 'split'):
            plt.errorbar(np.log10(rbin1_kin),sigpmt1,sigpmt1,\
                             color=colorpop1data,ecolor=colorpop1data,\
                             linewidth=mylinewidth)
            plt.errorbar(np.log10(rbin2_kin),sigpmt2,sigpmt2,\
                             color=colorpop2data,ecolor=colorpop2data,\
                             linewidth=mylinewidth)
            plt.fill_between(np.log10(rbin1),sigpmt1_int[5,:],\
                                 sigpmt1_int[6,:],\
                                 facecolor=colorpop1,alpha=alp3sig,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin1),sigpmt1_int[3,:],\
                                 sigpmt1_int[4,:],\
                                 facecolor=colorpop1,alpha=0.33,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin1),sigpmt1_int[1,:],\
                                 sigpmt1_int[2,:],\
                                 facecolor=colorpop1,alpha=0.66,\
                                 edgecolor='none')
            plt.plot(np.log10(rbin1),sigpmt1_int[0,:],colorpop1,\
                     linewidth=mylinewidth,\
                     label=r'Fit 1')
            plt.fill_between(np.log10(rbin2),sigpmt2_int[5,:],\
                                 sigpmt2_int[6,:],\
                                 facecolor=colorpop2,alpha=alp3sig,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin2),sigpmt2_int[3,:],\
                                 sigpmt2_int[4,:],\
                                 facecolor=colorpop2,alpha=0.33,\
                                 edgecolor='none')
            plt.fill_between(np.log10(rbin2),sigpmt2_int[1,:],\
                                 sigpmt2_int[2,:],\
                                 facecolor=colorpop2,alpha=0.66,\
                                 edgecolor='none')
            plt.plot(np.log10(rbin2),sigpmt2_int[0,:],color=colorpop2,\
                         linewidth=mylinewidth,\
                         label=r'Fit 2')

            if (plotminmax == 'yes'):
                plt.plot(np.log10(rbin1),sigpmt1minmax[:,0],'g',\
                         linewidth=4)
                plt.plot(np.log10(rbin1),sigpmt1minmax[:,1],'g',\
                         linewidth=2)
                plt.plot(np.log10(rbin2),sigpmt2minmax[:,0],'g',\
                         linewidth=4)
                plt.plot(np.log10(rbin2),sigpmt2minmax[:,1],'g',\
                         linewidth=2)

            plt.plot([np.log10(Rhalf1),np.log10(Rhalf1)],[0,100],\
                         color=colorpop1,\
                         alpha=0.5,\
                         linewidth=mylinewidth)
            plt.plot([np.log10(Rhalf2),np.log10(Rhalf2)],[0,100],\
                         color=colorpop2,\
                         alpha=0.5,\
                         linewidth=mylinewidth)
        plt.xlabel(r'${\rm Log}_{10}[R/{\rm kpc}]$',\
                       fontsize=myfontsize)
        plt.ylabel(r'$\sigma_{\rm pmt}[{\rm km\,s}^{-1}]$',\
                       fontsize=myfontsize)
        plt.ylim([0,y_sigLOSmax])

        if (rbin_inner < 0):
            plt.xlim([np.log10(np.min(rbin)),np.log10(np.max(rbin))])
        else:
            plt.xlim([np.log10(rbin_inner),np.log10(rbin_outer)])

        if (plotlegend == 'yes'):
            plt.legend(fontsize=mylegendfontsize)

        plt.savefig(outdir+'output_sigpmt.pdf',bbox_inches='tight')

    #And now the beta profile:
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    if (reduce_xticks == 'yes'):
        plt.locator_params(axis='x',nbins=4)

    if (singlesplit == 'single'):
        plt.fill_between(np.log10(rbin),bet_int[5,:],bet_int[6,:],\
                             facecolor='black',alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin),bet_int[3,:],bet_int[4,:],\
                             facecolor='black',alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin),bet_int[1,:],bet_int[2,:],\
                             facecolor='black',alpha=0.66,\
                             edgecolor='none')
        plt.plot(np.log10(rbin),bet_int[0,:],'k',linewidth=mylinewidth,\
                     label=r'Fit')
        if (plotminmax == 'yes'):
            plt.plot(np.log10(rbin),betminmax[:,0],'g',linewidth=4)
            plt.plot(np.log10(rbin),betminmax[:,1],'g',linewidth=2)
        plt.plot([np.log10(Rhalf),np.log10(Rhalf)],[-100,100],'b',alpha=0.5,\
                     linewidth=mylinewidth)
    elif (singlesplit == 'split'):
        plt.fill_between(np.log10(rbin1),bet1_int[5,:],bet1_int[6,:],\
                             facecolor=colorpop1,alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin1),bet1_int[3,:],bet1_int[4,:],\
                             facecolor=colorpop1,alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin1),bet1_int[1,:],bet1_int[2,:],\
                             facecolor=colorpop1,alpha=0.66,\
                             edgecolor='none')
        plt.plot(np.log10(rbin1),bet1_int[0,:],colorpop1,linewidth=mylinewidth,\
                     label=r'Fit 1')
        plt.fill_between(np.log10(rbin2),bet2_int[5,:],bet2_int[6,:],\
                             facecolor=colorpop2,alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin2),bet2_int[3,:],bet2_int[4,:],\
                             facecolor=colorpop2,alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin2),bet2_int[1,:],bet2_int[2,:],\
                             facecolor=colorpop2,alpha=0.66,\
                             edgecolor='none')
        plt.plot(np.log10(rbin2),bet2_int[0,:],color=colorpop2,\
                     linewidth=mylinewidth,\
                     label=r'Fit 2')
        if (plotminmax == 'yes'):
            plt.plot(np.log10(rbin1),bet1minmax[:,0],'g',\
                     linewidth=4)
            plt.plot(np.log10(rbin1),bet1minmax[:,1],'g',\
                     linewidth=2)
            plt.plot(np.log10(rbin2),bet2minmax[:,0],'g',\
                     linewidth=4)
            plt.plot(np.log10(rbin2),bet2minmax[:,1],'g',\
                     linewidth=2)
        plt.plot([np.log10(Rhalf1),np.log10(Rhalf1)],[-100,100],\
                  color=colorpop1,\
                  alpha=0.5,linewidth=mylinewidth)
        plt.plot([np.log10(Rhalf2),np.log10(Rhalf2)],[-100,100],\
                  color=colorpop2,\
                  alpha=0.5,linewidth=mylinewidth)

    #And true answer (mock data):
    if (overtrue == 'yes'):
        if (singlesplit == 'single'):
            if (betprofile == 'osipkov'):
                betatrue = ranal**2./(ranal**2. + ra**2.)
            elif (betprofile == 'triaxial'):
                betatrue = (bet_rs**bet_eta*bet_bet0+ranal**bet_eta*bet_betinf)/\
                    (ranal**bet_eta+bet_rs**bet_eta)
            elif (betprofile == 'tangential'):
                betatrue = np.zeros(len(ranal)) + ra
            plt.plot(np.log10(ranal),betatrue,'b--',linewidth=mylinewidth,\
                         label=r'True')
        elif (singlesplit == 'split'):
            betatrue = ranal**2./(ranal**2. + ra1**2.)
            plt.plot(np.log10(ranal),betatrue,'b--',linewidth=mylinewidth,\
                         label=r'True 1')
            betatrue = ranal**2./(ranal**2. + ra2**2.)
            plt.plot(np.log10(ranal),betatrue,'r--',linewidth=mylinewidth,\
                         label=r'True 2')
    
    plt.xlabel(r'${\rm Log}_{10}[r/{\rm kpc}]$',\
                   fontsize=myfontsize)
    plt.ylabel(r'$\beta$',\
                   fontsize=myfontsize)
    plt.ylim([np.min([bet0min,betinfmin]),np.max([bet0max,betinfmax])])
    if (rbin_inner < 0):
        plt.xlim([np.log10(np.min(rbin)),np.log10(np.max(rbin))])
    else:
        plt.xlim([np.log10(rbin_inner),np.log10(rbin_outer)])

    if (plotlegend == 'yes'):
        plt.legend(fontsize=mylegendfontsize,loc='upper left')
    plt.savefig(outdir+'output_beta.pdf',bbox_inches='tight')

    #Write the above data to a file for comparitive plotting later:
    if (singlesplit == 'single'):
        f = open(outdir+'output_bet.txt','w')
        for i in range(len(rbin)):
            f.write('%f %f %f %f %f %f %f %f\n' % \
                        (rbin[i],bet_int[0,i],bet_int[1,i],bet_int[2,i],bet_int[3,i],\
                             bet_int[4,i],bet_int[5,i],bet_int[6,i]))
        f.close()
    else:
        f = open(outdir+'output_bet1.txt','w')
        for i in range(len(rbin)):
            f.write('%f %f %f %f %f %f %f %f\n' % \
                        (rbin[i],bet1_int[0,i],bet1_int[1,i],bet1_int[2,i],bet1_int[3,i],\
                             bet1_int[4,i],bet1_int[5,i],bet1_int[6,i]))
        f.close()
        f = open(outdir+'output_bet2.txt','w')
        for i in range(len(rbin2)):
            f.write('%f %f %f %f %f %f %f %f\n' % \
                        (rbin2[i],bet2_int[0,i],bet2_int[1,i],bet2_int[2,i],bet2_int[3,i],\
                             bet2_int[4,i],bet2_int[5,i],bet2_int[6,i]))
        f.close()
    if (plotminmax == 'yes'):
        f = open(outdir+'output_betminmax.txt','w')
        for i in range(len(rbin)):
            f.write('%f %f %f\n' % \
                        (rbin[i],betminmax[i,0],betminmax[i,1]))
        f.close()

    #And the *symmetrised* betastar:
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    if (reduce_xticks == 'yes'):
        plt.locator_params(axis='x',nbins=4)

    if (singlesplit == 'single'):
        plt.fill_between(np.log10(rbin),betstar_int[5,:],betstar_int[6,:],\
                             facecolor='black',alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin),betstar_int[3,:],betstar_int[4,:],\
                             facecolor='black',alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin),betstar_int[1,:],betstar_int[2,:],\
                             facecolor='black',alpha=0.66,\
                             edgecolor='none')
        plt.plot(np.log10(rbin),betstar_int[0,:],'k',linewidth=mylinewidth,\
                     label=r'Fit')
        if (plotminmax == 'yes'):
            plt.plot(np.log10(rbin),betstarminmax[:,0],'g',\
                     linewidth=4)
            plt.plot(np.log10(rbin),betstarminmax[:,1],'g',\
                     linewidth=2)
        plt.plot([np.log10(Rhalf),np.log10(Rhalf)],[-100,100],'b',alpha=0.5,\
                     linewidth=mylinewidth)
    elif (singlesplit == 'split'):
        plt.fill_between(np.log10(rbin1),betstar1_int[5,:],betstar1_int[6,:],\
                             facecolor=colorpop1,alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin1),betstar1_int[3,:],betstar1_int[4,:],\
                             facecolor=colorpop1,alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin1),betstar1_int[1,:],betstar1_int[2,:],\
                             facecolor=colorpop1,alpha=0.66,\
                             edgecolor='none')
        plt.plot(np.log10(rbin1),betstar1_int[0,:],color=colorpop1,\
                     linewidth=mylinewidth,\
                     label=r'Fit 1')
        plt.fill_between(np.log10(rbin2),betstar2_int[5,:],betstar2_int[6,:],\
                             facecolor=colorpop2,alpha=alp3sig,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin2),betstar2_int[3,:],betstar2_int[4,:],\
                             facecolor=colorpop2,alpha=0.33,\
                             edgecolor='none')
        plt.fill_between(np.log10(rbin2),betstar2_int[1,:],betstar2_int[2,:],\
                             facecolor=colorpop2,alpha=0.66,\
                             edgecolor='none')
        plt.plot(np.log10(rbin2),betstar2_int[0,:],colorpop2,linewidth=mylinewidth,\
                     label=r'Fit 2')
        if (plotminmax == 'yes'):
            plt.plot(np.log10(rbin1),betstar1minmax[:,0],'g',\
                     linewidth=4)
            plt.plot(np.log10(rbin1),betstar1minmax[:,1],'g',\
                     linewidth=2)
            plt.plot(np.log10(rbin2),betstar2minmax[:,0],'g',\
                     linewidth=4)
            plt.plot(np.log10(rbin2),betstar2minmax[:,1],'g',\
                     linewidth=2)
        plt.plot([np.log10(Rhalf1),np.log10(Rhalf1)],[-100,100],\
                     color=colorpop1,alpha=0.5,\
                     linewidth=mylinewidth)
        plt.plot([np.log10(Rhalf2),np.log10(Rhalf2)],[-100,100],\
                     color=colorpop2,alpha=0.5,\
                     linewidth=mylinewidth)

    #And true answer (mock data):
    if (overtrue == 'yes'):
        if (singlesplit == 'single'):
            if (betprofile == 'osipkov'):
                betatrue = ranal**2./(ranal**2. + ra**2.)
            elif (betprofile == 'triaxial'):
                betatrue = (bet_rs**bet_eta*bet_bet0+ranal**bet_eta*bet_betinf)/\
                    (ranal**bet_eta+bet_rs**bet_eta)
            elif (betprofile == 'tangential'):
                betatrue = np.zeros(len(ranal)) + ra
            betatruestar = betatrue/(2.0-betatrue)
            plt.plot(np.log10(ranal),betatruestar,'b--',linewidth=mylinewidth,\
                         label=r'True')
        elif (singlesplit == 'split'):
            betatrue = ranal**2./(ranal**2. + ra1**2.)
            betatruestar = betatrue/(2.0-betatrue)
            plt.plot(np.log10(ranal),betatruestar,'b--',linewidth=mylinewidth,\
                         label=r'True 1')
            betatrue = ranal**2./(ranal**2. + ra2**2.)
            betatruestar = betatrue/(2.0-betatrue)
            plt.plot(np.log10(ranal),betatruestar,'r--',linewidth=mylinewidth,\
                         label=r'True 2')

    plt.xlabel(r'${\rm Log}_{10}[r/{\rm kpc}]$',\
                   fontsize=myfontsize)
    plt.ylabel(r'$\tilde{\beta}$',\
                   fontsize=myfontsize)
    plt.ylim([-1,1])
    if (rbin_inner < 0):
        plt.xlim([np.log10(np.min(rbin)),np.log10(np.max(rbin))])
    else:
        plt.xlim([np.log10(rbin_inner),np.log10(rbin_outer)])

    if (plotlegend == 'yes'):
        plt.legend(fontsize=mylegendfontsize,loc='upper left')
    plt.savefig(outdir+'output_betastar.pdf',bbox_inches='tight')

    #Write the above data to a file for comparitive plotting later:
    if (singlesplit == 'single'):
        f = open(outdir+'output_betstar.txt','w')
        for i in range(len(rbin)):
            f.write('%f %f %f %f %f %f %f %f\n' % \
                        (rbin[i],betstar_int[0,i],betstar_int[1,i],\
                             betstar_int[2,i],betstar_int[3,i],\
                             betstar_int[4,i],betstar_int[5,i],\
                             betstar_int[6,i]))
        f.close()
    else:
        f = open(outdir+'output_betstar1.txt','w')
        for i in range(len(rbin)):
            f.write('%f %f %f %f %f %f %f %f\n' % \
                        (rbin[i],betstar1_int[0,i],betstar1_int[1,i],\
                             betstar1_int[2,i],betstar1_int[3,i],\
                             betstar1_int[4,i],betstar1_int[5,i],\
                             betstar1_int[6,i]))
        f.close()
        f = open(outdir+'output_betstar2.txt','w')
        for i in range(len(rbin2)):
            f.write('%f %f %f %f %f %f %f %f\n' % \
                        (rbin2[i],betstar2_int[0,i],betstar2_int[1,i],\
                             betstar2_int[2,i],betstar2_int[3,i],\
                             betstar2_int[4,i],betstar2_int[5,i],\
                             betstar2_int[6,i]))
        f.close()
    if (plotminmax == 'yes'):
        f = open(outdir+'output_betstarminmax.txt','w')
        for i in range(len(rbin)):
            f.write('%f %f %f\n' % \
                        (rbin[i],betstarminmax[i,0],betstarminmax[i,1]))
        f.close()

    #And the mass profile: 
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    plt.loglog()

    plt.fill_between(rbin,M_int[5,:],M_int[6,:],\
                         facecolor='black',alpha=alp3sig,\
                         edgecolor='none')
    plt.fill_between(rbin,M_int[3,:],M_int[4,:],\
                         facecolor='black',alpha=0.33,\
                         edgecolor='none')
    plt.fill_between(rbin,M_int[1,:],M_int[2,:],\
                         facecolor='black',alpha=0.66,\
                         edgecolor='none')
    if (Mstar > 1.0):
        plt.plot(rbin,M_int[0,:],'k',linewidth=mylinewidth,\
                     label=r'Fit Dark Matter')
    else:
        plt.plot(rbin,M_int[0,:],'k',linewidth=mylinewidth,\
                     label=r'Fit')

    plt.fill_between(Mstar_rad,Mstar_int[5,:],Mstar_int[6,:],\
                         facecolor=colorpop2,alpha=alp3sig,\
                         edgecolor='none')
    plt.fill_between(Mstar_rad,Mstar_int[3,:],Mstar_int[4,:],\
                         facecolor=colorpop2,alpha=0.33,\
                         edgecolor='none')
    plt.fill_between(Mstar_rad,Mstar_int[1,:],Mstar_int[2,:],\
                         facecolor=colorpop2,alpha=0.66,\
                         edgecolor='none')
    if (Mstar > 1.0):
        plt.plot(Mstar_rad,Mstar_int[0,:],color=colorpop2,linewidth=mylinewidth,\
                     label=r'Fit Stars')

    if (plotminmax == 'yes'):
        plt.plot(rbin,Mminmax[:,0],'g',linewidth=4,\
                 label=r'minmax$[\beta_*]$')
        plt.plot(rbin,Mminmax[:,1],'g',linewidth=2)

    #Overplot true answer (mock data):
    if (overtrue == 'yes'):
        truemass = alpbetgammass(ranal,rho0,r0,alp,bet,gam)
        plt.plot(ranal,truemass,'b--',linewidth=mylinewidth,\
                 label=r'True')        

    if (singlesplit == 'single'):
        plt.plot([Rhalf,Rhalf],[1e-30,1e30],'b',alpha=0.5,\
                     linewidth=mylinewidth)
    elif (singlesplit == 'split'):
        plt.plot([Rhalf1,Rhalf1],[1e-30,1e30],color=colorpop1,\
                     alpha=0.5,\
                     linewidth=mylinewidth)
        plt.plot([Rhalf2,Rhalf2],[1e-30,1e30],color=colorpop2,\
                     alpha=0.5,\
                     linewidth=mylinewidth)

    plt.xlabel(r'$r\,[{\rm kpc}]$',\
                   fontsize=myfontsize)
    plt.ylabel(r'$M(<r)\,[{\rm M}_\odot]$',\
                   fontsize=myfontsize)

    plt.ylim([yMlow,yMhigh])

    if (rbin_inner < 0):
        plt.xlim([np.min(rbin),np.max(rbin)])
    else:
        plt.xlim([rbin_inner,rbin_outer])

    if (plotlegend_mass == 'yes'):
        plt.legend(fontsize=mylegendfontsize,loc='lower right')
    plt.savefig(outdir+'output_Mass.pdf',bbox_inches='tight')

    #Write the above data to a file for comparitive plotting later:
    f = open(outdir+'output_M.txt','w')
    for i in range(len(rbin)):
        f.write('%f %f %f %f %f %f %f %f\n' % \
         (rbin[i],M_int[0,i],M_int[1,i],M_int[2,i],M_int[3,i],\
              M_int[4,i], M_int[5,i], M_int[6,i]))
    f.close()

    #And the Mdyn/Mstar ratio:
    f = open(outdir+'output_MdynMstar.txt','w')
    for i in range(len(Mstar_rad)):
        f.write('%f %f %f %f %f %f %f %f\n' % \
         (Mstar_rad[i],Mdynrat_int[0,i],Mdynrat_int[1,i],Mdynrat_int[2,i],\
              Mdynrat_int[3,i],\
              Mdynrat_int[4,i], Mdynrat_int[5,i], Mdynrat_int[6,i]))
    f.close()

    #And nu_mass_r:
    f = open(outdir+'output_nu_mass_r.txt','w')
    for i in range(len(Mstar_rad)):
        f.write('%f %f %f %f %f %f %f %f\n' % \
         (Mstar_rad[i],nu_int[0,i],nu_int[1,i],nu_int[2,i],\
              nu_int[3,i],\
              nu_int[4,i], nu_int[5,i], nu_int[6,i]))
    f.close()

    if (plotminmax == 'yes'):
        f = open(outdir+'output_Mminmax.txt','w')
        for i in range(len(rbin)):
            f.write('%f %f %f\n' % \
                        (rbin[i],Mminmax[i,0],Mminmax[i,1]))
        f.close()

    #And the density profile:
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    plt.loglog()

    if (plotdenresults == 'yes'):
        plt.fill_between(rbin,rho_int[5,:],rho_int[6,:],\
                         facecolor='black',alpha=alp3sig,\
                         edgecolor='none')
        plt.fill_between(rbin,rho_int[3,:],rho_int[4,:],\
                         facecolor='black',alpha=0.33,\
                         edgecolor='none')
        plt.fill_between(rbin,rho_int[1,:],rho_int[2,:],\
                         facecolor='black',alpha=0.66,\
                         edgecolor='none')
        plt.plot(rbin,rho_int[0,:],'k',linewidth=mylinewidth,\
                     label=r'Fit')

    if (plotminmax == 'yes'):
        plt.plot(rbin,rhominmax[:,0],'g',linewidth=4)
        plt.plot(rbin,rhominmax[:,1],'g',linewidth=2)

    #Overplot true solution (for mock data): 
    if (overtrue == 'yes'):
        trueden = alpbetgamden(ranal,rho0,r0,alp,bet,gam)
        plt.plot(ranal,trueden,'b--',linewidth=mylinewidth,\
                 label=r'True')

    if (singlesplit == 'single'):
        plt.plot([Rhalf,Rhalf],[1e-30,1e30],'b',alpha=0.5,\
                     linewidth=mylinewidth)
    elif (singlesplit == 'split'):
        plt.plot([Rhalf1,Rhalf1],[1e-30,1e30],color=colorpop1,\
                     alpha=0.5,\
                     linewidth=mylinewidth)
        plt.plot([Rhalf2,Rhalf2],[1e-30,1e30],color=colorpop2,\
                     alpha=0.5,\
                     linewidth=mylinewidth)

    plt.xlabel(r'$r\,[{\rm kpc}]$',\
                   fontsize=myfontsize)
    plt.ylabel(r'$\rho\,[{\rm M}_\odot\,{\rm kpc}^{-3}]$',\
                   fontsize=myfontsize)

    if (rbin_inner < 0):
        plt.xlim([np.min(rbin),np.max(rbin)])
    else:
        plt.xlim([rbin_inner,rbin_outer])
    plt.ylim([1e4,5e9])

    #Overlay coreNFW:
    plot_evolve = 'no'
    plot_nfw = 'no'
    plot_evolve_low = 'no'
    plot_fullcore = 'no'
    plot_4sig = 'no'
    if (data_file_type == 'walker_style' or \
        data_file_type == 'LeoT' or \
        data_file_type == 'walker_style_vs' or \
        data_file_type == 'sluggs'):
        rhocrit = 135.05
        eta = 1.75
        kappa = 0.04
        oden = 200
        h = 0.7
        if (whichdata == 'Fornax'):
            M200 = [15.46e9]
            tSF = 14.0
            compare_file = \
                './Data/for_dsphmassespar6louie_massprofile.res'
            evolve_file = \
                './Data/Fornax_den_evolved_2sig.txt'
            evolve_file_low = \
                './Data/Fornax_den_evolved_2sig_lowrez.txt'
            evolve_file_nfw = \
                './Data/Fornax_den_evolved_2sig_NFW.txt'
            plot_evolve = 'yes'
            plot_nfw = 'yes'
            plot_evolve_low = 'no'
        elif (whichdata == 'Draco'):
            M200 = [1.31e9]
            tSF = 3.8
            compare_file = \
                './Data/dra_hecto_dsphmassespar6louie_massprofile.res'
            evolve_file = \
                './Data/Draco_den_evolved_median.txt'
            evolve_file_low = \
                './Data/Draco_den_evolved_median_lowrez.txt'
            evolve_file_fullcore = \
                './Data/Draco_den_evolved_median_fullcoreNFW.txt'
            plot_evolve = 'yes'
            plot_fullcore = 'no'
            plot_evolve_low = 'no'
        elif (whichdata == 'UMi'):
            M200 = [1.125e9]
            tSF = 4.8
        elif (whichdata == 'Carina'):
            M200 = [0.74e9]
            tSF = 11.8
            evolve_file = \
                './Data/Carina_den_evolved_median.txt'
            plot_evolve = 'yes'
        elif (whichdata == 'LeoI'):
            M200 = [4.13e9]
            tSF = 12.8
            evolve_file = \
                './Data/LeoI_den_evolved_median.txt'
            plot_evolve = 'yes'
        elif (whichdata == 'LeoII'):
            M200 = [1.52e9]
            tSF = 7.8
            evolve_file = \
                './Data/LeoII_den_evolved_median.txt'
            plot_evolve = 'yes'
        elif (whichdata == 'Sextans'):
            M200 = [1.18e9]
            tSF = 6.8
        elif (whichdata == 'Sculptor'):
            M200 = [5.19e9]
            tSF = 3.8
        elif (whichdata == 'LeoT'):
            M200 = [1e8,5e8,1e9,5e9,1e10]
            tSF = 14.0
        elif (whichdata == 'DracoMock' or whichdata == 'DracoMockSplit'):
            M200 = [1.31e9]
            tSF = 3.8
            evolve_file = \
                './Data/Draco_den_evolved_median_split_lowrez_14Gyrs.txt'
            plot_evolve = 'yes'
        elif (whichdata == 'Crater'):
            M200 = [1.0e6,1.0e7,1.0e8]
            tSF = 0.01
        elif (whichdata == 'NGC1407'):
            M200 = [1.0e13,5.0e13,1.0e14]
            tSF = 0.01
 
        for ii in range(len(M200)):
            c = cosmo_cfunc(M200[ii],h)
            logc_scat = 0.1
            clowp = 10.0**(np.log10(c)-logc_scat)
            chighp = 10.0**(np.log10(c)+logc_scat)
            rhoanal = \
                corenfw_den(ranal,M200[ii],c,oden,tSF,Rhalf,\
                            rhocrit,eta,kappa)
            rhoanal_NFW = \
                corenfw_den(ranal,M200[ii],c,oden,0.01,Rhalf,\
                            rhocrit,eta,kappa)
            rhoanal_low = \
                corenfw_den(ranal,M200[ii],clowp,oden,\
                            tSF,Rhalf,rhocrit,eta,kappa)
            rhoanal_NFW_low = \
                corenfw_den(ranal,M200[ii],clowp,oden,\
                            0.01,Rhalf,rhocrit,eta,kappa)
            rhoanal_high = \
                corenfw_den(ranal,M200[ii],chighp,oden,\
                            tSF,Rhalf,rhocrit,eta,kappa)
            rhoanal_NFW_high = corenfw_den(ranal,M200[ii],chighp,oden,\
                                               0.01,Rhalf,rhocrit,eta,kappa)
            if (ii == 0):
                plt.plot(ranal,rhoanal,color='red',linewidth=mylinewidth,\
                             label=r'coreNFW')
                plt.plot(ranal,rhoanal_NFW,'--',\
                             color='red',linewidth=mylinewidth,\
                             label=r'NFW')
                plt.fill_between(ranal,rhoanal_low,rhoanal_high,\
                                     facecolor='red',alpha=0.33,\
                                     edgecolor='none')
                plt.fill_between(ranal,rhoanal_NFW_low,rhoanal_NFW_high,\
                                     facecolor='red',alpha=0.33,\
                                     edgecolor='none')
            else:
                plt.plot(ranal,rhoanal,color='red',linewidth=mylinewidth)
                plt.plot(ranal,rhoanal_NFW,'--',\
                             color='red',linewidth=mylinewidth)
                plt.fill_between(ranal,rhoanal_low,rhoanal_high,\
                                     facecolor='red',alpha=0.33,\
                                     edgecolor='none')
                plt.fill_between(ranal,rhoanal_NFW_low,rhoanal_NFW_high,\
                                     facecolor='red',alpha=0.33,\
                                     edgecolor='none')
            
        #Compare with Walker:
        if (walker_compare == 'yes'):
            f = open(compare_file,'r')
            data = np.genfromtxt(f)
            f.close()
            plt.plot(data[:,0]/1000.,data[:,8]*1000.**3.0,\
                     linewidth=mylinewidth,\
                     color='blue',label=r'Walker')
            plt.plot(data[:,0]/1000.,data[:,6]*1000.**3.0,'--',\
                     linewidth=mylinewidth,\
                     color='blue',label=r'Walker')
            plt.plot(data[:,0]/1000.,data[:,10]*1000.**3.0,'--',\
                     linewidth=mylinewidth,\
                     color='blue',label=r'Walker')

        if (plot_evolve == 'yes'):
            f = open(evolve_file,'r')
            data = np.genfromtxt(f)
            f.close()
            plt.plot(data[:,0],data[:,1],\
                     linewidth=mylinewidth,\
                     color='blue',label=r'coreNFW+tides')
        if (plot_fullcore == 'yes'):
            f = open(evolve_file_fullcore,'r')
            data = np.genfromtxt(f)
            f.close()
            plt.plot(data[:,0],data[:,1],'--',\
                     linewidth=mylinewidth,\
                     color='blue',label=r'Full-coreNFW+tides')
        if (plot_evolve_low == 'yes'):
            f = open(evolve_file_low,'r')
            data = np.genfromtxt(f)
            f.close()
            plt.plot(data[:,0],data[:,1],':',\
                     linewidth=mylinewidth,\
                     color='blue',label=r'Lowrez-coreNFW+tides')
            mylegendfontsize = 15
        if (plot_nfw == 'yes'):
            f = open(evolve_file_nfw,'r')
            data = np.genfromtxt(f)
            f.close()
            plt.plot(data[:,0],data[:,1],'--',\
                     linewidth=mylinewidth,\
                     color='blue',label=r'NFW+tides')

    if (plotlegend == 'yes'):
        plt.legend(fontsize=mylegendfontsize,loc='lower left')
    plt.savefig(outdir+'output_rho.pdf',bbox_inches='tight')

    #Write the above data to a file for comparitive plotting later:
    f = open(outdir+'output_rho.txt','w')
    for i in range(len(rbin)):
        f.write('%f %f %f %f %f %f %f %f\n' % \
         (rbin[i],rho_int[0,i],rho_int[1,i],rho_int[2,i],rho_int[3,i],\
              rho_int[4,i],rho_int[5,i],rho_int[6,i]))
    f.close()

    #And the coreNFW parameters:
    f = open(outdir+'output_M200c200_chain.txt','w')
    for i in range(len(M200store)):
        f.write('%f %f %f %f %f %f %f\n' % \
                (M200store[i],cstore[i],nstore[i],rcstore[i],rtstore[i],\
                 delstore[i],vmaxstore[i]))
    f.close()

    if (plotminmax == 'yes'):
        f = open(outdir+'output_rhominmax.txt','w')
        for i in range(len(rbin)):
            f.write('%f %f %f\n' % \
                        (rbin[i],rhominmax[i,0],rhominmax[i,1]))
        f.close()

    #And the density exponent:
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    if (reduce_xticks == 'yes'):
        plt.locator_params(axis='x',nbins=4)

    plt.fill_between(np.log10(rbin),dlnrhodlnr_int[5,:],dlnrhodlnr_int[6,:],\
                         facecolor='black',alpha=alp3sig,\
                         edgecolor='none')
    plt.fill_between(np.log10(rbin),dlnrhodlnr_int[3,:],dlnrhodlnr_int[4,:],\
                         facecolor='black',alpha=0.33,\
                         edgecolor='none')
    plt.fill_between(np.log10(rbin),dlnrhodlnr_int[1,:],dlnrhodlnr_int[2,:],\
                         facecolor='black',alpha=0.66,\
                         edgecolor='none')
    plt.plot(np.log10(rbin),dlnrhodlnr_int[0,:],'k',linewidth=mylinewidth,\
                 label=r'Fit')

    if (plotminmax == 'yes'):
        plt.plot(np.log10(rbin),dlnrhodlnrminmax[:,0],'g',\
                 linewidth=4)
        plt.plot(np.log10(rbin),dlnrhodlnrminmax[:,1],'g',\
                 linewidth=2)

    #And overplot true model (if mock):
    if (overtrue == 'yes'):
        truedlnrhodlnr = alpbetgamdlnrhodlnr(ranal,rho0,r0,alp,bet,gam)
        plt.plot(np.log10(ranal),truedlnrhodlnr,'b--',linewidth=mylinewidth,\
                 label=r'True')
    if (singlesplit == 'single'):
        plt.plot([np.log10(Rhalf),np.log10(Rhalf)],[-100,100],'b',alpha=0.5,\
                     linewidth=mylinewidth)
    elif (singlesplit == 'split'):
        plt.plot([np.log10(Rhalf1),np.log10(Rhalf1)],[-100,100],\
                     color=colorpop1,alpha=0.5,\
                     linewidth=mylinewidth)
        plt.plot([np.log10(Rhalf2),np.log10(Rhalf2)],[-100,100],\
                     color=colorpop2,alpha=0.5,\
                     linewidth=mylinewidth)

    plt.xlabel(r'${\rm Log}_{10}[r/{\rm kpc}]$',\
                   fontsize=myfontsize)
    plt.ylabel(r'${\rm dln}\rho/{\rm dln}r$',\
                   fontsize=myfontsize)

    if (rbin_inner < 0):
        plt.xlim([np.log10(np.min(rbin)),np.log10(np.max(rbin))])
    else:
        plt.xlim([np.log10(rbin_inner),np.log10(rbin_outer)])

    plt.ylim([-4,0])

    if (plotlegend == 'yes'):
        plt.legend(fontsize=mylegendfontsize)
    plt.savefig(outdir+'output_dlnrhodlnr.pdf',bbox_inches='tight')

    #Write the above data to a file for comparitive plotting later:
    f = open(outdir+'output_dlnrhodlnr.txt','w')
    for i in range(len(rbin)):
        f.write('%f %f %f %f %f %f %f %f\n' % \
         (rbin[i],dlnrhodlnr_int[0,i],dlnrhodlnr_int[1,i],\
              dlnrhodlnr_int[2,i],dlnrhodlnr_int[3,i],\
              dlnrhodlnr_int[4,i],dlnrhodlnr_int[5,i],\
              dlnrhodlnr_int[6,i]))
    f.close()
    if (plotminmax == 'yes'):
        f = open(outdir+'output_dlnrhodlnrminmax.txt','w')
        for i in range(len(rbin)):
            f.write('%f %f %f\n' % \
                        (rbin[i],dlnrhodlnrminmax[i,0],dlnrhodlnrminmax[i,1]))
        f.close()

    #And the J-factor data:
    if (calc_Jfac == 'yes'):
        f = open(outdir+'output_Jfac.txt','w')
        for i in range(len(Jstore)):
            f.write('%f\n' % Jstore[i])
        f.close()

    if (rotation == 'yes'):
        #And make a plot of the Arot parameter if rotation is activated:
        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)
        nbin = 25
        n, bins, patches = plt.hist(Arotstore,bins=nbin,\
                       range=(Arotmin,Arotmax),\
                       normed=1, facecolor='b', \
                       histtype='bar',alpha=0.5)
        plt.xlabel(r'$A_{\rm rot}$',\
                       fontsize=myfontsize)
        plt.ylabel(r'$N$',fontsize=myfontsize)
        plt.xlim([Arotmin,Arotmax])
        plt.savefig(outdir+'output_Arot.pdf',bbox_inches='tight')

    if (virialshape == 'yes'):
        #And make a plot of the Virial shape parameters, if 
        #activated:
        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)
        nbin = 25

        if (data_file_phot_vs != ''):
            mph = np.sum(data_phot_vs[:,10])*vsfac
            print 'Effective number of photometric stars:', mph
        else:
            mph = 1.0

        if (singlesplit == 'single'):
            n, bins, patches = plt.hist(vs1store/mph,bins=nbin,\
                       range=(np.min(vs1store)/mph,\
                              np.max(vs1store)/mph),\
                       normed=1, facecolor='b', \
                       histtype='bar',alpha=0.5, \
                       label='vs_1')
            plt.errorbar([vs1bin/mph],[0.5*np.max(n)],\
                             xerr=vs1err/mph,fmt='ob')
        else:
            n, bins, patches = plt.hist(vs1_1store/mph,bins=nbin,\
                       range=(np.min(vs1_1store)/mph,\
                              np.max(vs1_1store)/mph),\
                       normed=1, facecolor=colorpop1, \
                       histtype='bar',alpha=0.5, \
                       label='vs_1_1')
            n, bins, patches = plt.hist(vs1_2store/mph,bins=nbin,\
                       range=(np.min(vs1_2store)/mph,\
                              np.max(vs1_2store)/mph),\
                       normed=1, facecolor=colorpop2, \
                       histtype='bar',alpha=0.5, \
                       label='vs_1_2')
            plt.errorbar([vs1bin1/mph],[0.5*np.max(n)],\
                             xerr=vs1err1/mph,fmt='o',color=colorpop1data)
            plt.errorbar([vs1bin2/mph],[0.5*np.max(n)],\
                             xerr=vs1err2/mph,fmt='o',color=colorpop2data)

        plt.xlabel(r'$v_{s1}/%d\,[{\rm km}^4\,{\rm s}^{-4}]$' % (vsfac),\
                   fontsize=myfontsize)
        plt.ylabel(r'$N$',fontsize=myfontsize)
        plt.ylim([0,np.max(n)])
        plt.savefig(outdir+'output_vs1.pdf',bbox_inches='tight')

        fig = plt.figure(figsize=(figx,figy))
        ax = fig.add_subplot(111)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(mylinewidth)
        plt.xticks(fontsize=myfontsize)
        plt.yticks(fontsize=myfontsize)
        nbin = 25

        if (singlesplit == 'single'):
            n, bins, patches = plt.hist(vs2store/mph,bins=nbin,\
                       range=(np.min(vs2store)/mph,\
                              np.max(vs2store)/mph),\
                       normed=1, facecolor='r', \
                       histtype='bar',alpha=0.5)
            plt.errorbar([vs2bin/mph],[0.5*np.max(n)],\
                     xerr=vs2err/mph,fmt='or')
        else:
            n, bins, patches = plt.hist(vs2_1store/mph,bins=nbin,\
                       range=(np.min(vs2_1store)/mph,\
                              np.max(vs2_1store)/mph),\
                       normed=1, facecolor=colorpop1, \
                       histtype='bar',alpha=0.5)
            n, bins, patches = plt.hist(vs2_2store/mph,bins=nbin,\
                       range=(np.min(vs2_2store)/mph,\
                              np.max(vs2_2store)/mph),\
                       normed=1, facecolor=colorpop2, \
                       histtype='bar',alpha=0.5)
            plt.errorbar([vs2bin1/mph],[0.5*np.max(n)],\
                     xerr=vs2err1/mph,fmt='o',color=colorpop1data)
            plt.errorbar([vs2bin2/mph],[0.5*np.max(n)],\
                     xerr=vs2err2/mph,fmt='o',color=colorpop2data)

        plt.xlabel(r'$v_{s2}/%d\,[{\rm km}^4\,{\rm s}^{-4}\,{\rm kpc}^2]$' % (vsfac),\
                   fontsize=myfontsize)
        plt.ylabel(r'$N$',fontsize=myfontsize)
        plt.ylim([0,np.max(n)])
        plt.savefig(outdir+'output_vs2.pdf',bbox_inches='tight')

    #Plot SIDM model historgrams here:
    nbin = 15
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    
    n, bins, patches = plt.hist(np.log10(M200store),bins=nbin,\
                                range=(logM200low,logM200high),\
                                normed=1, facecolor='b', \
                                histtype='bar',alpha=0.5)

    plt.xlabel(r'${\rm Log}_{10}[M_{200}/{\rm M}_\odot]$',\
               fontsize=myfontsize)
    plt.ylabel(r'$N$',fontsize=myfontsize)
    
    plt.xlim([logM200low,logM200high])
    plt.ylim([0,1])
    
    plt.savefig(outdir+'output_M200.pdf',bbox_inches='tight')

    nbin = 15
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)

    vmaxlow = 5.0
    vmaxhigh = 50.0
    
    n, bins, patches = plt.hist(np.log10(vmaxstore),bins=nbin,\
                                range=(vmaxlow,vmaxhigh),\
                                normed=1, facecolor='b', \
                                histtype='bar',alpha=0.5)

    plt.xlabel(r'$v_{\rm max}\,[{\rm km/s}]$',\
               fontsize=myfontsize)
    plt.ylabel(r'$N$',fontsize=myfontsize)

    plt.xlim([vmaxlow,vmaxhigh])
    plt.ylim([0,1])
    
    plt.savefig(outdir+'output_vmax.pdf',bbox_inches='tight')
    
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    
    n, bins, patches = plt.hist(cstore,bins=nbin,\
                                range=(clow,chigh),\
                                normed=1, facecolor='b', \
                                histtype='bar',alpha=0.5)

    plt.xlabel(r'$c$',fontsize=myfontsize)
    plt.ylabel(r'$N$',fontsize=myfontsize)

    plt.xlim([clow,chigh])
    plt.ylim([0,0.25])

    plt.savefig(outdir+'output_c.pdf',bbox_inches='tight')
  
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize2)
    plt.yticks(fontsize=myfontsize2)
    if (use_rclinear == 'no'):
        ax.xaxis.set_ticks(np.arange(logrclow,logrchigh+1.0,1.0))
        n, bins, patches = plt.hist(np.log10(rcstore),bins=nbin,\
                                    range=(logrclow,logrchigh),\
                                    normed=1, facecolor='k', \
                                    histtype='bar',alpha=0.5)
    else:
        n, bins, patches = plt.hist(rcstore,bins=nbin,\
                                    range=(rclow,rchigh),\
                                    normed=1, facecolor='k', \
                                    histtype='bar',alpha=0.5)
        
    #Overlay true answer:
    if (rcmock > 0.0):
        if (use_rclinear == 'no'):
            plt.axvline(x=np.log10(rcmock),color='blue',\
                        linewidth=mylinewidth)
        else:
            plt.axvline(x=rcmock,color='blue',\
                        linewidth=mylinewidth)

    if (use_rclinear == 'no'):
        plt.xlabel(r'${\rm Log}_{10}[r_c/{\rm kpc}]$',fontsize=myfontsize2)
        plt.xlim([logrclow,logrchigh])
    else:
        plt.xlabel(r'$r_c({\rm kpc})$',fontsize=myfontsize2)
        plt.xlim([0,1.0])    
    plt.ylabel(r'$N$',fontsize=myfontsize2) 
    plt.ylim([0,2])
    plt.savefig(outdir+'output_rc.pdf',bbox_inches='tight')

    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)

    if (use_rtlinear == 'no'):
        n, bins, patches = plt.hist(np.log10(rtstore),bins=nbin,\
                                        range=(logrtlow,\
                                               logrthigh),\
                                    normed=1, facecolor='k', \
                                    histtype='bar',alpha=0.5)
    else:
        n, bins, patches = plt.hist(rtstore,bins=nbin,\
                                        range=(rtlow,\
                                               rthigh),\
                                    normed=1, facecolor='k', \
                                    histtype='bar',alpha=0.5)

    if (use_rtlinear == 'no'):
        plt.xlabel(r'${\rm Log}_{10}[r_t/{\rm kpc}]$',fontsize=myfontsize)
    else:
        plt.xlabel(r'$r_t({\rm kpc})$',fontsize=myfontsize)
    plt.ylabel(r'$N$',fontsize=myfontsize)
    plt.savefig(outdir+'output_rt.pdf',bbox_inches='tight')
        
    sigmstore = np.zeros(len(rcstore))
    for i in range(len(rcstore)):
        sigmstore[i] = sidm_novel3(rcstore[i],M200store[i],cstore[i],\
                                   oden,rhocrit,rtstore[i])
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize2)
    plt.yticks(fontsize=myfontsize2)

    n, bins, patches = plt.hist(sigmstore,bins=nbin,\
                  range=(0.0,1.0),\
                  normed=True,facecolor='k', \
                  histtype='bar',alpha=0.5)

    #Overlay true answer:
    if (sigmmock > 0.0):
        plt.axvline(x=sigmmock,color='blue',\
                        linewidth=mylinewidth)

    plt.xlim([0.0,1.0])
    plt.xlabel(r'$\sigma/m\,({\rm cm}^2/{\rm g})$',fontsize=myfontsize2)
    plt.ylabel(r'$N$',fontsize=myfontsize2)
    plt.savefig(outdir+'output_sigm.pdf',bbox_inches='tight')

    #Calculate 95% and 99% intervals:
    index = np.argsort(rcstore,axis=0)
    print 'rc 95% upper:', \
        rcstore[index[np.int(97.5/100. * len(rcstore))]]
    print 'rc 99% upper:', \
        rcstore[index[np.int(99.5/100. * len(rcstore))]]
    index = np.argsort(sigmstore,axis=0)
    print 'sigm 95% upper:', \
        sigmstore[index[np.int(97.5/100. * len(sigmstore))]]
    print 'sigm 99% upper:', \
        sigmstore[index[np.int(99.5/100. * len(sigmstore))]]

    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    
    n, bins, patches = plt.hist(nstore,bins=nbin,\
                                range=(nlow,nhigh),\
                                normed=1, facecolor='b', \
                                histtype='bar',alpha=0.5)

    plt.xlabel(r'$n$',\
               fontsize=myfontsize)
    plt.ylabel(r'$N$',fontsize=myfontsize)
    
    plt.xlim([nlow,nhigh])
    plt.ylim([0,0.25])
    plt.savefig(outdir+'output_n.pdf',bbox_inches='tight')

    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=20, width=2, which='major')
    ax.tick_params('both', length=10, width=1, which='minor')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)

    n, bins, patches = plt.hist(delstore,bins=nbin,\
                                range=(dellow,delhigh),\
                                normed=1, facecolor='b', \
                                histtype='bar',alpha=0.5)

    plt.xlabel(r'$\delta$',\
               fontsize=myfontsize)
    plt.ylabel(r'$N$',fontsize=myfontsize)

    plt.xlim([dellow,delhigh])
    plt.ylim([0,0.25])
    plt.savefig(outdir+'output_del.pdf',bbox_inches='tight')

    #Calculate M200 +/- 68%:
    M200med, M200sixlow, M200sixhi,\
        M200ninelow, M200ninehi, \
        M200nineninehi, M200nineninelow = calcmedquartnine(M200store)
    print '*******************************'
    print 'M200 -/+ 68% :: ', M200med, M200sixlow, M200sixhi

    #And same for vmax:
    vmaxmed, vmaxsixlow, vmaxsixhi,\
        vmaxninelow, vmaxninehi, \
        vmaxnineninehi, vmaxnineninelow = calcmedquartnine(vmaxstore)
    print '*******************************'
    print 'vmax -/+ 68% :: ', vmaxmed, vmaxsixlow, vmaxsixhi            
    
    #And the same for rt: 
    rtmed, rtsixlow, rtsixhi,\
        rtninelow, rtninehi, \
        rtnineninehi, rtnineninelow = calcmedquartnine(rtstore)
    print '*******************************'
    print 'rt -/+ 68% :: ', rtmed, rtsixlow, rtsixhi
        
