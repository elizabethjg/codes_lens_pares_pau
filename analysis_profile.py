import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/elizabeth/lens_codes_v3.7')
from models_profiles import *
from fit_profiles_curvefit import *
from fit_profiles_curvefit import Delta_Sigma_fit
from models_profiles import *

folder = '../profiles_new/'   
# folder = '../profiles/'   

def plt_profile_wofit(samp):
    

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    cosmo_as = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)
    '''
    
    h = profile[1].header
    p = profile[1].data
    '''
    print(h['N_LENSES'])
    
    ### compute dilution
    bines = np.logspace(np.log10(h['RIN']),np.log10(h['ROUT']),num=len(p)+1)
    area = np.pi*np.diff(bines**2)
    
    ngal = p.NGAL_w
    
    d = ngal/area
    
    fcl = ((d - np.mean(d[-2:]))*area)/ngal
    
    bcorr = 1./(1-fcl)

    p.DSigma_T = bcorr*p.DSigma_T
    p.DSigma_X = bcorr*p.DSigma_X
    p.error_DSigma_T = bcorr*p.error_DSigma_T
    p.error_DSigma_X = bcorr*p.error_DSigma_X
    
    zmean = h['z_mean']    
    
    CovDST  = cov.COV_ST.reshape(len(p),len(p))
    CovDSX  = cov.COV_SX.reshape(len(p),len(p))
    
    mr = (p.Rp > 0.3)*(p.Rp < 2.)
    
    # FIT MONOPOLE
    nfw     = Delta_Sigma_fit(p.Rp[mr],p.DSigma_T[mr],np.sqrt(np.diag(CovDST))[mr],zmean,cosmo_as)
    
    
    mass = str(np.round(np.log10(nfw.M200),2))
    cfit = str(np.round((nfw.c200),2))
    
    f, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(6,6),sharex = True)
    f.subplots_adjust(hspace=0)
                    
    ax[0].plot(p.Rp,p.DSigma_T,'C1o')
    ax[0].errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,ecolor='C1',fmt='None')
    ax[0].plot(nfw.xplot,nfw.yplot,'C3',label='fitted nfw $\log M_{200}=$'+mass+' $c_{200} = $'+cfit)
    ax[0].fill_between(p.Rp,p.DSigma_T+p.error_DSigma_T,p.DSigma_T-p.error_DSigma_T,color='C1',alpha=0.3)
    ax[0].fill_between(p.Rp,p.DSigma_T+np.sqrt(np.diag(CovDST))*bcorr,p.DSigma_T-np.sqrt(np.diag(CovDST))*bcorr,color='C2',alpha=0.3)
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'$\Delta\Sigma_{T} [M_{\odot}pc^{-2} h ]$')
    ax[0].set_ylim(0.05,100)
    ax[0].set_xlim(h['RIN']/1000.,h['ROUT']/1000.)
    ax[0].yaxis.set_ticks([0.1,1,10,100])
    ax[0].set_yticklabels([0.1,1,10,100])
    ax[0].legend()

    ax[1].plot(p.Rp,p.DSigma_X,'C7x')
    ax[1].errorbar(p.Rp,p.DSigma_X,yerr=p.error_DSigma_X,ecolor='C7',fmt='None')
    ax[1].fill_between(p.Rp,p.DSigma_X+p.error_DSigma_X,p.DSigma_X-p.error_DSigma_X,color='C7',alpha=0.3)
    ax[1].fill_between(p.Rp,p.DSigma_X+np.sqrt(np.diag(CovDSX))*bcorr,p.DSigma_X-np.sqrt(np.diag(CovDSX))*bcorr,color='C2',alpha=0.3)
    ax[1].set_ylabel(r'$\Delta\Sigma_{\times} [M_{\odot}pc^{-2} h ]$')
    ax[1].set_xlabel(r'$R [Mpc/h]$')
    ax[1].plot([0,h['ROUT']/1000.],[0,0],'k--')
    ax[1].set_ylim(-20,20)
    
    f.savefig(folder+'plots/profile_'+samp+'.png',bbox_inches='tight')
    

def plt_profile_align(samp):
    

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[1].header
    p   = profile[1].data
    
    cosmo_as = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)
    '''
    
    h = profile[1].header
    p = profile[1].data
    '''

    ### compute dilution
    bines = np.logspace(np.log10(h['RIN']),np.log10(h['ROUT']),num=len(p)+1)
    area = np.pi*np.diff(bines**2)
    
    ngal = p.NGAL_w
    
    d = ngal/area
    
    fcl = ((d - np.mean(d[-2:]))*area)/ngal
    
    bcorr = 1./(1-fcl)

    p.DSigma_T = bcorr*p.DSigma_T
    p.DSigma_X = bcorr*p.DSigma_X
    p.error_DSigma_T = bcorr*p.error_DSigma_T
    p.error_DSigma_X = bcorr*p.error_DSigma_X
    
    zmean = h['z_mean']    
    
    ndots = p.shape[0]
    
    
    # FIT MONOPOLE
    nfw     = Delta_Sigma_fit(p.Rp,p.DSigma_T,p.error_DSigma_T,zmean,cosmo_as)
    
    
    mass = str(np.round(np.log10(nfw.M200),2))
    cfit = str(np.round((nfw.c200),2))
    
    f, ax = plt.subplots(2, 1, figsize=(6,6),sharex = True,sharey=True)
    f.subplots_adjust(hspace=0)
                    
    ax[0].plot(p.Rp,p.DSigma_T,'C1o')
    ax[0].errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,ecolor='C1',fmt='None')
    ax[0].plot([0,h['ROUT']/1000.],[0,0],'C3')
    ax[0].fill_between(p.Rp,p.DSigma_T+p.error_DSigma_T,p.DSigma_T-p.error_DSigma_T,color='C1',alpha=0.3)
    ax[0].set_xscale('log')
    # ax[0].set_yscale('log')
    ax[0].set_ylabel(r'$\Delta\Sigma_{T} [M_{\odot}pc^{-2} h ]$')
    ax[0].set_xlim(h['RIN']/1000.,h['ROUT']/1000.)
    # ax[0].yaxis.set_ticks([0.1,1,10,100])
    # ax[0].set_yticklabels([0.1,1,10,100])
    ax[0].legend()

    ax[1].plot(p.Rp,p.DSigma_X,'C7x')
    ax[1].errorbar(p.Rp,p.DSigma_X,yerr=p.error_DSigma_X,ecolor='C7',fmt='None')
    ax[1].fill_between(p.Rp,p.DSigma_X+p.error_DSigma_X,p.DSigma_X-p.error_DSigma_X,color='C7',alpha=0.3)
    ax[1].set_ylabel(r'$\Delta\Sigma_{\times} [M_{\odot}pc^{-2} h ]$')
    ax[1].set_xlabel(r'$R [Mpc/h]$')
    ax[1].plot([0,h['ROUT']/1000.],[0,0],'k--')
    ax[1].set_ylim(-100,100)
    
    f.savefig(folder+'plots/profile_'+samp+'.png',bbox_inches='tight')


def plt_profile_fit_2h(samp,lsamp,
                       axDS = plt,axC = plt,
                       RIN=300,ROUT=10000,
                       ylabel = True, plot = True):
    

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data

    CovDST  = cov.COV_ST.reshape(len(p),len(p))
    CovDSX  = cov.COV_SX.reshape(len(p),len(p))
    
    cosmo_as = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

    error_DSX = np.sqrt(np.diag(CovDSX))
    error_DST = np.sqrt(np.diag(CovDST))

    ### compute dilution
    bines = np.logspace(np.log10(h['RIN']),np.log10(h['ROUT']),num=len(p)+1)
    area = np.pi*np.diff(bines**2)
    
    ngal = p.NGAL_w
    
    d = ngal/area
    
    fcl = ((d - np.mean(d[-2:]))*area)/ngal
    
    bcorr = 1./(1-fcl)

    p.DSigma_T = bcorr*p.DSigma_T
    p.DSigma_X = bcorr*p.DSigma_X
    p.error_DSigma_T = bcorr*p.error_DSigma_T
    p.error_DSigma_X = bcorr*p.error_DSigma_X

    
    zmean = h['z_mean']    
    
    ndots = p.shape[0]
    
    
    # FIT MONOPOLE
    fitted = fits.open(folder+'fitresults_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)
    fitpar = fitted[0].header
    lgM   = fitted[1].data.logM
    
    rplot = np.logspace(np.log10(h['RIN']/1000.),np.log10(h['ROUT']/1000.),20)
    
    nfw  = Delta_Sigma_fit(p.Rp,p.DSigma_T,p.error_DSigma_T,zmean,cosmo_as)
    # ds2h   = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='2h')    
    ds1h   = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='1h')    
    
    ds = ds1h #+ds2h
    
    mass = str(np.round(fitpar['lM200'],2))
    cfit = str(np.round(fitpar['c200'],2))
                        
    axDS.plot(p.Rp,p.DSigma_T,'C1o')
    axDS.errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,ecolor='C1',fmt='None')
    axDS.plot(rplot,ds,'C3',label=lsamp+' $\log M_{200}=$'+mass+' $c_{200} = $'+cfit)
    axDS.plot(rplot,ds1h,'C4')
    # axDS.plot(rplot,ds2h,'C4--')
    axDS.fill_between(p.Rp,p.DSigma_T+error_DST,p.DSigma_T-error_DST,color='C1',alpha=0.4)
    axDS.set_xscale('log')
    axDS.set_yscale('log')
    if ylabel:
        axDS.set_ylabel(r'$\Delta\Sigma_{T} [M_{\odot}pc^{-2} h ]$')
    axDS.set_xlabel(r'$R [Mpc/h]$')
    axDS.set_ylim(0.005,100)
    axDS.set_xlim(h['RIN']/1000.,h['ROUT']/1000.)
    axDS.yaxis.set_ticks([0.01,0.1,1,10,100])
    axDS.set_yticklabels([0.01,0.1,1,10,100])
    axDS.axvline(RIN/1000.,color='C7')
    axDS.axvline(ROUT/1000.,color='C7')    
    axDS.legend(frameon=False,loc=1)
    
    
    axC.plot(lgM,label=lsamp,alpha=0.5)
    axC.axhline(np.median(lgM[2500:]),alpha=0.5)
    axC.axhline(np.percentile(lgM[2500:], [16,50,84])[0],ls='--')
    axC.axhline(np.percentile(lgM[2500:], [16,50,84])[2],ls='--')
    axC.axvline(2500)
    axC.legend(frameon=False,loc=1)
    axC.set_xlabel('$N_{it}$')
    if ylabel:
        axC.set_ylabel('$\log M_{200}$')
    
    return np.percentile(lgM[2500:], [16,50,84])


def dilution(samp):
    

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[1].header
    p   = profile[1].data

    bines = np.logspace(np.log10(h['RIN']/1000.),np.log10(h['ROUT']/1000.),num=len(p)+1)
    area = np.pi*np.diff(bines**2)
    
    ngal = p.NGAL_w

    d = ngal/area

    fcl = ((d - np.mean(d[-2:]))*area)/ngal
    
    # plt.plot(p.Rp,1./(1-fcl))
    # plt.plot(p.Rp,1./(1-fcl),'o')
    
    return p.Rp,d,fcl,h['N_LENSES']
    
def test():
    
    samp1 = 'new_w1__photo_z_2nd_run_mag_i'
    samp2 = 'new_w2__photo_z_2nd_run_mag_i'
    samp3 = 'new_w3__photo_z_2nd_run_mag_i'
    
    r,d1,fcl1,n1 = dilution(samp1)
    r,d2,fcl2,n2 = dilution(samp2)
    r,d3,fcl3,n3 = dilution(samp3)
    
    plt.plot(r,d1/n1,label='W1')
    plt.plot(r,d2/n2,label='W2')
    plt.plot(r,d3/n3,label='W3')
    plt.legend()
    plt.xscale('log')
    plt.xlabel('R [Mpc]')
    plt.ylabel('n')


fDS, axDS = plt.subplots(5,2, figsize=(12,14),sharex = True,sharey = True)
fDS.subplots_adjust(hspace=0,wspace=0)

fC, axC = plt.subplots(5,2, figsize=(12,14),sharex = True,sharey = True)
fC.subplots_adjust(hspace=0,wspace=0)

axDS = axDS.flatten()
axC = axC.flatten()

pcat = '_zspec'
pcat = '_photo_z_2nd_run_mag_i'
pcat = '_photo_z_2nd_run_mag_i_best'

samp =  ['wc_all_'+pcat,'wc_w3_'+pcat,
         'wc_LrM_all_'+pcat,'wc_LrM_w3_'+pcat,
         'wc_Lrm_all_'+pcat,'wc_Lrm_w3_'+pcat,
         'wc_zM_all_'+pcat,'wc_zM_w3_'+pcat,
         'wc_zm_all_'+pcat,'wc_zm_w3_'+pcat]

lsamp = ['all -','',
         'HLratio -','',
         'LLratio -','',
         'Hz -','',
         'Lz -','']

ylabel = [True,False]*5

for j in range(len(axDS)):
    plt_profile_fit_2h(samp[j],lsamp[j],axDS[j],axC[j],ylabel=ylabel[j])


fDS.savefig('../profile'+pcat+'.png',bbox_inches='tight')
fC.savefig('../chains'+pcat+'.png',bbox_inches='tight')



'''
plt.figure()    
plt.plot(w1_sources.RAJ2000,w1_sources.DECJ2000,',')
plt.plot(L1[4],L1[5],'.')
plt.xlabel('R.A.')
plt.xlabel('Dec.')
plt.savefig('/home/elizabeth/PARES-PAU/W1.png')

plt.figure()
plt.plot(w2_sources.RAJ2000,w2_sources.DECJ2000,',')
plt.plot(L2[4],L2[5],'.')
plt.xlabel('R.A.')
plt.xlabel('Dec.')
plt.savefig('/home/elizabeth/PARES-PAU/W2.png')

plt.figure()
plt.plot(w3_sources.RAJ2000,w3_sources.DECJ2000,',')
plt.plot(L3[4],L3[5],'.')
plt.xlabel('R.A.')
plt.xlabel('Dec.')
plt.savefig('/home/elizabeth/PARES-PAU/W3.png')
'''
