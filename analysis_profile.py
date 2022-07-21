import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/elizabeth/lens_codes_v3.7')
from models_profiles import *
from fit_profiles_curvefit import *
from fit_profiles_curvefit import Delta_Sigma_fit
from models_profiles import *
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 12})

cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

folder = '../profiles3/'   
# folder = '../profiles/'   
meanmag = np.array([-20.03361442, -20.57903032, -21.2230643 , -21.84807599,-22.48726666])
lM200  = np.array([11.32221929, 11.54406804, 11.92427929, 12.22530928, 12.67117284])

def color_plot():
    
    from medianas import separate_medianas
    
    pcat = '_photo_z_2nd_run_mag_i'
    
    L1 = np.loadtxt('../catlogoscon5log10h/Pares-PAUS_W1-Photo_z_calibrate'+pcat).T
    field = np.ones(len(L1[1]))*1                             
    L1 = np.vstack((L1,field))                                
                                                            
    L2 = np.loadtxt('../catlogoscon5log10h/Pares-PAUS_W2-Photo_z_calibrate'+pcat).T
    field = np.ones(len(L2[1]))*2                             
    L2 = np.vstack((L2,field))                                
                                    
    L3 = np.loadtxt('../catlogoscon5log10h/Pares-PAUS_W3-Photo_z_calibrate'+pcat).T
    field = np.ones(len(L3[1]))*3
    L3 = np.vstack((L3,field))
    
    L = np.vstack((L1.T,L2.T,L3.T)).T
    
    M1 = L[8]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    M2 = L[-2]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    Mtot = -2.5*np.log10(10**(-0.4*M1)+10**(-0.4*M2))
    M1i = L[7]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    M2i = L[-3]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    Mtoti = -2.5*np.log10(10**(-0.4*M1i)+10**(-0.4*M2i))
    
    label_x = '$M_{r_{par}}$'
    label_y = '$M_{i_{par}}$'
    
    separate_medianas(Mtot,Mtot-Mtoti,label_x = label_x, label_y = label_y, out_plot = '../final_plots/color_mag.pdf')


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
                       RIN=300,ROUT=10000, fytpe = '',
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
    Mmean = h['M_mean']    
    
    ndots = p.shape[0]
    
    
    # FIT MONOPOLE
    fitted = fits.open(folder+'fitresults'+ftype+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)
    fitpar = fitted[0].header
    lgM   = fitted[1].data.logM
    
    rplot = np.logspace(np.log10(h['RIN']/1000.),np.log10(h['ROUT']/1000.),20)
    
    nfw  = Delta_Sigma_fit(p.Rp,p.DSigma_T,p.error_DSigma_T,zmean,cosmo_as)
    
    mass = '$'+str(np.round(fitpar['lM200'],2))+'$'
    cfit = '$'+str(np.round(fitpar['c200'],2))+'$'
    
    
    if plot:
        
        ds2h   = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='2h')    
        ds1h   = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='1h')    
        
        ds = ds1h+ds2h
        
        axDS.plot(p.Rp,p.DSigma_T,'C1o')
        axDS.errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,ecolor='C1',fmt='None')
        axDS.plot(rplot,ds,'C3',label=lsamp+' $\log M_{200}= $'+mass+' $c_{200} = $'+cfit)
        axDS.plot(rplot,ds1h,'C4')
        axDS.plot(rplot,ds2h,'C4--')
        # axDS.fill_between(p.Rp,p.DSigma_T+error_DST,p.DSigma_T-error_DST,color='C1',alpha=0.4)
        axDS.set_xscale('log')
        axDS.set_yscale('log')
        if ylabel:
            axDS.set_ylabel(r'$\Delta\Sigma_{T} [M_{\odot}pc^{-2} h ]$')
        axDS.set_xlabel(r'$R [Mpc/h]$')
        axDS.set_ylim(0.095,200)
        axDS.set_xlim(h['RIN']/1000.,h['ROUT']/1000.)
        axDS.yaxis.set_ticks([0.1,1,10,100])
        axDS.set_yticklabels([0.1,1,10,100])
        axDS.axvline(RIN/1000.,color='C7',alpha=0.5,ls='--')
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
    
    return np.append(np.percentile(lgM[2500:], [16,50,84]),Mmean)


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


# pcats = ['_zspec',
        # '_zspec_best',
        # '_photo_z_2nd_run_mag_i',
        # '_photo_z_2nd_run_mag_i_best']

ftype = '_boost'

pcat = '_photo_z_2nd_run_mag_i'
best = '_photo_z_2nd_run_mag_i_best'


lMfit = []

fDS, axDS = plt.subplots(9,2, figsize=(10,18),sharex = True,sharey = True)
fDS.subplots_adjust(hspace=0,wspace=0)

fC, axC = plt.subplots(9,2, figsize=(12,14),sharex = True,sharey = True)
fC.subplots_adjust(hspace=0,wspace=0)

axDS = axDS.flatten()
axC = axC.flatten()



# samp =  ['mh_all_'+pcat,'mh_w3_'+pcat,
        # 'mh_LrM_all_'+pcat,'mh_LrM_w3_'+pcat,
        # 'mh_Lrm_all_'+pcat,'mh_Lrm_w3_'+pcat,
        # 'mh_zM_all_'+pcat,'mh_zM_w3_'+pcat,
        # 'mh_zm_all_'+pcat,'mh_zm_w3_'+pcat]
        
# samp =  ['mh_all_'+pcat,'mh_w3_'+pcat,
        # 'mh_Mm_all_'+pcat,'mh_Mm_w3_'+pcat,
        # 'mh_MM_all_'+pcat,'mh_MM_w3_'+pcat,
        # 'mh_red_all_'+pcat,'mh_red_w3_'+pcat,
        # 'mh_blue_all_'+pcat,'mh_blue_w3_'+pcat]

# samp =  ['mh_all_'+pcat,'mh_w3_'+pcat,
        # 'mh_LrM_all_'+pcat,'mh_LrM_w3_'+pcat,
        # 'mh_Lrm_all_'+pcat,'mh_Lrm_w3_'+pcat,
        # 'mh_zM_all_'+pcat,'mh_zM_w3_'+pcat,
        # 'mh_zm_all_'+pcat,'mh_zm_w3_'+pcat,
        # 'mh_MM_all_'+pcat,'mh_MM_w3_'+pcat,
        # 'mh_Mm_all_'+pcat,'mh_Mm_w3_'+pcat,
        # 'mh_red_all_'+pcat,'mh_red_w3_'+pcat,
        # 'mh_blue_all_'+pcat,'mh_blue_w3_'+pcat]
        
samp =  ['mh_all_'+pcat,'mh_all_'+best,
        'mh_LrM_all_'+pcat,'mh_LrM_all_'+best,
        'mh_Lrm_all_'+pcat,'mh_Lrm_all_'+best,
        'mh_zM_all_'+pcat,'mh_zM_all_'+best,
        'mh_zm_all_'+pcat,'mh_zm_all_'+best,
        'mh_MM_all_'+pcat,'mh_MM_all_'+best,
        'mh_Mm_all_'+pcat,'mh_Mm_all_'+best,
        'mh_red_all_'+pcat,'mh_red_all_'+best,
        'mh_blue_all_'+pcat,'mh_blue_all_'+best]

lsamp = ['$all$ - ','$all$ - ',
         '$L_2/L_1 > 0.8$ - ','$L_2/L_1 > 0.8$ - ',
         '$L_2/L_1 < 0.8$ - ','$L_2/L_1 < 0.8$ - ',
         '$z > 0.4$ - ','$z > 0.4$ - ',
         '$z < 0.4$ - ','$z < 0.4$ - ',
         '$M_{r_{par}} > -21.0$ -','$M_{r_{par}} > -21.0$ -',
         '$M_{r_{par}} < -21.0$ - ','$M_{r_{par}} < -21.0$ - ',
         r'$red\,\,pairs$ - ',r'$red\,\,pairs$ - ',
         r'$blue\,\,pairs$ - ',r'$blue\,\,pairs$ - ']
        
# lsamp = ['all -','',
        # 'M < -21.0 -','',
        # 'M > -21.0 -','',
        # 'red pairs -','',
        # 'blue pairs -','']

ylabel = [True,False]*len(samp)

for j in range(len(samp)):
    lMfit += [plt_profile_fit_2h(samp[j],
              lsamp[j],axDS[j],axC[j],
              fytpe = ftype, ylabel=ylabel[j])]
    # lMfit += [plt_profile_fit_2h(samp[j],lsamp[j],plot=False,fytpe = ftype)]



fDS.savefig('../final_plots/profile.pdf',bbox_inches='tight')
fC.savefig('../final_plots/chains2'+pcat+ftype+'.png',bbox_inches='tight')

'''

csamp = ['k',
         'palevioletred','palevioletred',
         'royalblue','royalblue',
         'gold','gold',
         'C9','C9']

         
mark = ['o']+['^','v']*4


plt.figure()
plt.title(pcat)
for j in np.arange(len(csamp)):
    plt.errorbar(lMfit[2*j][1],lMfit[2*j+1][1],
                     xerr=np.array([np.diff(lMfit[2*j])[:-1]]).T,
                     yerr=np.array([np.diff(lMfit[2*j+1])[:-1]]).T,
                     fmt=csamp[j],label=lsamp[j],marker=mark[j])
plt.legend(frameon=False,loc=4)
plt.plot([11.5,13.3],[11.5,13.3],'C7--')
plt.xlabel('$\log M_{200}$')
plt.ylabel('$\log M_{200}(W3)$')
plt.savefig('../compare2_w3'+pcat+ftype+'.png',bbox_inches='tight')

plt.figure()
plt.title('All - '+pcat)
for j in np.arange(len(csamp)):
    # plt.plot(meanmag,lM200,'k')
    plt.plot(meanmag-5.*np.log10(0.69),lM200,'k--')
    plt.errorbar(lMfit[2*j][-1],lMfit[2*j][1],
                     yerr=np.array([np.diff(lMfit[2*j])[:-1]]).T,
                     fmt=csamp[j],label=lsamp[j],marker=mark[j])
plt.legend(frameon=False,loc=3)
plt.xlabel(r'$\langle M_r \rangle$')
plt.ylabel('$\log M_{200}$')
plt.axis([-22.1,-19.3,10.2,13.2])
plt.savefig('../Mag_lM200_'+pcat+ftype+'.png',bbox_inches='tight')

plt.figure()
plt.title('W3 - '+pcat)
for j in np.arange(len(csamp)):
    # plt.plot(meanmag,lM200,'k')
    plt.plot(meanmag-5.*np.log10(0.69),lM200,'k--')
    plt.errorbar(lMfit[2*j][-1],lMfit[2*j+1][1],
                     yerr=np.array([np.diff(lMfit[2*j+1])[:-1]]).T,
                     fmt=csamp[j],label=lsamp[j],marker=mark[j])
plt.legend(frameon=False,loc=3)
plt.xlabel(r'$\langle M_r \rangle$')
plt.ylabel('$\log M_{200}$')
plt.axis([-22.1,-19.3,10.2,13.2])
plt.savefig('../Mag_lM200_w3_'+pcat+ftype+'.png',bbox_inches='tight')
    
'''
