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
from matplotlib import cbook
import pandas as pd
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
# matplotlib.rcParams.update({'font.size': 14})

cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

folder = '../profiles_pk/'   
# folder = '../profiles/'   
meanmag = np.array([-20.03361442, -20.57903032, -21.2230643 , -21.84807599,-22.48726666])
lM200  = np.array([11.32221929, 11.54406804, 11.92427929, 12.22530928, 12.67117284])

meanmag_p1 = np.array([-21.44801166, -22.24667956])
lM200_p1   = np.array([12.1487476 , 12.76705489])
elM200_p1  = np.array([0.22447506, 0.29605648])


pcat = 'Total'
best = 'Gold'
vane = '_vane'
ftype = '2_boost'

def compute_density(ra,dec):
    
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)
    
    area = 0.7*(np.sin(dec.max())-np.sin(dec.min()))*(ra.max()-ra.min())*(180/np.pi)**2
    
    return len(ra)/(area*3600.)

def color_plot_vane():
    
    from medianas import separate_medianas
    
    # LOAD PAUS CATALOGS
    
    w1_deep = pd.read_csv('../pareconk/w1_DEEPZ_IATE.csv')
    w1 = pd.read_csv('../pareconk/Photo_z_calibrate_iterative_W1_no_k_corr.csv')
    w3 = pd.read_csv('../pareconk/Photo_z_calibrate_iterative_W3_no_k_corr.csv')
    w1_best = pd.read_csv('../pareconk/Photo_z_calibrate_iterative_W1_best_pz_no_k_corr.csv')
    w3_best = pd.read_csv('../pareconk/Photo_z_calibrate_iterative_W3_best_pz_no_k_corr.csv')
    
    mr  = np.append(w1.mag_r_no_k,w3.mag_r_no_k)
    mrv = w1_deep.mag_r_K
    mrb = np.append(w1_best.mag_r_no_k,w3_best.mag_r_no_k)
    zp  = np.append(w1.zb,w3.zb)
    zpb  = np.append(w1_best.zb,w3_best.zb)
    zpv  = w1_deep.zb

    mp  = (mr > -23.5) & (zp < 0.6) & (zp > 0.2)
    mpb = (mrb > -23.5) & (zpb < 0.6) & (zpb > 0.2)
    mpv = (mrv > -23.5) & (zpv < 0.6) & (zpv > 0.2)
    
    mr  = mr[mp]
    mrv  = mrv[mpv]
    mrb  = mrb[mpb]
    
    zp   = zp[mp]
    zpv  = zpv[mpv]
    zpb  = zpb[mpb]
    
    
    mi  = np.append(w1.mag_i_no_k,w3.mag_i_no_k)[mp]
    mib = np.append(w1_best.mag_i_no_k,w3_best.mag_i_no_k)[mpb]
    miv = w1_deep.mag_i_K[mpv]
    
    Mr  = mr - 5.*np.log10(np.array(cosmo.luminosity_distance(zp))*1.e6) + 5
    Mi  = mi - 5.*np.log10(np.array(cosmo.luminosity_distance(zp))*1.e6) + 5
    
    Mrb = mrb - 5.*np.log10(np.array(cosmo.luminosity_distance(zpb))*1.e6) + 5
    Mib = mib - 5.*np.log10(np.array(cosmo.luminosity_distance(zpb))*1.e6) + 5

    Mrv = mrv - 5.*np.log10(np.array(cosmo.luminosity_distance(zpv))*1.e6) + 5
    Miv = miv - 5.*np.log10(np.array(cosmo.luminosity_distance(zpv))*1.e6) + 5

    colorp  = Mr - Mi
    colorpb = Mrb - Mib
    colorpv = Mrv - Miv
        
    # LOAD PAIRS CATALOGUES
    
    L1 = np.loadtxt('../pareconk/Pares-PAUS_W1-'+pcat).T                                                            
    # L2 = np.loadtxt('../pareconk/Pares-PAUS_W2-'+pcat).T                                    
    L3 = np.loadtxt('../pareconk/Pares-PAUS_W3-'+pcat).T
    
    print('density W1 = ',compute_density(L1[1],L1[2]))
    # print('density W2 = ',compute_density(L2[1],L2[2]))
    print('density W3 = ',compute_density(L3[1],L3[2]))
    
    # L = np.vstack((L1.T,L2.T,L3.T)).T
    L = np.vstack((L1.T,L3.T)).T
    
    M1 = L[8]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    M2 = L[13]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    Mtot = -2.5*np.log10(10**(-0.4*M1)+10**(-0.4*M2))
    M1i = L[7]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    M2i = L[12]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    Mtoti = -2.5*np.log10(10**(-0.4*M1i)+10**(-0.4*M2i))
    
    rp = L[-2]
    Vr = L[-1]*cvel*1.e-3
    
    color = Mtot-Mtoti

    Lv = np.loadtxt('../pareconk/Pares-PAUS_W1-Vane').T
    M1v = Lv[8]-5.*np.log10(np.array(cosmo.luminosity_distance(Lv[3]))*1.e6)+5
    M2v = Lv[13]-5.*np.log10(np.array(cosmo.luminosity_distance(Lv[3]))*1.e6)+5
    Mtotv = -2.5*np.log10(10**(-0.4*M1v)+10**(-0.4*M2v))
    M1iv = Lv[7]-5.*np.log10(np.array(cosmo.luminosity_distance(Lv[3]))*1.e6)+5
    M2iv = Lv[12]-5.*np.log10(np.array(cosmo.luminosity_distance(Lv[3]))*1.e6)+5
    Mtotiv = -2.5*np.log10(10**(-0.4*M1iv)+10**(-0.4*M2iv))
    colorv = Mtotv-Mtotiv
                                                                
    L1b = np.loadtxt('../pareconk/Pares-PAUS_W1-'+best).T                                                            
    # L2b = np.loadtxt('../pareconk/Pares-PAUS_W2-'+best).T                                    
    L3b = np.loadtxt('../pareconk/Pares-PAUS_W3-'+best).T
    
    # Lb = np.vstack((L1b.T,L2b.T,L3b.T)).T
    Lb = np.vstack((L1b.T,L3b.T)).T
    
    M1b = Lb[8]-5.*np.log10(np.array(cosmo.luminosity_distance(Lb[3]))*1.e6)+5
    M2b = Lb[13]-5.*np.log10(np.array(cosmo.luminosity_distance(Lb[3]))*1.e6)+5
    Mtotb = -2.5*np.log10(10**(-0.4*M1b)+10**(-0.4*M2b))
    M1ib = Lb[7]-5.*np.log10(np.array(cosmo.luminosity_distance(Lb[3]))*1.e6)+5
    M2ib = Lb[12]-5.*np.log10(np.array(cosmo.luminosity_distance(Lb[3]))*1.e6)+5
    Mtotib = -2.5*np.log10(10**(-0.4*M1ib)+10**(-0.4*M2ib))
    
    rpb = Lb[-2]
    Vrb = Lb[-1]*cvel*1.e-3
    
    
    colorb = Mtotb-Mtotib

    label_x = '$M_{r}$'
    label_y = '$M_{r} - M_{i}$'
    
    Lratio = 10.**(-0.4*(L[-3]-L[8]))
    Lratiob = 10.**(-0.4*(Lb[-3]-Lb[8]))
    
    f, ax = plt.subplots(3,3, figsize=(12,12))
    f.subplots_adjust(hspace=0)
    
    ax = ax.flatten()

    for j in range(len(ax)):
            ax[j].set_ylabel('$n$')
    
    # ax[0].hist(Mtot_mice,20,color='C7',label='MICE sample',lw=2,alpha=0.4,density=True)
    ax[0].plot(-24,0.1,'w,',label='Total sample')
    ax[0].hist(Mr,np.linspace(-24,-17.5,20),color='C7',alpha=0.5,label='PAUS galaxies',lw=2,density=True)
    ax[0].axvline(np.median(Mr),color='C7',lw=2)
    ax[0].hist(M1,np.linspace(-24,-19.5,20),color='k',label='Galaxy pairs',lw=2,histtype='step',density=True)
    ax[0].hist(M1,np.linspace(-24,-19.5,20),color='C4',lw=2,histtype='step',density=True)
    ax[0].hist(M2,np.linspace(-24,-17.5,20),color='C4',lw=2,histtype='step',density=True,ls='--')
    ax[0].axvline(np.median(M1),color='C4',lw=2)
    ax[0].axvline(np.median(M2),color='C4',lw=2,ls='--')
    ax[0].legend(frameon=False,loc=2,fontsize=12)    
    
    ax[3].plot(-24,0.1,'w,',label='Gold sample')
    ax[3].legend(frameon=False,loc=2,fontsize=12)    
    ax[3].set_ylim([0,0.52])
    ax[3].hist(Mrb,np.linspace(-24,-17.5,20),color='C7',alpha=0.5,label='Total sample',lw=2,density=True)
    ax[3].axvline(np.median(Mrb),color='C7',lw=2)
    ax[3].hist(M1b,np.linspace(-24,-19.5,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    ax[3].hist(M2b,np.linspace(-24,-17.5,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True,ls='--')
    ax[3].axvline(np.median(M1b),color='C1',lw=2)
    ax[3].axvline(np.median(M2b),color='C1',lw=2,ls='--')
    ax[3].set_xlabel('$M_{r}$',fontsize=11)

    ax[6].plot(-24,0.1,'w,',label='Deep-z sample')
    ax[6].legend(frameon=False,loc=2,fontsize=12)    
    ax[6].hist(Mrv,np.linspace(-24,-17.5,20),color='C7',alpha=0.5,label='Total sample',lw=2,density=True)
    ax[6].axvline(np.median(Mrv),color='C7',lw=2)
    ax[6].hist(M1v,np.linspace(-24,-19.5,20),color='C2',label='Gold sample',lw=2,histtype='step',density=True)
    ax[6].hist(M2v,np.linspace(-24,-17.5,20),color='C2',label='Gold sample',lw=2,histtype='step',density=True,ls='--')
    ax[6].axvline(np.median(M1v),color='C2',lw=2)
    ax[6].axvline(np.median(M2v),color='C2',lw=2,ls='--')
    ax[6].set_xlabel('$M_{r}$',fontsize=11)
    
    ax[1].axvline(np.median(M1-M1i),color='k',lw=2, label='Brightest')
    ax[1].axvline(np.median(M2-M2i),color='k',lw=2,ls='--', label='Faintest')
    ax[1].legend(frameon=False,loc=2,fontsize=12)    
    ax[1].hist(colorp,np.linspace(-0.,0.6,20),color='C7',alpha=0.5,label='PAUS galaxies',lw=2,density=True)
    ax[1].axvline(np.median(colorp),color='C7',lw=2)
    ax[1].hist(M1-M1i,np.linspace(-0.,0.6,20),color='C4',lw=2,histtype='step',density=True)
    ax[1].hist(M2-M2i,np.linspace(-0.,0.6,20),color='C4',lw=2,histtype='step',density=True,ls='--')
    ax[1].axvline(np.median(M1-M1i),color='C4',lw=2)
    ax[1].axvline(np.median(M2-M2i),color='C4',lw=2,ls='--')
    
    
    ax[4].hist(colorpb,np.linspace(-0.,0.6,20),color='C7',alpha=0.5,label='Gold sample',lw=2,density=True)
    ax[4].axvline(np.median(colorpb),color='C7',lw=2)
    ax[4].hist(M1b-M1ib,np.linspace(-0.,0.6,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    ax[4].hist(M2b-M2ib,np.linspace(-0.,0.6,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True,ls='--')
    ax[4].axvline(np.median(M1b-M1ib),color='C1',lw=2)
    ax[4].axvline(np.median(M2b-M2ib),color='C1',lw=2,ls='--')
    ax[4].set_xlabel('$M_{r}-M_{i}$')

    ax[7].hist(colorpv,np.linspace(-0.,0.6,20),color='C7',alpha=0.5,label='Gold sample',lw=2,density=True)
    ax[7].axvline(np.median(colorpv),color='C7',lw=2)
    ax[7].hist(M1v-M1iv,np.linspace(-0.,0.6,20),color='C2',label='Gold sample',lw=2,histtype='step',density=True)
    ax[7].hist(M2v-M2iv,np.linspace(-0.,0.6,20),color='C2',label='Gold sample',lw=2,histtype='step',density=True,ls='--')
    ax[7].axvline(np.median(M1v-M1iv),color='C2',lw=2)
    ax[7].axvline(np.median(M2v-M2iv),color='C2',lw=2,ls='--')
    ax[7].set_xlabel('$M_{r}-M_{i}$')

    ax[2].hist(zp,np.linspace(0.2,0.6,20),color='C7',label='Total sample',lw=2,alpha=0.5,density=True)
    ax[2].axvline(np.median(zp),color='C7',lw=2)
    ax[2].hist(L[3],np.linspace(0.2,0.6,20),color='C4',label='Total sample',lw=2,histtype='step',density=True)
    ax[2].axvline(np.median(L[3]),color='C4',lw=2)

    ax[5].hist(zpb,np.linspace(0.2,0.6,20),color='C7',label='Total sample',lw=2,alpha=0.5,density=True)
    ax[5].axvline(np.median(zpb),color='C7',lw=2)
    ax[5].hist(Lb[3],np.linspace(0.2,0.6,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    ax[5].axvline(np.median(Lb[3]),color='C1',lw=2)
    ax[5].set_xlabel(r'$z$')

    ax[8].hist(zpv,np.linspace(0.2,0.6,20),color='C7',label='Total sample',lw=2,alpha=0.5,density=True)
    ax[8].axvline(np.median(zpv),color='C7',lw=2)
    ax[8].hist(Lv[3],np.linspace(0.2,0.6,20),color='C2',label='Gold sample',lw=2,histtype='step',density=True)
    ax[8].axvline(np.median(Lv[3]),color='C2',lw=2)
    ax[8].set_xlabel(r'$z$')
    
    f.savefig('../Mdist_withdeep.pdf',bbox_inches='tight')


    # separate_medianas(Mtotb[colorb>0],colorb[colorb>0],label_x = label_x, label_y = label_y, out_plot = '../final_plots/color_mag_gold.pdf')
    # separate_medianas(Mtotv[colorv>0],colorv[colorv>0],label_x = label_x, label_y = label_y, out_plot = '../color_mag_withdeep.pdf')


def color_plot():
    
    from medianas import separate_medianas
    
    # LOAD PAUS CATALOGS
    
    w1 = pd.read_csv('../pareconk/Photo_z_calibrate_iterative_W1_no_k_corr.csv')
    w3 = pd.read_csv('../pareconk/Photo_z_calibrate_iterative_W3_no_k_corr.csv')
    w1_best = pd.read_csv('../pareconk/Photo_z_calibrate_iterative_W1_best_pz_no_k_corr.csv')
    w3_best = pd.read_csv('../pareconk/Photo_z_calibrate_iterative_W3_best_pz_no_k_corr.csv')
    
    mr  = np.append(w1.mag_r_no_k,w3.mag_r_no_k)
    mrb = np.append(w1_best.mag_r_no_k,w3_best.mag_r_no_k)
    zp  = np.append(w1.zb,w3.zb)
    zpb  = np.append(w1_best.zb,w3_best.zb)

    mp  = (mr > -23.5) & (zp < 0.6) & (zp > 0.2)
    mpb = (mrb > -23.5) & (zpb < 0.6) & (zpb > 0.2)
    
    mr  = mr[mp]
    zp  = zp[mp]
    zpb = zpb[mpb]
    
    mrb  = mrb[mpb]
    mi  = np.append(w1.mag_i_no_k,w3.mag_i_no_k)[mp]
    
    mib = np.append(w1_best.mag_i_no_k,w3_best.mag_i_no_k)[mpb]
    
    Mr  = mr - 5.*np.log10(np.array(cosmo.luminosity_distance(zp))*1.e6) + 5
    Mi  = mi - 5.*np.log10(np.array(cosmo.luminosity_distance(zp))*1.e6) + 5
    Mrb = mrb - 5.*np.log10(np.array(cosmo.luminosity_distance(zpb))*1.e6) + 5
    Mib = mib - 5.*np.log10(np.array(cosmo.luminosity_distance(zpb))*1.e6) + 5

    colorp  = Mr - Mi
    colorpb = Mrb - Mib
        
    
    # LOAD MICE CATALOG
    
    M1 = np.loadtxt('../MICE-pairs/True-pairs-0.2-0.6-50-paralentes-R1.dat').T
    M2 = np.loadtxt('../MICE-pairs/True-pairs-0.2-0.6-50-paralentes-R2.dat').T
    M3 = np.loadtxt('../MICE-pairs/True-pairs-0.2-0.6-50-paralentes-R3.dat').T
    M4 = np.loadtxt('../MICE-pairs/True-pairs-0.2-0.6-50-paralentes-R4.dat').T
    
    M = np.vstack((M1.T,M2.T,M3.T,M4.T)).T
    
    M1_mice = np.zeros(len(M[1]))
    M2_mice = np.zeros(len(M[1]))
    
    M1_mice[M[8]<M[17]] = M[8][M[8]<M[17]]
    M1_mice[M[8]>M[17]] = M[17][M[8]>M[17]]
    M2_mice[M[8]<M[17]] = M[17][M[8]<M[17]]
    M2_mice[M[8]>M[17]] = M[8][M[8]>M[17]]

    print('density M1 = ',compute_density(M1[1],M1[2]))
    print('density M2 = ',compute_density(M2[1],M2[2]))
    print('density M3 = ',compute_density(M3[1],M3[2]))
    print('density M4 = ',compute_density(M4[1],M4[2]))

    
    Mtot_mice   = -2.5*np.log10(10**(-0.4*M[8])+10**(-0.4*M[17]))
    Mtot_i_mice = -2.5*np.log10(10**(-0.4*M[10])+10**(-0.4*M[19]))
    z_mice      = M[3]
    Lratio_mice = 10.**(-0.4*(M2_mice-M1_mice))
    color_mice  = Mtot_mice - Mtot_i_mice
    
    L1 = np.loadtxt('../pareconk/Pares-PAUS_W1-'+pcat).T                                                            
    # L2 = np.loadtxt('../pareconk/Pares-PAUS_W2-'+pcat).T                                    
    L3 = np.loadtxt('../pareconk/Pares-PAUS_W3-'+pcat).T
    
    print('density W1 = ',compute_density(L1[1],L1[2]))
    # print('density W2 = ',compute_density(L2[1],L2[2]))
    print('density W3 = ',compute_density(L3[1],L3[2]))
    
    # L = np.vstack((L1.T,L2.T,L3.T)).T
    L = np.vstack((L1.T,L3.T)).T
    
    M1 = L[8]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    M2 = L[13]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    Mtot = -2.5*np.log10(10**(-0.4*M1)+10**(-0.4*M2))
    M1i = L[7]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    M2i = L[12]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    Mtoti = -2.5*np.log10(10**(-0.4*M1i)+10**(-0.4*M2i))
    
    rp = L[-2]
    Vr = L[-1]*cvel*1.e-3
    
    color = Mtot-Mtoti

    L1b = np.loadtxt('../pareconk/Pares-PAUS_W1-'+best).T                                                            
    # L2b = np.loadtxt('../pareconk/Pares-PAUS_W2-'+best).T                                    
    L3b = np.loadtxt('../pareconk/Pares-PAUS_W3-'+best).T
    
    # Lb = np.vstack((L1b.T,L2b.T,L3b.T)).T
    Lb = np.vstack((L1b.T,L3b.T)).T
    
    M1b = Lb[8]-5.*np.log10(np.array(cosmo.luminosity_distance(Lb[3]))*1.e6)+5
    M2b = Lb[13]-5.*np.log10(np.array(cosmo.luminosity_distance(Lb[3]))*1.e6)+5
    Mtotb = -2.5*np.log10(10**(-0.4*M1b)+10**(-0.4*M2b))
    M1ib = Lb[7]-5.*np.log10(np.array(cosmo.luminosity_distance(Lb[3]))*1.e6)+5
    M2ib = Lb[12]-5.*np.log10(np.array(cosmo.luminosity_distance(Lb[3]))*1.e6)+5
    Mtotib = -2.5*np.log10(10**(-0.4*M1ib)+10**(-0.4*M2ib))
    
    rpb = Lb[-2]
    Vrb = Lb[-1]*cvel*1.e-3
    
    
    colorb = Mtotb-Mtotib

    
    Lratio = 10.**(-0.4*(L[-3]-L[8]))
    Lratiob = 10.**(-0.4*(Lb[-3]-Lb[8]))
    
    f, ax = plt.subplots(2,3, figsize=(12,6))
    f.subplots_adjust(hspace=0,wspace=0)
    
    ax = ax.flatten()
    

    ax[0].set_ylabel('$n$')
    ax[3].set_ylabel('$n$')
    
    # ax[0].hist(Mtot_mice,20,color='C7',label='MICE sample',lw=2,alpha=0.4,density=True)
    ax[0].plot(-24,0.1,'w,',label='Total sample')
    ax[0].hist(Mr,np.linspace(-24,-19.5,20),color='C7',alpha=0.5,label='PAUS galaxies',lw=2,density=True)
    ax[0].axvline(np.median(Mr),color='C7',lw=2)
    ax[0].hist(M1,np.linspace(-24,-19.5,20),color='k',label='Galaxy pairs',lw=2,histtype='step',density=True)
    ax[0].hist(M1,np.linspace(-24,-19.5,20),color='C4',lw=2,histtype='step',density=True)
    ax[0].hist(M2,np.linspace(-24,-19.5,20),color='C4',lw=2,histtype='step',density=True,ls='--')
    ax[0].axvline(np.median(M1),color='C4',lw=2)
    ax[0].axvline(np.median(M2),color='C4',lw=2,ls='--')
    ax[0].legend(frameon=False,loc=2,fontsize=12)    
    
    ax[3].plot(-24,0.1,'w,',label='Gold sample')
    ax[3].legend(frameon=False,loc=2,fontsize=12)    
    ax[3].set_ylim([0,0.52])
    ax[3].hist(Mrb,np.linspace(-24,-19.5,20),color='C7',alpha=0.5,label='Total sample',lw=2,density=True)
    ax[3].axvline(np.median(Mrb),color='C7',lw=2)
    ax[3].hist(M1b,np.linspace(-24,-19.5,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    ax[3].hist(M2b,np.linspace(-24,-19.5,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True,ls='--')
    ax[3].axvline(np.median(M1b),color='C1',lw=2)
    ax[3].axvline(np.median(M2b),color='C1',lw=2,ls='--')
    ax[3].set_xlabel('$M_{r}$',fontsize=11)
    
    ax[1].axvline(np.median(M1-M1i),color='k',lw=2, label='Brightest')
    ax[1].axvline(np.median(M2-M2i),color='k',lw=2,ls='--', label='Faintest')
    ax[1].legend(frameon=False,loc=2,fontsize=12)    
    ax[1].hist(colorp,np.linspace(-0.,0.6,20),color='C7',alpha=0.5,label='PAUS galaxies',lw=2,density=True)
    ax[1].axvline(np.median(colorp),color='C7',lw=2)
    ax[1].hist(M1-M1i,np.linspace(-0.,0.6,20),color='C4',lw=2,histtype='step',density=True)
    ax[1].hist(M2-M2i,np.linspace(-0.,0.6,20),color='C4',lw=2,histtype='step',density=True,ls='--')
    ax[1].axvline(np.median(M1-M1i),color='C4',lw=2)
    ax[1].axvline(np.median(M2-M2i),color='C4',lw=2,ls='--')
    
    
    ax[4].hist(colorpb,np.linspace(-0.,0.6,20),color='C7',alpha=0.5,label='Gold sample',lw=2,density=True)
    ax[4].axvline(np.median(colorpb),color='C7',lw=2)
    ax[4].hist(M1b-M1ib,np.linspace(-0.,0.6,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    ax[4].hist(M2b-M2ib,np.linspace(-0.,0.6,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True,ls='--')
    ax[4].axvline(np.median(M1b-M1ib),color='C1',lw=2)
    ax[4].axvline(np.median(M2b-M2ib),color='C1',lw=2,ls='--')
    ax[4].set_xlabel('$M_{r}-M_{i}$')

    ax[2].hist(zp,np.linspace(0.2,0.6,20),color='C7',label='Total sample',lw=2,alpha=0.5,density=True)
    ax[2].axvline(np.median(zp),color='C7',lw=2)
    ax[2].hist(L[3],np.linspace(0.2,0.6,20),color='C4',label='Total sample',lw=2,histtype='step',density=True)
    ax[2].axvline(np.median(L[3]),color='C4',lw=2)

    ax[5].hist(zpb,np.linspace(0.2,0.6,20),color='C7',label='Total sample',lw=2,alpha=0.5,density=True)
    ax[5].axvline(np.median(zpb),color='C7',lw=2)
    ax[5].hist(Lb[3],np.linspace(0.2,0.6,20),color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    ax[5].axvline(np.median(Lb[3]),color='C1',lw=2)
    ax[5].set_xlabel(r'$z$')
    
    for x in ax:
        # x.set_yticklabels([])
        x.set_yticks([])
    
    f.savefig('../final_plots/Mdist.pdf',bbox_inches='tight')

    f, ax = plt.subplots(1,3, figsize=(13,3))
    f.subplots_adjust(wspace=0)
    ax = ax.flatten()

    for x in ax:
        # x.set_yticklabels([])
        x.set_yticks([])
    
    
    ax[0].hist(rp,20,color='C4',label='Total sample',lw=2,histtype='step',density=True)
    ax[0].hist(rpb,20,color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    ax[0].axvline(np.median(rp),color='C4',lw=2)
    ax[0].axvline(np.median(rpb),color='C1',lw=2)

    # ax[1].hist(Vr,20,color='C4',label='Total sample',lw=2,histtype='step',density=True)
    # ax[1].hist(Vrb,20,color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    # ax[1].axvline(np.median(Vr),color='C4',lw=2)
    # ax[1].axvline(np.median(Vrb),color='C1',lw=2)

    ax[1].hist(Mtot,20,color='C4',label='Total sample',lw=2,histtype='step',density=True)
    ax[1].hist(Mtotb,20,color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    ax[1].axvline(np.median(Mtot),color='C4',lw=2)
    ax[1].axvline(np.median(Mtotb),color='C1',lw=2)

    ax[2].hist(Lratio,20,color='C4',label='Total sample',lw=2,histtype='step',density=True)
    ax[2].hist(Lratiob,20,color='C1',label='Gold sample',lw=2,histtype='step',density=True)
    ax[2].axvline(np.median(Lratio),color='C4',lw=2)
    ax[2].axvline(np.median(Lratiob),color='C1',lw=2)
    ax[0].set_xlabel('$r_p [kpc]$')
    ax[1].set_xlabel('$M^\mathrm{pair}_r$')
    ax[2].set_xlabel('$L_2/L_1$')
    ax[0].set_ylabel('$n$')
    ax[2].legend(frameon=True,loc=1,fontsize=12)

    f.savefig('../final_plots/Lratio.pdf',bbox_inches='tight')

    label_x = r'$M^\mathrm{pair}_{r}$'
    label_y = r'$M^\mathrm{pair}_{r} - M^\mathrm{pair}_{i}$'

    # separate_medianas(Mtotb[colorb>0],colorb[colorb>0],label_x = label_x, label_y = label_y, out_plot = '../final_plots/color_mag_gold.pdf')
    separate_medianas(Mtot[color>0],color[color>0],label_x = label_x, label_y = label_y, out_plot = '../final_plots/color_mag.pdf')


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
                       axDS = plt,
                       RIN=300,ROUT=10000, fytpe = '',
                       ylabel = True, chist = 'C0',plot = True):
    

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    print(h['N_LENSES'])
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
    Mmean = h['M_mean_r']    
    
    ndots = p.shape[0]
    
    
    # FIT MONOPOLE
    fitted = fits.open(folder+'fitresults'+ftype+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)
    fitpar = fitted[0].header
    lgM   = fitted[1].data.logM
    
    rplot = np.logspace(np.log10((h['RIN']-20.)/1000.),np.log10((h['ROUT'])/1000.),20)
    
    nfw  = Delta_Sigma_fit(p.Rp,p.DSigma_T,p.error_DSigma_T,zmean,cosmo_as)
    
    fmass = np.percentile(lgM[2500:], [16,50,84])
    
    mass = str(np.round(fmass[1],1))+'^{+'+str(np.round(np.diff(fmass)[1],1))+'}'+'_{-'+str(np.round(np.diff(fmass)[0],1))+'}'
    cfit = str(np.round(fitpar['c200'],1))+'$'
    
    
    if plot:
        
        ds1h   = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='1h')    
        ds2h   = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='2h')    
        #ds2h   = ds1h
        
        
        ds = ds1h+ds2h
        
        maskr = (p.Rp > (RIN/1000.))*(p.Rp < (ROUT/1000.))
        dsfit = Delta_Sigma_NFW_2h(p.Rp[maskr],zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='1h')    
        
        chi2 = str(np.round(chi_red(dsfit,p.DSigma_T[maskr],p.error_DSigma_T[maskr],1),1))
        
        axDS.plot(0,0,'w.',label=lsamp)
        axDS.plot(p.Rp,p.DSigma_T,'C1o')
        axDS.errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,ecolor='C1',fmt='None')
        axDS.plot(rplot,ds,'C3',label='$\log M_{200}= '+mass+'\,\,c_{200} = '+cfit)
        axDS.plot(0,0,'w.',label=r'$\chi^2_{red} = $'+chi2)
        axDS.plot(rplot,ds1h,'C4')
        axDS.plot(rplot,ds2h,'C4--')
        axDS.set_xscale('log')
        axDS.set_yscale('log')
        if ylabel:
            axDS.set_ylabel(r'$\Delta\Sigma_{T} [M_{\odot}pc^{-2} h ]$')
        axDS.set_xlabel(r'$R [Mpc/h]$')
        axDS.set_ylim(0.095,1500)
        axDS.set_xlim((h['RIN']-20)/1000.,(h['ROUT']+3e3)/1000.)
        axDS.yaxis.set_ticks([0.1,1,10,100])
        axDS.set_yticklabels([0.1,1,10,100])
        axDS.xaxis.set_ticks([0.1,1,10])
        axDS.set_xticklabels([0.1,1,10])
        axDS.axvline(RIN/1000.,color='C7',alpha=0.5,ls='--')
        axDS.legend(frameon=False,loc=1,fontsize=10)
        
    return [fmass,Mmean,lgM[2500:]]


def plt_mcmc(samp,lsamp,axC = plt,axC2 = plt,
             RIN=300,ROUT=10000, fytpe = '',
            ylabel = True, legend=True,
            lw=1,chist = 'C0'):
    

    p_name = 'profile_'+samp+'.fits'
    
    # FIT MONOPOLE
    fitted = fits.open(folder+'fitresults'+ftype+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)
    fitpar = fitted[0].header
    lgM   = fitted[1].data.logM
    
    axC2.plot(lgM,chist,label=lsamp,lw=lw,alpha=0.2)
    axC.hist(lgM[2500:],np.linspace(10.7,13.5,50),histtype='step',color=chist,label=lsamp,lw=lw)
    axC.axvline(np.median(lgM[2500:]),alpha=1.,color=chist)
    # axC.axvline(np.percentile(lgM[2500:], [16,50,84])[0],ls='--',color=chist)
    # axC.axvline(np.percentile(lgM[2500:], [16,50,84])[2],ls='--',color=chist)
    
    axC.axvspan(np.percentile(lgM[2500:], [16,50,84])[0],np.percentile(lgM[2500:], [16,50,84])[2],alpha=0.1,color=chist)
    
    axC.set_xlabel('$\log M_{200}$')
    
    axC2.set_xlabel('$N$')
    if legend:
        axC.legend(loc=2,fontsize=12)
        axC2.legend(loc=2,fontsize=12)
    if ylabel:
        axC.set_ylabel('$P(\log M_{200})$')   
        axC2.set_ylabel('$\log M_{200}$')


def dilution(samp):
    

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data

    bines = np.logspace(np.log10(h['RIN']/1000.),np.log10(h['ROUT']/1000.),num=len(p)+1)
    area = np.pi*np.diff(bines**2)
    
    ngal = p.NGAL_w

    d = ngal/area

    fcl = ((d - np.mean(d[-3:-1]))*area)/ngal
    
    # plt.plot(p.Rp,1./(1-fcl))
    # plt.plot(p.Rp,1./(1-fcl),'o')
    
    return p.Rp,d,fcl,h['N_LENSES']
    
def fcl_plot(samples,lsamps,csamps,marker,ax=plot):
    
    for j in range(len(samples)):
    
        samp = samples[j]
      
        r,d,fcl,n = dilution(samp)
        if j == 0:
            ax.plot(r,1./(1-fcl),marker[j],color=csamps[j],label=lsamps[j],markersize=10)
        else:
            ax.plot(r,1./(1-fcl),marker[j],color=csamps[j],label=lsamps[j],markersize=6)
        ax.plot(r,1./(1-fcl),csamps[j],alpha=0.5)

    ax.set_xscale('log')
    ax.set_xlabel('$R [Mpc/h]$')
    ax.set_ylabel('$(1-f_{s})^{-1}$')


# pcats = ['_zspec',
        # '_zspec_best',
        # '_photo_z_2nd_run_mag_i',
        # '_photo_z_2nd_run_mag_i_best']

def make_plot_profile():
    
    lMfit = []
    
    fDS, axDS = plt.subplots(5,4, figsize=(14.5,17),sharex = True,sharey = True)
    fDS.subplots_adjust(hspace=0,wspace=0)
    
    fC, axC = plt.subplots(5,2, figsize=(6,12),sharex = True,sharey = True)
    fC.subplots_adjust(hspace=0,wspace=0)
    fC2, axC2 = plt.subplots(5,2, figsize=(14,12),sharex = True,sharey = True)
    fC2.subplots_adjust(hspace=0,wspace=0)
    
    axC[0,1].axis('off')
    axC2[0,1].axis('off')
    axDS[0,2].axis('off')
    axDS[0,3].axis('off')
    
    axDS = axDS.flatten()
    axC = axC.flatten()
    axC2 = axC2.flatten()
    
    axDS = np.append(axDS[:2],axDS[4:])
    axC  = np.append(axC[:1],axC[2:])
    axC2 = np.append(axC2[:1],axC2[2:])
    
    
            
    samp =  ['pk_w1w3_'+pcat,'pk_w1w3_'+best,
            'pk_Mm2_w1w3_'+pcat,'pk_Mm2_w1w3_'+best,
            'pk_MM2_w1w3_'+pcat,'pk_MM2_w1w3_'+best,
            'pk_zm_w1w3_'+pcat,'pk_zm_w1w3_'+best,
            'pk_zM_w1w3_'+pcat,'pk_zM_w1w3_'+best,
            'pk_Lrm2_w1w3_'+pcat,'pk_Lrm2_w1w3_'+best,
            'pk_LrM2_w1w3_'+pcat,'pk_LrM2_w1w3_'+best,
            'pk_red_w1w3_'+pcat,'pk_red_w1w3_'+best,
            'pk_blue_w1w3_'+pcat,'pk_blue_w1w3_'+best]
    
    lsamp1 = ['Total sample','Gold sample']*len(samp)
    
    lsamp2 = ['all pairs','all pairs', 
              r'$M^\mathrm{pair}_r < -22.5$',r'$M^\mathrm{pair}_r < -22.5$', 
              r'$M^\mathrm{pair}_r \geq -22.5$',r'$M^\mathrm{pair}_r \geq -22.5$', 
              '$z < 0.4$','$z < 0.4$', 
              '$z \geq 0.4$','$z \geq 0.4$',
              '$L_2/L_1 < 0.5$','$L_2/L_1 < 0.5$',
              '$L_2/L_1 \geq 0.5$','$L_2/L_1 \geq 0.5$',
              r'$red\,\,pairs$',r'$red\,\,pairs$',
              r'$blue\,\,pairs$',r'$blue\,\,pairs$']
              
    lsamp = ['Total sample','Gold sample']+lsamp2[2:]
            
    ind = [0,0,1,2,1,2,3,4,3,4,5,6,5,6,7,8,7,8]
    chist = ['C4','C1']*9
    lw = [1,1]+[1,1,2,2]*4
    
    ylabel = [True,False]+[True,False,False,False]*4
    legend = [True,True]+[True,False,True,False]*4
    
    for j in range(len(samp)):
        lMfit += [plt_profile_fit_2h(samp[j],
                  lsamp1[j]+' - '+lsamp2[j],axDS[j],
                  fytpe = ftype, ylabel=ylabel[j],chist=chist[j])]
                
        plt_mcmc(samp[j],lsamp[j],axC[ind[j]], axC2[ind[j]],
                fytpe = ftype, ylabel=ylabel[j],
                legend=legend[j],chist=chist[j],lw=lw[j])
        # lMfit += [plt_profile_fit_2h(samp[j],lsamp[j],plot=False,fytpe = ftype)]
    
    
    axC[0].set_ylim(0,4500)
    axC[0].set_xlim(10.7,13.5)
    axC2[0].set_ylim(10.7,13.5)
    fDS.savefig('../final_plots/profile.pdf',bbox_inches='tight')
    fC.savefig('../final_plots/chains2'+pcat+ftype+'.pdf',bbox_inches='tight')


def make_fcl_plot():

    samples =  ['pk_w1w3_'+pcat,
            'pk_Mm2_w1w3_'+pcat,'pk_MM2_w1w3_'+pcat,
            'pk_zm_w1w3_'+pcat,'pk_zM_w1w3_'+pcat,
            'pk_Lrm2_w1w3_'+pcat,'pk_LrM2_w1w3_'+pcat,
            'pk_blue_w1w3_'+pcat,'pk_red_w1w3_'+pcat,]

    samples_gold =  ['pk_w1w3_'+best,'pk_Mm2_w1w3_'+best,
            'pk_MM2_w1w3_'+best,'pk_zm_w1w3_'+best,
            'pk_zM_w1w3_'+best,'pk_Lrm2_w1w3_'+best,
            'pk_LrM2_w1w3_'+best,
            'pk_blue_w1w3_'+best,'pk_red_w1w3_'+best]
    
    
    csamp = ['k',
            'gold','gold',
            'royalblue','royalblue',
            'C9','C9',
            'palevioletred','palevioletred',]
    
    lsamp = ['all pairs',
            r'$M^\mathrm{pair}_r < -22.5$',
            r'$M^\mathrm{pair}_r \geq -22.5$',
            '$z < 0.4$',
            '$z \geq 0.4$',
            '$L_2/L_1 < 0.5$',
            '$L_2/L_1 \geq 0.5$',
            r'$blue\,\,pairs$',
            r'$red\,\,pairs$']

    mark = ['o']+['v','^']*4

    f, ax = plt.subplots(2,1, figsize=(6,8),sharex = True,sharey = True)
    f.subplots_adjust(hspace=0,wspace=0)

    ax[0].text(0.13,1.27,'Total sample') 
    ax[1].text(0.13,1.27,'Gold sample') 
    ax[1].legend(frameon=False,loc=1)
    fcl_plot(samples,lsamp,csamp,mark,ax=ax[0])
    fcl_plot(samples_gold,lsamp,csamp,mark,ax=ax[1])
    
    ax[0].legend(frameon=False,ncol=2,loc=1,fontsize=11)
    ax[0].set_ylim([0.97,1.3])
    ax[0].grid(True, linestyle=":", color='C7')
    ax[1].grid(True, linestyle=":", color='C7')
    f.savefig('../final_plots/contamination.pdf',bbox_inches='tight')

def make_mag_mass_plot():
    
    def linear(x,m,n):
        return x*m+n
    
    
    from scipy.optimize import curve_fit

    # MICE MASS-TO-LIGHT RELATION
    
    popt, pcov = curve_fit(linear, meanmag, lM200)
    m_mice,n_mice = popt
    print('MICE sample fit')
    print(m_mice,n_mice)    

    Lsun    = 10**(-0.4*4.65)
    LViola  = ((10**(-0.4*meanmag)) / (Lsun*10**11.5))
    lMviola = np.log10(((0.95*LViola)**1.16)*1.e14)

    m_v = -0.4*1.16
    em_v = np.abs(-0.4*0.13)
    n_v = np.log10(0.95)+14.-1.16*np.log10(Lsun*10**11.5)
    en_v = np.abs(0.14/(0.95*np.log(10.))) + np.abs(0.13*np.log10(Lsun*10**11.5))
    print('Viola et. al.')
    print(m_v,n_v)    
    print(em_v,en_v)    

    m_mice,n_mice = m_v,n_v
    
    # SAMPLES AND LABELS
    
    samples =  ['pk_w1w3_'+pcat,
            'pk_Mm2_w1w3_'+pcat,'pk_MM2_w1w3'+pcat,
            'pk_zm_w1w3_'+pcat,'pk_zM_w1w3_'+pcat,
            'pk_Lrm2_w1w3_'+pcat,'pk_LrM2_w1w3_'+pcat,
            'pk_blue_w1w3_'+pcat,'pk_red_w1w3_'+pcat]

    samples_gold =  ['pk_w1w3_'+best,
            'pk_Mm2_w1w3_'+best,'pk_MM2_w1w3'+best,
            'pk_zm_w1w3_'+best,'pk_zM_w1w3_'+best,
            'pk_Lrm2_w1w3_'+best,'pk_LrM2_w1w3_'+best,
            'pk_blue_w1w3_'+best,'pk_red_w1w3_'+best]
    
    
    csamp = ['k',
            'gold','gold',
            'royalblue','royalblue',
            'C9','C9',
            'palevioletred','palevioletred']
    
    lsamp = ['all pairs',
            r'$M^{pair}_r < -22.5$',
            r'$M^\mathrm{pair}_r \geq -22.5$',
            '$z < 0.4$',
            '$z \geq 0.4$',
            '$L_2/L_1 < 0.5$',
            '$L_2/L_1 \geq 0.5$',
            r'$blue\,\,pairs$',
            r'$red\,\,pairs$']

    mark = ['o'] + ['v','^']*4
    
    # EXTRACTING DATA
    
    lMfit = []
    lMfit_gold = []
    
    for j in range(len(samples)):
        lMfit += [plt_profile_fit_2h(samples[j],lsamp[j],plot=False,fytpe = ftype)]
        lMfit_gold += [plt_profile_fit_2h(samples_gold[j],lsamp[j],plot=False,fytpe = ftype)]


    # MAKE PLOTS

    fig, axall = plt.subplots(2,1, figsize=(11,10),sharex = True,gridspec_kw={'height_ratios': [2.5, 1]})
    fig.subplots_adjust(hspace=0,wspace=0)

    ax  = axall[0]
    axr = axall[1]


    ax.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt='k',marker='o',label='Total sample')
    ax.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt='k',marker='o',label='Gold sample', mfc='w',mew=3)
    axr.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt='k',marker='o',label='Total sample')
    axr.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt='k',marker='o',label='Gold sample', mfc='w',mew=3)

    
    
    MM      = []
    LM      = []
    eLM      = []
    MM_gold = []
    LM_gold = []
    eLM_gold = []
    
    for j in np.arange(len(csamp)):
        
        MM       += [lMfit[j][1]]
        LM       += [lMfit[j][0][1]]
        eLM      += [np.diff(lMfit[j][0]).mean()]
        MM_gold  += [lMfit_gold[j][1]]
        LM_gold  += [lMfit_gold[j][0][1]]
        eLM_gold += [np.diff(lMfit_gold[j][0]).mean()]
        
        
        # vp = ax.violinplot(lMfit[j][2],
                              # positions=[lMfit[j][1]],
                              # showextrema=False, 
                              # showmedians=False)
        # vp['bodies'][0].set_facecolor(csamp[j])
        # vp['bodies'][0].set_alpha(0.1)
        # vp = ax.violinplot(lMfit_gold[j][2],
                              # positions=[lMfit_gold[j][1]],
                              # showextrema=False, 
                              # showmedians=False)
        # vp['bodies'][0].set_facecolor(csamp[j])
        # vp['bodies'][0].set_alpha(0.1)


        # MAKE SCATTER PLOT

        ax.errorbar(lMfit[j][1],lMfit[j][0][1],markersize=15,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt=csamp[j],marker=mark[j])
        ax.errorbar(lMfit_gold[j][1],lMfit_gold[j][0][1],markersize=15,
                     yerr=np.array([np.diff(lMfit_gold[j][0])]).T,
                     fmt=csamp[j],marker=mark[j], mfc='w',mew=3)
                     

        # MAKE MICE RATIO PLOT

        err = 10**(lMfit[j][0][1]-(lMfit[j][1]*m_mice+n_mice))*np.log(10.)*np.array([np.diff(lMfit[j][0])]).T
        axr.errorbar(lMfit[j][1],10**(lMfit[j][0][1]-(lMfit[j][1]*m_mice+n_mice)),
                        yerr=err,markersize=10,
                        fmt=csamp[j],marker=mark[j])
        
        err = 10**(lMfit_gold[j][0][1]-(lMfit_gold[j][1]*m_mice+n_mice))*np.log(10.)*np.array([np.diff(lMfit_gold[j][0])]).T
        
        axr.errorbar(lMfit_gold[j][1],10**(lMfit_gold[j][0][1]-(lMfit_gold[j][1]*m_mice+n_mice)),
                        yerr=err,markersize=10,
                        fmt=csamp[j],marker=mark[j],mfc='w',mew=3)


        # JUST FOR THE LABELS

        if j > 0:
            ax.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt=csamp[j],marker=mark[j],label=lsamp[j])
            axr.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt=csamp[j],marker=mark[j],label=lsamp[j])
                     
        # if j == 6:
            # ax.errorbar(0,0,markersize=5,
                     # yerr=np.array([np.diff(lMfit[j][0])]).T,
                     # fmt='w',marker=',',label=' ')
    


    MM = np.array(MM)
    MM_gold = np.array(MM_gold)
    LM_gold = np.array(LM_gold)
    eLM_gold = np.array(eLM_gold)
    LM = np.array(LM)
    eLM = np.array(eLM)
    
    
    # FIT MASS-TO-LIGHT RELATION

    mlim = -20.
    
    mask = (MM < mlim)
    popt, pcov = curve_fit(linear, MM[mask], LM[mask],sigma=eLM[mask])
    m,n = popt
    em,en = np.abs(np.sqrt(np.diag(pcov)))
    print('Total sample fit')
    print(m,n)
    print(em,en)

    mask = MM_gold < mlim
    popt, pcov = curve_fit(linear, MM_gold[mask], LM_gold[mask],sigma=eLM_gold[mask])
    m_gold,n_gold = popt
    em_gold,en_gold = np.abs(np.sqrt(np.diag(pcov)))
    print('Gold sample fit')
    print(m_gold,n_gold)
    print(em_gold,en_gold)
    
    mask = np.append(MM,MM_gold) < mlim
    popt, pcov = curve_fit(linear, np.append(MM,MM_gold)[mask], np.append(LM,LM_gold)[mask],sigma=np.append(eLM,eLM_gold)[mask])
    m_all,n_all = popt
    em_all,en_all = np.abs(np.sqrt(np.diag(pcov)))
    print('All sample fit')
    print(m_all,n_all)    
    print(em_all,en_all)    

    ax.plot(meanmag,meanmag*m_v+n_v,'C8',
            label=r'Viola et. al. (2015), $\alpha = $'+str(np.round(m_v,2))+'$\pm$'+str(np.round(em_v,2))+r', $\beta = $'+str(np.round(n_v,1))+'$\pm$'+str(np.round(en_v,1)))
            
    ax.fill_between(meanmag, 
                    linear(meanmag,m_v+em_v,n_v+en_v), 
                    linear(meanmag,m_v-em_v,n_v-en_v), 
                    interpolate=True, 
                    color='C8',alpha=0.3)

    
    ax.plot(meanmag,meanmag*m+n,'C4',
            label=r'Total sample, $\alpha = $'+str(np.round(m,2))+'$\pm$'+str(np.round(em,2))+r', $\beta = $'+str(np.round(n,1))+'$\pm$'+str(np.round(en,1)))
            
    ax.fill_between(meanmag, 
                    linear(meanmag,m+em,n+en), 
                    linear(meanmag,m-em,n-en), 
                    interpolate=True, 
                    color='C4',alpha=0.3)

    ax.plot(meanmag,meanmag*m_gold+n_gold,'C1',
            label=r'Gold sample, $\alpha = $'+str(np.round(m_gold,2))+'$\pm$'+str(np.round(em_gold,2))+r', $\beta = $'+str(np.round(n_gold,1))+'$\pm$'+str(np.round(en_gold,1)))
            
    ax.fill_between(meanmag, 
                    linear(meanmag,m_gold+em_gold,n_gold+en_gold), 
                    linear(meanmag,m_gold-em_gold,n_gold-en_gold), 
                    interpolate=True, 
                    color='C1',alpha=0.3)

    ax.plot(meanmag,meanmag*m_all+n_all,'C3',
            label=r'All sample, $\alpha = $'+str(np.round(m_all,2))+'$\pm$'+str(np.round(em_all,2))+r', $\beta = $'+str(np.round(n_gold,1))+'$\pm$'+str(np.round(en_all,1)))
            
    ax.fill_between(meanmag, 
                    linear(meanmag,m_all+em_all,n_all+en_all), 
                    linear(meanmag,m_all-em_all,n_all-en_all), 
                    interpolate=True, 
                    color='C3',alpha=0.3)

    ax.plot(meanmag,lM200,'C7',label='MICE - True pairs')
    # ax.plot(0,0,'w,',label=' ')
    
    ax.legend(frameon=False,loc=3,ncol=3,fontsize=11)
    ax.set_xlabel(r'$\langle M_r \rangle$')
    ax.set_ylabel(r'$\log (M_{200}/M_\odot h^{-1})$')
    ax.axis([-21.7,-20.4,10.8,13.2])

    # f = lambda x: 10**(linear(x,m+em,n+en) - linear(x,m_mice,n_mice))
    # g = lambda x: 10**(linear(x,m-em,n-en) - linear(x,m_mice,n_mice))
    # axr.plot(meanmag,10**(linear(meanmag,m,n) - linear(meanmag,m_mice,n_mice)),'C4',label='Total sample')
    # axr.fill_between(meanmag, f(meanmag), g(meanmag), where=f(meanmag)<=g(meanmag), interpolate=True, color='C4',alpha=0.3)
    # axr.fill_between(meanmag, f(meanmag), g(meanmag), where=f(meanmag)>=g(meanmag), interpolate=True, color='C4',alpha=0.3)

    # f = lambda x: 10**(linear(x,m_gold+em_gold,n_gold+en_gold) - linear(x,m_mice,n_mice))
    # g = lambda x: 10**(linear(x,m_gold-em_gold,n_gold-en_gold) - linear(x,m_mice,n_mice))   
    # axr.plot(meanmag,10**((meanmag*m_gold+n_gold) - (meanmag*m_mice+n_mice)),'C1',label='Gold sample')
    # axr.fill_between(meanmag, f(meanmag), g(meanmag), where=f(meanmag)<=g(meanmag), interpolate=True, color='C1',alpha=0.3)
    # axr.fill_between(meanmag, f(meanmag), g(meanmag), where=f(meanmag)>=g(meanmag), interpolate=True, color='C1',alpha=0.3)

    # f = lambda x: 10**(linear(x,m_all+em_all,n_all+en_all) - linear(x,m_mice,n_mice))
    # g = lambda x: 10**(linear(x,m_all-em_all,n_all-en_all) - linear(x,m_mice,n_mice))    
    # axr.plot(meanmag,10**(linear(meanmag,m_all,n_all) - linear(meanmag,m_mice,n_mice)),'C3',label='All sample')
    # axr.fill_between(meanmag, f(meanmag), g(meanmag), where=f(meanmag)<=g(meanmag), interpolate=True, color='C3',alpha=0.3)
    # axr.fill_between(meanmag, f(meanmag), g(meanmag), where=f(meanmag)>=g(meanmag), interpolate=True, color='C3',alpha=0.3)


    # err = 10**(lM200_p1 - linear(meanmag_p1,m_mice,n_mice))*np.log(10.)*elM200_p1

    # axr.errorbar(meanmag_p1,10**(lM200_p1 - linear(meanmag_p1,m_mice,n_mice)),
                 # yerr=err,markersize=10,label='MICE - PAUS\,1',
                 # fmt='C7',marker='x')
    
    axr.plot(meanmag,np.ones(len(meanmag)),'k--')
    
    axr.grid()
    # axr.legend(loc=2,ncol=5,fontsize=11)
    axr.set_xlabel(r'$\langle M_r \rangle$')
    axr.set_ylabel(r'$M_{200}/M^{Viola}_{200}$')
    axr.set_ylim(0,2.5)
    
    fig.savefig('../final_plots/mass_mag.pdf',bbox_inches='tight')

    
def make_lum_mass_plot():
    
    meanmag_mice = np.array([-20.03361442, -20.57903032, -21.2230643 , -21.84807599,-22.48726666])
    lM200_mice  = np.array([11.32221929, 11.54406804, 11.92427929, 12.22530928, 12.67117284])
    
    meanmag  = np.arange(-24,-20.,0.2)
    
    MGirardi = (np.array([14.12,14.35,14.41,14.43,14.26,14.38,13.65,13.86,14.66,14.45]))+np.log10(0.7)
    eMGirardi = np.array([0.25,0.14,0.23,0.19,0.24,0.20,0.10,0.23,0.10,0.14])
    LGirardi = (np.array([11.85,12.14,12.10,12.05,11.75,11.79,11.41,11.81,12.20,12.10]))+np.log10(0.7**2)
    
    L_Mandel    = np.log10(np.array([0.075,0.19,0.47,1.1,2.1,3.2,5.0])*1.2e10)

    M_Mandel    = np.log10(np.array([3.16,2.2,4.1,14.9,54.,157.,578.])*1.e11*0.77)
    eM_MandelM  = (np.array([40.,7.9,4.1,5.0,14.,46.,180.])*1.e11*0.77)/(10**M_Mandel*np.log(10.))
    eM_Mandelm  = (np.array([3.13,2.1,3.0,4.6,14.,43.,114.])*1.e11*0.77)/(10**M_Mandel*np.log(10.))
    ML_Mandel   = np.array([246,69,51,79,151,285,674])*0.77*(1./0.7)
    eML_MandelM = np.array([3113,242,51,27,38,83,78])*0.77*(1./0.7)
    eML_Mandelm = np.array([243,66,37,24,38,78,203])*0.77*(1./0.7)

    M_Mandel2     = np.log10(np.array([0.31,4.0,4.4,7.1,26.,21.,100.])*1.e11*0.77)
    eM_MandelM2   = (np.array([6.44,2.7,3.4,2.8,15.,48.,208.])*1.e11*0.77)/(10**M_Mandel*np.log(10.))
    eM_Mandelm2   = (np.array([0.30,3.9,2.6,3.0,13.,21.,100.])*1.e11*0.77)/(10**M_Mandel*np.log(10.))
    ML_Mandel2   = np.array([26,129,57,41,74,40,122])*0.77*(1./0.7)
    eML_MandelM2 = np.array([507,88,44,16,43,90,253])*0.77*(1./0.7)
    eML_Mandelm2 = np.array([24,128,33,17,37,40,121])*0.77*(1./0.7)
    
    Msun = 4.65
    lM0  = np.log10(1.e14)
    L0   = (10**11.0)
    Lsun = (10**(-0.4*Msun))
    
    

    def linear_viola(x,m,ln):
        return ln+m*(x-11.5)+14.
        
    def linear(x,m,ln):
        return ln+m*(x-np.log10(L0))+14.
    
    
    from scipy.optimize import curve_fit

    # MICE MASS-TO-LIGHT RELATION
    
    meanmag      = -0.4*(meanmag-Msun) 
    meanmag_mice = -0.4*(meanmag_mice-Msun) 
    
    popt, pcov = curve_fit(linear, meanmag_mice, lM200_mice)
    m_mice,ln_mice = popt
    n_mice  = 10**ln_mice
    print('MICE sample fit')
    print(m_mice,n_mice)    


    dL0 = (10**11.5)/L0
    m_v = 1.16
    em_v = 0.13
    n_v = 0.95
    en_v = 0.14
    ln_v = np.log10(n_v)
    eln_v = en_v/(np.log(10.)*n_v)
    
    print('Viola et. al.')
    print(m_v,n_v)    
    print(em_v,en_v)    

        
    # SAMPLES AND LABELS
    
    samples =  ['pk_w1w3_'+pcat,
            'pk_Mm2_w1w3_'+pcat,'pk_MM2_w1w3_'+pcat,
            'pk_zm_w1w3_'+pcat,'pk_zM_w1w3_'+pcat,
            'pk_Lrm2_w1w3_'+pcat,'pk_LrM2_w1w3_'+pcat,
            'pk_blue_w1w3_'+pcat,'pk_red_w1w3_'+pcat]

    samples_gold =  ['pk_w1w3_'+best,
            'pk_Mm2_w1w3_'+best,'pk_MM2_w1w3_'+best,
            'pk_zm_w1w3_'+best,'pk_zM_w1w3_'+best,
            'pk_Lrm2_w1w3_'+best,'pk_LrM2_w1w3_'+best,
            'pk_blue_w1w3_'+best,'pk_red_w1w3_'+best]
    
    
    csamp = ['k',
            'gold','gold',
            'royalblue','royalblue',
            'C9','C9',
            'palevioletred','palevioletred']
    
    lsamp = ['all pairs',
            r'$M^\mathrm{pair}_r < -22.5$',
            r'$M^\mathrm{pair}_r \geq -22.5$',
            '$z < 0.4$',
            '$z \geq 0.4$',
            '$L_2/L_1 < 0.5$',
            '$L_2/L_1 \geq 0.5$',
            r'$blue\,\,pairs$',
            r'$red\,\,pairs$']

    mark = ['o'] + ['v','^']*4
    
    # EXTRACTING DATA
    
    lMfit = []
    lMfit_gold = []
    
    for j in range(len(samples)):
        lMfit += [plt_profile_fit_2h(samples[j],lsamp[j],plot=False,fytpe = ftype)]
        lMfit_gold += [plt_profile_fit_2h(samples_gold[j],lsamp[j],plot=False,fytpe = ftype)]


    # MAKE PLOTS

    fig, axall = plt.subplots(3,1, figsize=(10,10),sharex = True,gridspec_kw={'height_ratios': [4, 1,1]})
    fig.subplots_adjust(hspace=0,wspace=0)
    # fig, ax = plt.subplots(1,1, figsize=(10,7))

    ax  = axall[0]
    axrt = axall[1]
    axrg = axall[2]


    # ax.errorbar(0,0,markersize=5,
                     # yerr=np.array([np.diff(lMfit[j][0])]).T,
                     # fmt='w',marker=',',label=' ')

    ax.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt='k',marker='o',label='Total sample')
    ax.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt='k',marker='o',label='Gold sample', mfc='w',mew=3)
                     
    axrt.plot(0,0,'w,',label='Total sample')
    axrg.plot(0,0,'w,',label='Gold sample')
    axrt.legend(frameon=False,loc=4)
    axrg.legend(frameon=False,loc=4)
    # axr.errorbar(0,0,markersize=10,
                     # yerr=np.array([np.diff(lMfit[j][0])]).T,
                     # fmt='k',marker='o',label='Total sample')
    # axr.errorbar(0,0,markersize=10,
                     # yerr=np.array([np.diff(lMfit[j][0])]).T,
                     # fmt='k',marker='o',label='Gold sample', mfc='w',mew=3)

    
    
    MM      = []
    LM      = []
    eLM      = []
    MM_gold = []
    LM_gold = []
    eLM_gold = []
    
    for j in np.arange(len(csamp)):
        
        
        
        MM       += [-0.4*(lMfit[j][1]-Msun)]
        LM       += [lMfit[j][0][1]]
        eLM      += [np.diff(lMfit[j][0])]
        MM_gold  += [-0.4*(lMfit_gold[j][1]-Msun)]
        LM_gold  += [lMfit_gold[j][0][1]]
        eLM_gold += [np.diff(lMfit_gold[j][0])]
        

        # MAKE SCATTER PLOT
        if j < 3:
            ax.errorbar(MM[-1],lMfit[j][0][1],markersize=15,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt=csamp[j],marker=mark[j])
            ax.errorbar(MM_gold[-1],lMfit_gold[j][0][1],markersize=15,
                     yerr=np.array([np.diff(lMfit_gold[j][0])]).T,
                     fmt=csamp[j],marker=mark[j], mfc='w',mew=3)
                     
        
        # JUST FOR THE LABELS

        if j > 0:
            ax.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt=csamp[j],marker=mark[j],label=lsamp[j])
            # axr.errorbar(0,0,markersize=10,
                     # yerr=np.array([np.diff(lMfit[j][0])]).T,
                     # fmt=csamp[j],marker=mark[j],label=lsamp[j])


    MM = np.array(MM)
    MM_gold = np.array(MM_gold)
    LM_gold = np.array(LM_gold)
    eLM_gold = np.array(eLM_gold)
    LM = np.array(LM)
    eLM = np.array(eLM)
    
    
    # FIT MASS-TO-LIGHT RELATION

    # mlim = 100.
    # mask = (MM < mlim)
    
    mask = [False,True,True,False,False,False,False,False,False]
    popt, pcov = curve_fit(linear, MM[mask], LM[mask])
    popt1, pcov1 = curve_fit(linear, MM[mask], LM[mask] + [eLM[mask][0][0], -1.*eLM[mask][1][1]])
    popt2, pcov2 = curve_fit(linear, MM[mask], LM[mask] - [eLM[mask][0][1], -1.*eLM[mask][1][0]])
    
    m,ln = popt
    em,eln = np.array([popt1 - popt, popt - popt2]).T
    n  = 10**ln
    en = np.abs(ln*np.log(10))*eln
    
    print('Total sample fit')
    print(m,n)
    print(em,en)

    popt, pcov = curve_fit(linear, MM_gold[mask], LM_gold[mask])
    popt1, pcov1 = curve_fit(linear, MM_gold[mask], LM_gold[mask] + [eLM_gold[mask][0][0], -1.*eLM_gold[mask][1][1]])
    popt2, pcov2 = curve_fit(linear, MM_gold[mask], LM_gold[mask] - [eLM_gold[mask][0][1], -1.*eLM_gold[mask][1][0]])
    
    m_gold,ln_gold = popt
    em_gold,eln_gold = np.array([popt1 - popt, popt - popt2]).T
    n_gold  = 10**ln_gold
    en_gold = np.abs(ln_gold*np.log(10))*eln_gold

    print('Gold sample fit')
    print(m_gold,n_gold)
    print(em_gold,en_gold)
    
    mask = mask*2
    eLM_all = np.vstack((eLM,eLM_gold))
    popt, pcov = curve_fit(linear, np.append(MM,MM_gold)[mask], np.append(LM,LM_gold)[mask],sigma=np.average(eLM_all,axis=1)[mask])
    m_all,ln_all = popt
    em_all,eln_all = np.abs(np.sqrt(np.diag(pcov)))
    n_all  = 10**ln_all
    en_all = np.abs(ln_all*np.log(10))*eln_all
       
    print('All sample fit')
    print(m_all,n_all)    
    print(em_all,en_all)    


    popt, pcov = curve_fit(linear, L_Mandel, M_Mandel,sigma=eM_MandelM)
    m_mandel,ln_mandel = popt
    em_mandel,eln_mandel = np.abs(np.sqrt(np.diag(pcov)))
    n_mandel  = 10**ln_mandel
    en_mandel = np.abs(ln_mandel*np.log(10))*eln_mandel

    print('Mandel sample fit')
    print(m_mandel,n_mandel)
    print(em_mandel,en_mandel)

    for j in np.arange(len(csamp)):
        
        # MAKE MICE RATIO PLOT

        err = np.abs(10**(LM[j]-linear(MM[j],m,ln))*np.log(10.))*eLM[j]
        axrt.errorbar(MM[j],10**(-1.*LM[j]+linear(MM[j],m,ln)),
                        yerr=err.reshape(2,1),markersize=10,
                        fmt=csamp[j],marker=mark[j])
        
        err = np.abs(10**(LM_gold[j]-linear(MM_gold[j],m_gold,ln_gold))*np.log(10.))*eLM_gold[j]
        axrg.errorbar(MM_gold[j],10**(-1.*LM_gold[j]+linear(MM_gold[j],m_gold,ln_gold)),
                        yerr=err.reshape(2,1),markersize=10,
                        fmt=csamp[j],marker=mark[j],mfc='w',mew=3)
    
    axrt.plot([10.3,11.2],[1,1],'C7--')
    axrg.plot([10.3,11.2],[1,1],'C7--')


    ax.plot(meanmag,linear_viola(meanmag,m_v,ln_v),'C8',
            # label=r'Group/Clusters, $\alpha = $'+str(np.round(m_v,2))+'$\pm$'+str(np.round(em_v,2))+r', $\beta = $'+str(np.round(n_v/dL0**m_v,2))+'$\pm$'+str(np.round(en_v/dL0**m_v,2)))
            label=r'Group/Clusters, $\alpha = $'+str(np.round(m_v,2))+'$\pm$'+str(np.round(em_v,2))+r'$, \beta = $'+str(np.round(ln_v,2))+'$\pm$'+str(np.round(eln_v,2)))

    ax.plot(meanmag,linear(meanmag,m_mice,ln_mice),'C7',
            label=r'MICE, $\alpha = $'+str(np.round(m_mice,2))+r', $\beta = $'+str(np.round(ln_mice,2)))
            
    ax.fill_between(meanmag, 
                    linear_viola(meanmag,m_v+em_v,ln_v+eln_v), 
                    linear_viola(meanmag,m_v-em_v,ln_v-eln_v), 
                    interpolate=True, 
                    color='C8',alpha=0.3)
    
    ax.plot(meanmag,linear(meanmag,m,ln),'C4',
            label=r'Total sample, $\alpha = $'+str(np.round(m,1))+r'$^{+'+str(np.round(em[0],1))+r'}_{-'+str(np.round(em[1],1))+'}$'+r'$, \beta = $'+str(np.round(ln,2))+r'$^{+'+str(np.round(eln[0],2))+r'}_{-'+str(np.round(eln[1],2))+'}$')

            
    ax.fill_between(meanmag, 
                    linear(meanmag,m-max(em),ln-max(eln)), 
                    linear(meanmag,m+max(em),ln+max(eln)), 
                    interpolate=True, 
                    color='C4',alpha=0.3)

    ax.plot(meanmag,linear(meanmag,m_gold,ln_gold),'C1',
            label=r'Gold sample, $\alpha = $'+str(np.round(m_gold,1))+r'$^{+'+str(np.round(em_gold[0],1))+r'}_{-'+str(np.round(em_gold[1],1))+'}$'+r', $\beta = $'+str(np.round(ln_gold,2))+r'$^{+'+str(np.round(eln_gold[0],2))+r'}_{-'+str(np.round(eln_gold[1],2))+'}$')
            
    ax.fill_between(meanmag, 
                    linear(meanmag,m_gold+max(em_gold),ln_gold+max(eln_gold)), 
                    linear(meanmag,m_gold-max(em_gold),ln_gold-max(eln_gold)), 
                    interpolate=True, 
                    color='C1',alpha=0.3)

    # ax.plot(meanmag,linear(meanmag,m_all,ln_all),'C3',
            # label=r'All pair sub-samples, $\alpha = $'+str(np.round(m_all,2))+'$\pm$'+str(np.round(em_all,2)))
            
    # ax.fill_between(meanmag, 
                    # linear(meanmag,m_all+em_all,ln_all+eln_all), 
                    # linear(meanmag,m_all-em_all,ln_all-eln_all), 
                    # interpolate=True, 
                    # color='C3',alpha=0.3)

    # ax.plot(meanmag,lM200,'C7',label='MICE - True pairs')
    # ax.plot(0,0,'w,',label=' ')
    
    ax.legend(loc=4,ncol=3,fontsize=11)
    ax.set_ylabel(r'$\log (M_{200}/M_\odot h^{-1})$')
    axrt.set_ylabel(r'$M_{200}/M_{fit}$')
    axrg.set_ylabel(r'$M_{200}/M_{fit}$')
    ax.axis([10.3,11.2,10.4,13.4])
    axrt.axis([10.3,11.2,-1.5,2.7])
    axrg.axis([10.3,11.2,-1.5,2.7])
    
    # axr.plot(meanmag,np.ones(len(meanmag)),'k--')
    
    # axr.grid()
    # axr.legend(loc=2,ncol=5,fontsize=11)
    axrg.set_xlabel(r'$ \log \langle L_{TOT}/(h^{-2} L_\odot) \rangle$')
    # axr.set_ylabel(r'$M^{Viola}_{200}/M_{200}$')
    # axr.set_ylim(0,7.)
    fig.savefig('../final_plots/mass_lum.pdf',bbox_inches='tight')   
    
    # INNER PLOT
    # axin = axall[0].inset_axes([0.6,0.1,0.37,0.45]) 
 
    fig, axin = plt.subplots(1,1, figsize=(6,6))
 
    axin.plot(LGirardi,MGirardi,'o',color='olive')    
    axin.errorbar(LGirardi,MGirardi,markersize=4,
                     yerr=eMGirardi,fmt='none',
                     color='olive',marker='o')
                     
    axin.plot(MM_gold,LM_gold,'C1o')    
    axin.errorbar(MM_gold,LM_gold,markersize=4,
                     yerr=(eLM_gold).T,fmt='none',
                     color='C1',marker='o')
    axin.plot(MM,LM,'C4o')    
    axin.errorbar(MM,LM,markersize=4,
                     yerr=eLM.T,fmt='none',
                     color='C4',marker='o')
 
    
    meanmag = np.arange(9.2,14,0.5)
        
    axin.plot(meanmag,linear_viola(meanmag,m_v,ln_v),'C8--')
    axin.plot(meanmag[meanmag >= 10.5],linear_viola(meanmag,m_v,ln_v)[meanmag >= 10.5],'C8')
            
    axin.fill_between(meanmag, 
                    linear_viola(meanmag,m_v+em_v,ln_v+eln_v), 
                    linear_viola(meanmag,m_v-em_v,ln_v-eln_v), 
                    interpolate=True, 
                    color='C8',alpha=0.3)

    
    axin.plot(meanmag[(meanmag < 10.5)*(meanmag > 10.)],linear(meanmag,m,ln)[(meanmag < 10.5)*(meanmag > 10.)],'C4')
    axin.plot(meanmag,linear(meanmag,m,ln),'C4--')
            
    axin.fill_between(meanmag, 
                    linear(meanmag,m+em[1],ln+eln[1]), 
                    linear(meanmag,m-em[0],ln-eln[0]), 
                    interpolate=True, 
                    color='C4',alpha=0.3)

    axin.plot(meanmag,linear(meanmag,m_gold,ln_gold),'C1',
            label=r'Gold sample, $\alpha = $'+str(np.round(m_gold,2))+'$\pm$'+str(np.round(em_gold,2)))
            
    axin.fill_between(meanmag, 
                    linear(meanmag,m_gold+em_gold[1],ln_gold+eln_gold[1]), 
                    linear(meanmag,m_gold-em_gold[0],ln_gold-eln_gold[0]), 
                    interpolate=True, 
                    color='C1',alpha=0.3)

    # axin.plot(meanmag[(meanmag <= 10.5)*(meanmag >= 10.)],linear(meanmag,m_all,ln_all)[(meanmag <= 10.5)*(meanmag >= 10.)],'C3')
    # axin.plot(meanmag,linear(meanmag,m_all,ln_all),'C3--')
            # label=r'All pair sub-samples, $\alpha = $'+str(np.round(m_all,2))+'$\pm$'+str(np.round(em_all,2)))
            
    # axin.fill_between(meanmag, 
                    # linear(meanmag,m_all+em_all,ln_all+eln_all), 
                    # linear(meanmag,m_all-em_all,ln_all-eln_all), 
                    # interpolate=True, 
                    # color='C3',alpha=0.3)   
                    
    axin.plot(np.log10(0.991e10),np.log10(166.*0.991e10*0.6),'ko')
    axin.errorbar(np.log10(0.991e10),np.log10(166.*0.991e10*0.6),markersize=4,
                     yerr=32/(np.log(10)*166.),fmt='none',
                     color='k',marker='o')
                    
    axin.axis([9.8,12.5,11,15])
    
    axin.axvspan(9.8,10.1,color='C7',alpha=0.2)
    axin.axvspan(10.3,11.05,color='C3',alpha=0.2)
    axin.axvspan(10.45,12.5,color='C8',alpha=0.2)
    
    axin.text(9.85,13.7,'Galaxies',rotation='vertical')
    axin.text(10.31,13.7,'Pairs',rotation='vertical')
    axin.text(10.55,13.7,'Groups/Clusters',rotation='vertical')

    axin.set_xlabel(r'$ \log \langle L_{TOT}/(h^{-2} L_\odot) \rangle$')
    axin.set_ylabel(r'$\log (M_{200}/M_\odot h^{-1})$')
    
    fig.savefig('../final_plots/mass_lum_ext.pdf',bbox_inches='tight')  
    
    # MRATIO PLOT 
    
    fig, axin = plt.subplots(1,1, figsize=(6,6))
 
    # axin.plot(LGirardi,10**(MGirardi-LGirardi),'o',color='olive')    
    # axin.errorbar(LGirardi,10**(MGirardi-LGirardi),markersize=4,
                     # yerr=eMGirardi*np.log(10.)*10**(MGirardi-LGirardi),
                     # fmt='none',
                     # color='olive',marker='o')
                     
    axin.plot(MM_gold[:3],10**(LM_gold-MM_gold)[:3],'C1o')    
    axin.errorbar(MM_gold[:3],10**(LM_gold-MM_gold)[:3],markersize=4,
                     yerr=np.array([(eLM_gold.T[0]*10**(LM_gold-MM_gold)*np.log(10.))[:3],(eLM_gold.T[1]*10**(LM_gold-MM_gold)*np.log(10.))[:3]]),
                     fmt='none',
                     color='C1',marker='o')
                     
    axin.plot(MM[:3],10**(LM-MM)[:3],'C4o')    
    axin.errorbar(MM[:3],10**(LM-MM)[:3],markersize=4,
                     yerr=np.array([(eLM.T[0]*10**(LM-MM)*np.log(10.))[:3],(eLM.T[1]*10**(LM-MM)*np.log(10.))[:3]]),
                     fmt='none',
                     color='C4',marker='o')
 
    
    meanmag = np.arange(9.,14,0.1)
        
    axin.plot(meanmag,10**(linear_viola(meanmag,m_v,ln_v)-meanmag),'C8--')
    axin.plot(meanmag[meanmag >= 10.5],10**(linear_viola(meanmag,m_v,ln_v)-meanmag)[meanmag >= 10.5],'C8')
            
    axin.fill_between(meanmag, 
                    10**(linear_viola(meanmag,m_v+em_v,ln_v+eln_v)-meanmag), 
                    10**(linear_viola(meanmag,m_v-em_v,ln_v-eln_v)-meanmag), 
                    interpolate=True, 
                    color='C8',alpha=0.3)

    

    # axin.plot(meanmag[(meanmag <= 10.5)*(meanmag >= 10.)],10**(linear(meanmag,m_all,ln_all)-meanmag)[(meanmag <= 10.5)*(meanmag >= 10.)],'C3')
    # axin.plot(meanmag,10**(linear(meanmag,m_all,ln_all)-meanmag),'C3--')
            
    # axin.fill_between(meanmag, 
                    # 10**(linear(meanmag,m_all+em_all,ln_all+eln_all)-meanmag), 
                    # 10**(linear(meanmag,m_all-em_all,ln_all-eln_all)-meanmag), 
                    # interpolate=True, 
                    # color='C3',alpha=0.3)   

    axin.plot(meanmag[(meanmag <= 11.1)*(meanmag >= 10.4)],10**(linear(meanmag,m_gold,ln_gold)-meanmag)[(meanmag <= 11.1)*(meanmag >= 10.4)],'C1')
    axin.plot(meanmag,10**(linear(meanmag,m_gold,ln_gold)-meanmag),'C1--')
            
    axin.fill_between(meanmag, 
                    10**(linear(meanmag,m_gold+em_gold[1],ln_gold+eln_gold[1])-meanmag), 
                    10**(linear(meanmag,m_gold-em_gold[0],ln_gold-eln_gold[0])-meanmag), 
                    interpolate=True, 
                    color='C1',alpha=0.3)   

    axin.plot(meanmag[(meanmag <= 11.1)*(meanmag >= 10.4)],10**(linear(meanmag,m,ln)-meanmag)[(meanmag <= 11.1)*(meanmag >= 10.4)],'C4')
    axin.plot(meanmag,10**(linear(meanmag,m,ln)-meanmag),'C4--')
            
    axin.fill_between(meanmag, 
                    10**(linear(meanmag,m+em[1],ln+eln[1])-meanmag), 
                    10**(linear(meanmag,m-em[0],ln-eln[0])-meanmag), 
                    interpolate=True, 
                    color='C4',alpha=0.3)   
                    
    # axin.plot(np.log10(0.991e10),166.*0.6,'ko')
    # axin.errorbar(np.log10(0.991e10),166.*0.6,markersize=4,
                     # yerr=32*0.6,fmt='none',
                     # color='k',marker='o')

    axin.plot(L_Mandel,ML_Mandel,'ko')
    axin.errorbar(L_Mandel,ML_Mandel,markersize=4,
                     yerr=np.array([eML_MandelM,eML_Mandelm]),fmt='none',
                     color='k',marker='o')

    axin.plot(L_Mandel,ML_Mandel2,'kx')
    axin.errorbar(L_Mandel,ML_Mandel2,markersize=4,
                     yerr=np.array([eML_MandelM2,eML_Mandelm2]),fmt='none',
                     color='k',marker='o')
                    
    axin.axis([9.2,12.5,0,800])
    
    # axin.axhspan((166-32)*0.6,(166+32)*0.6,color='C7',alpha=0.2)
    # axin.axvspan(9.2,10.1,color='C7',alpha=0.2)
    # axin.axvspan(10.3,11.05,color='C3',alpha=0.2)
    # axin.axvspan(10.45,12.5,color='C8',alpha=0.2)
    
    # axin.plot([9.2,10.1],[166*0.6]*2,'k')
    # axin.plot([9.2,12.5],[166*0.6]*2,'k--')
    axin.text(9.3,500,'Galaxies',rotation='vertical')
    axin.text(10.31,500,'Pairs',rotation='vertical')
    axin.text(10.8,500,'Groups/Clusters',rotation='vertical')

    axin.set_xlabel(r'$ \log \langle L_{TOT}/(h^{-2} L_\odot) \rangle$')
    axin.set_ylabel(r'$M_{200}/L_{TOT} h$')
    
    
    fig.savefig('../final_plots/mlratio.pdf',bbox_inches='tight')   

         
def compare_MICE_mags():

    # LOAD MICE CAT
    
    MICE = fits.open('../MICE_cat.fits')[1].data
    
    z_h  = MICE.z_cgal
    mi_h = MICE.sdss_i_true - 0.8 * (np.arctan(1.5 * z_h) - 0.1489)
    mr_h = MICE.sdss_r_true - 0.8 * (np.arctan(1.5 * z_h) - 0.1489)
    
    Dl_h = np.array(cosmo.luminosity_distance(z_h).value)*1.e6 
    
    Mr_h=mr_h+5.0-5.0*np.log10(Dl_h)
    Mi_h=mi_h+5.0-5.0*np.log10(Dl_h)
    
    maskh  = (z_h > 0.2)*(z_h < 0.6)*(Mr_h < -19.5)
    maskh2 = (z_h > 0.2)*(z_h < 0.6)*(Mr_h < -19.5)*(MICE.nsats <3)*(MICE.nsats >1)

    # LOAD MICE PAIRS 
    
    M1 = np.loadtxt('../MICE-pairs/True-pairs-0.2-0.6-50-paralentes-R1.dat').T
    M2 = np.loadtxt('../MICE-pairs/True-pairs-0.2-0.6-50-paralentes-R2.dat').T
    M3 = np.loadtxt('../MICE-pairs/True-pairs-0.2-0.6-50-paralentes-R3.dat').T
    M4 = np.loadtxt('../MICE-pairs/True-pairs-0.2-0.6-50-paralentes-R4.dat').T
    
    M = np.vstack((M1.T,M2.T,M3.T,M4.T)).T
    
    M1_mice = np.zeros(len(M[1]))
    M2_mice = np.zeros(len(M[1]))
    
    M1_mice[M[8]<M[17]] = M[8][M[8]<M[17]]
    M1_mice[M[8]>M[17]] = M[17][M[8]>M[17]]
    M2_mice[M[8]<M[17]] = M[17][M[8]<M[17]]
    M2_mice[M[8]>M[17]] = M[8][M[8]>M[17]]
    
    Mtot_mice   = -2.5*np.log10(10**(-0.4*M[8])+10**(-0.4*M[17]))
    Mtot_i_mice = -2.5*np.log10(10**(-0.4*M[10])+10**(-0.4*M[19]))
    z_mice      = M[3]
    Lratio_mice = 10.**(-0.4*(M2_mice-M1_mice))
    color_mice  = Mtot_mice - Mtot_i_mice
    
    # LOAD PAUS PAIRS 
    
    L1 = np.loadtxt('../pareconk/Pares-PAUS_W1-'+pcat).T                                                            
    # L2 = np.loadtxt('../pareconk/Pares-PAUS_W2-'+pcat).T                                    
    L3 = np.loadtxt('../pareconk/Pares-PAUS_W3-'+pcat).T
    
    L = np.vstack((L1.T,L3.T)).T
    
    M1 = L[8]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    M2 = L[-1]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    Mtot = -2.5*np.log10(10**(-0.4*M1)+10**(-0.4*M2))
    M1i = L[7]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    M2i = L[-2]-5.*np.log10(np.array(cosmo.luminosity_distance(L[3]))*1.e6)+5
    Mtoti = -2.5*np.log10(10**(-0.4*M1i)+10**(-0.4*M2i))
    
    color = Mtot-Mtoti
    
    plt.figure()
    plt.plot(Mr_h[maskh],(mr_h-mi_h)[maskh],'C7.',label='MICE halos')
    plt.plot(Mr_h[maskh2],(mr_h-mi_h)[maskh2],'C0.',label='MICE halos - Nsats =  2',alpha=0.5)
    plt.plot(Mtot_mice,color_mice,'C2x',label='MICE pairs',alpha=0.3)
    plt.plot(Mtot,color,'C3x',label='PAUS pairs',alpha=0.5)
    plt.xlabel('$M_r$')
    plt.ylabel('$m_r - m_i$')
    plt.legend()
    plt.savefig('../color_mag.png')
    
    plt.figure()
    plt.hist((mr_h-mi_h)[maskh],np.linspace(0,1.2,50),edgecolor='C7',label='MICE halos',histtype='step',lw=3,density=True)
    plt.hist((mr_h-mi_h)[maskh2],np.linspace(0,1.2,50),edgecolor='C0',label='MICE halos - Nsats =  2',histtype='step',lw=3,density=True)
    plt.hist(color_mice,np.linspace(0,1.2,50),edgecolor='C2',label='MICE pairs',histtype='step',lw=3,density=True)
    plt.hist(color,np.linspace(0,1.2,50),edgecolor='C3',label='PAUS pairs',histtype='step',lw=3,density=True)
    plt.ylabel('$n$')
    plt.xlabel('$m_r - m_i$')
    plt.legend()
    plt.savefig('../color_hist.png')
    
def make_lum_mass_plot_wiso():
    
    meanmag_mice = np.array([-20.03361442, -20.57903032, -21.2230643 , -21.84807599,-22.48726666])
    lM200_mice  = np.array([11.32221929, 11.54406804, 11.92427929, 12.22530928, 12.67117284])
    
    meanmag  = np.arange(-24,-20.,0.2)
    
    MGirardi = (np.array([14.12,14.35,14.41,14.43,14.26,14.38,13.65,13.86,14.66,14.45]))+np.log10(0.7)
    eMGirardi = np.array([0.25,0.14,0.23,0.19,0.24,0.20,0.10,0.23,0.10,0.14])
    LGirardi = (np.array([11.85,12.14,12.10,12.05,11.75,11.79,11.41,11.81,12.20,12.10]))+np.log10(0.7**2)
    
    L_Mandel    = np.log10(np.array([0.075,0.19,0.47,1.1,2.1,3.2,5.0])*1.2e10)

    M_Mandel    = np.log10(np.array([3.16,2.2,4.1,14.9,54.,157.,578.])*1.e11*0.77)
    eM_MandelM  = (np.array([40.,7.9,4.1,5.0,14.,46.,180.])*1.e11*0.77)/(10**M_Mandel*np.log(10.))
    eM_Mandelm  = (np.array([3.13,2.1,3.0,4.6,14.,43.,114.])*1.e11*0.77)/(10**M_Mandel*np.log(10.))
    ML_Mandel   = np.array([246,69,51,79,151,285,674])*0.77*(1./0.7)
    eML_MandelM = np.array([3113,242,51,27,38,83,78])*0.77*(1./0.7)
    eML_Mandelm = np.array([243,66,37,24,38,78,203])*0.77*(1./0.7)

    M_Mandel2     = np.log10(np.array([0.31,4.0,4.4,7.1,26.,21.,100.])*1.e11*0.77)
    eM_MandelM2   = (np.array([6.44,2.7,3.4,2.8,15.,48.,208.])*1.e11*0.77)/(10**M_Mandel*np.log(10.))
    eM_Mandelm2   = (np.array([0.30,3.9,2.6,3.0,13.,21.,100.])*1.e11*0.77)/(10**M_Mandel*np.log(10.))
    ML_Mandel2   = np.array([26,129,57,41,74,40,122])*0.77*(1./0.7)
    eML_MandelM2 = np.array([507,88,44,16,43,90,253])*0.77*(1./0.7)
    eML_Mandelm2 = np.array([24,128,33,17,37,40,121])*0.77*(1./0.7)
    
    Msun = 4.65
    lM0  = np.log10(1.e14)
    L0   = (10**11.0)
    Lsun = (10**(-0.4*Msun))
    
    

    def linear_viola(x,m,ln):
        return ln+m*(x-11.5)+14.
        
    def linear(x,m,ln):
        return ln+m*(x-np.log10(L0))+14.
    
    
    from scipy.optimize import curve_fit

    # MICE MASS-TO-LIGHT RELATION
    
    meanmag      = -0.4*(meanmag-Msun) 
    meanmag_mice = -0.4*(meanmag_mice-Msun) 
    
    popt, pcov = curve_fit(linear, meanmag_mice, lM200_mice)
    m_mice,ln_mice = popt
    n_mice  = 10**ln_mice
    print('MICE sample fit')
    print(m_mice,n_mice)    


    dL0 = (10**11.5)/L0
    m_v = 1.16
    em_v = 0.13
    n_v = 0.95
    en_v = 0.14
    ln_v = np.log10(n_v)
    eln_v = en_v/(np.log(10.)*n_v)
    
    print('Viola et. al.')
    print(m_v,n_v)    
    print(em_v,en_v)    

        
    # SAMPLES AND LABELS
    
    samples =  ['pk_w1w3_'+pcat,
            'pk_Mm2_w1w3_'+pcat,'pk_MM2_w1w3_'+pcat,
            'pk_zm_w1w3_'+pcat,'pk_zM_w1w3_'+pcat,
            'pk_Lrm2_w1w3_'+pcat,'pk_LrM2_w1w3_'+pcat,
            'pk_blue_w1w3_'+pcat,'pk_red_w1w3_'+pcat,
            'pk_w1w3_'+pcat+'-7radios','pk_w1w3_'+pcat+'-10radios']

    samples_gold =  ['pk_w1w3_'+best,
            'pk_Mm2_w1w3_'+best,'pk_MM2_w1w3_'+best,
            'pk_zm_w1w3_'+best,'pk_zM_w1w3_'+best,
            'pk_Lrm2_w1w3_'+best,'pk_LrM2_w1w3_'+best,
            'pk_blue_w1w3_'+best,'pk_red_w1w3_'+best,
            'pk_w1w3_'+best+'-7radios','pk_w1w3_'+best+'-10radios']
    
    
    csamp = ['k',
            'gold','gold',
            'royalblue','royalblue',
            'C9','C9',
            'palevioletred','palevioletred',
            'k','k']
    
    lsamp = ['all pairs',
            r'$M^\mathrm{pair}_r < -22.5$',
            r'$M^\mathrm{pair}_r \geq -22.5$',
            '$z < 0.4$',
            '$z \geq 0.4$',
            '$L_2/L_1 < 0.5$',
            '$L_2/L_1 \geq 0.5$',
            r'$blue\,\,pairs$',
            r'$red\,\,pairs$',
            'iso 7',
            'iso 10']

    mark = ['o'] + ['v','^']*4 + ['s','D']
    
    # EXTRACTING DATA
    
    lMfit = []
    lMfit_gold = []
    
    for j in range(len(samples)):
        lMfit += [plt_profile_fit_2h(samples[j],lsamp[j],plot=False,fytpe = ftype)]
        lMfit_gold += [plt_profile_fit_2h(samples_gold[j],lsamp[j],plot=False,fytpe = ftype)]


    # MAKE PLOTS

    # fig, axall = plt.subplots(2,1, figsize=(10,10),sharex = True,gridspec_kw={'height_ratios': [3, 1]})
    # fig.subplots_adjust(hspace=0,wspace=0)
    fig, ax = plt.subplots(1,1, figsize=(10,7))

    # ax  = axall[0]
    # axr = axall[1]


    # ax.errorbar(0,0,markersize=5,
                     # yerr=np.array([np.diff(lMfit[j][0])]).T,
                     # fmt='w',marker=',',label=' ')

    ax.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt='k',marker='o',label='Total sample')
    ax.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt='k',marker='o',label='Gold sample', mfc='w',mew=3)
    # axr.errorbar(0,0,markersize=10,
                     # yerr=np.array([np.diff(lMfit[j][0])]).T,
                     # fmt='k',marker='o',label='Total sample')
    # axr.errorbar(0,0,markersize=10,
                     # yerr=np.array([np.diff(lMfit[j][0])]).T,
                     # fmt='k',marker='o',label='Gold sample', mfc='w',mew=3)

    
    
    MM      = []
    LM      = []
    eLM      = []
    MM_gold = []
    LM_gold = []
    eLM_gold = []
    
    for j in np.arange(len(csamp)):
        
        
        
        MM       += [-0.4*(lMfit[j][1]-Msun)]
        LM       += [lMfit[j][0][1]]
        eLM      += [np.diff(lMfit[j][0])]
        MM_gold  += [-0.4*(lMfit_gold[j][1]-Msun)]
        LM_gold  += [lMfit_gold[j][0][1]]
        eLM_gold += [np.diff(lMfit_gold[j][0])]
        

        # MAKE SCATTER PLOT

        ax.errorbar(MM[-1],lMfit[j][0][1],markersize=15,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt=csamp[j],marker=mark[j])
        ax.errorbar(MM_gold[-1],lMfit_gold[j][0][1],markersize=15,
                     yerr=np.array([np.diff(lMfit_gold[j][0])]).T,
                     fmt=csamp[j],marker=mark[j], mfc='w',mew=3)
                     
        '''
        # MAKE MICE RATIO PLOT

        err = np.abs(10**(lMfit[j][0][1]-linear_viola(MM[-1],m_v,ln_v))*np.log(10.))*np.array([np.diff(lMfit[j][0])]).T
        axr.errorbar(MM[-1],10**(-1.*lMfit[j][0][1]+linear_viola(MM[-1],m_v,ln_v)),
                        yerr=err,markersize=10,
                        fmt=csamp[j],marker=mark[j])
        
        err = np.abs(10**(lMfit_gold[j][0][1]-linear_viola(MM[-1],m_v,ln_v))*np.log(10.))*np.array([np.diff(lMfit_gold[j][0])]).T
        axr.errorbar(MM_gold[-1],10**(-1.*lMfit_gold[j][0][1]+linear_viola(MM[-1],m_v,ln_v)),
                        yerr=err,markersize=10,
                        fmt=csamp[j],marker=mark[j],mfc='w',mew=3)

        '''
        # JUST FOR THE LABELS

        if j > 0:
            ax.errorbar(0,0,markersize=10,
                     yerr=np.array([np.diff(lMfit[j][0])]).T,
                     fmt=csamp[j],marker=mark[j],label=lsamp[j])
            # axr.errorbar(0,0,markersize=10,
                     # yerr=np.array([np.diff(lMfit[j][0])]).T,
                     # fmt=csamp[j],marker=mark[j],label=lsamp[j])


    MM = np.array(MM)
    MM_gold = np.array(MM_gold)
    LM_gold = np.array(LM_gold)
    eLM_gold = np.array(eLM_gold)
    LM = np.array(LM)
    eLM = np.array(eLM)
    
    
    # FIT MASS-TO-LIGHT RELATION

    # mlim = 100.
    # mask = (MM < mlim)
    
    mask = [False,True,True,False,False,False,False,False,False,False,False]
    popt, pcov = curve_fit(linear, MM[mask], LM[mask])
    popt1, pcov1 = curve_fit(linear, MM[mask], LM[mask] + [eLM[mask][0][0], -1.*eLM[mask][1][1]])
    popt2, pcov2 = curve_fit(linear, MM[mask], LM[mask] - [eLM[mask][0][1], -1.*eLM[mask][1][0]])
    
    m,ln = popt
    em,eln = np.array([popt1 - popt, popt - popt2]).T
    n  = 10**ln
    en = np.abs(ln*np.log(10))*eln
    
    print('Total sample fit')
    print(m,n)
    print(em,en)

    popt, pcov = curve_fit(linear, MM_gold[mask], LM_gold[mask])
    popt1, pcov1 = curve_fit(linear, MM_gold[mask], LM_gold[mask] + [eLM_gold[mask][0][0], -1.*eLM_gold[mask][1][1]])
    popt2, pcov2 = curve_fit(linear, MM_gold[mask], LM_gold[mask] - [eLM_gold[mask][0][1], -1.*eLM_gold[mask][1][0]])
    
    m_gold,ln_gold = popt
    em_gold,eln_gold = np.array([popt1 - popt, popt - popt2]).T
    n_gold  = 10**ln_gold
    en_gold = np.abs(ln_gold*np.log(10))*eln_gold

    print('Gold sample fit')
    print(m_gold,n_gold)
    print(em_gold,en_gold)
    
    mask = mask*2
    eLM_all = np.vstack((eLM,eLM_gold))
    popt, pcov = curve_fit(linear, np.append(MM,MM_gold)[mask], np.append(LM,LM_gold)[mask],sigma=np.average(eLM_all,axis=1)[mask])
    m_all,ln_all = popt
    em_all,eln_all = np.abs(np.sqrt(np.diag(pcov)))
    n_all  = 10**ln_all
    en_all = np.abs(ln_all*np.log(10))*eln_all
       
    print('All sample fit')
    print(m_all,n_all)    
    print(em_all,en_all)    


    popt, pcov = curve_fit(linear, L_Mandel, M_Mandel,sigma=eM_MandelM)
    m_mandel,ln_mandel = popt
    em_mandel,eln_mandel = np.abs(np.sqrt(np.diag(pcov)))
    n_mandel  = 10**ln_mandel
    en_mandel = np.abs(ln_mandel*np.log(10))*eln_mandel

    print('Mandel sample fit')
    print(m_mandel,n_mandel)
    print(em_mandel,en_mandel)


    ax.plot(meanmag,linear_viola(meanmag,m_v,ln_v),'C8',
            # label=r'Group/Clusters, $\alpha = $'+str(np.round(m_v,2))+'$\pm$'+str(np.round(em_v,2))+r', $\beta = $'+str(np.round(n_v/dL0**m_v,2))+'$\pm$'+str(np.round(en_v/dL0**m_v,2)))
            label=r'Group/Clusters, $\alpha = $'+str(np.round(m_v,2))+'$\pm$'+str(np.round(em_v,2))+r'$, \beta = $'+str(np.round(ln_v,2))+'$\pm$'+str(np.round(eln_v,2)))

    ax.plot(meanmag,linear(meanmag,m_mice,ln_mice),'C7',
            label=r'MICE, $\alpha = $'+str(np.round(m_mice,2))+r', $\beta = $'+str(np.round(ln_mice,2)))
            
    ax.fill_between(meanmag, 
                    linear_viola(meanmag,m_v+em_v,ln_v+eln_v), 
                    linear_viola(meanmag,m_v-em_v,ln_v-eln_v), 
                    interpolate=True, 
                    color='C8',alpha=0.3)
    
    ax.plot(meanmag,linear(meanmag,m,ln),'C4',
            label=r'Total sample, $\alpha = $'+str(np.round(m,2))+r'$^{+'+str(np.round(em[0],2))+r'}_{-'+str(np.round(em[1],2))+'}$'+r'$, \beta = $'+str(np.round(ln,2))+r'$^{+'+str(np.round(eln[0],2))+r'}_{-'+str(np.round(eln[1],2))+'}$')

            
    ax.fill_between(meanmag, 
                    linear(meanmag,m-max(em),ln-max(eln)), 
                    linear(meanmag,m+max(em),ln+max(eln)), 
                    interpolate=True, 
                    color='C4',alpha=0.3)

    ax.plot(meanmag,linear(meanmag,m_gold,ln_gold),'C1',
            label=r'Gold sample, $\alpha = $'+str(np.round(m_gold,2))+r'$^{+'+str(np.round(em_gold[0],2))+r'}_{-'+str(np.round(em_gold[1],2))+'}$'+r', $\beta = $'+str(np.round(ln_gold,2))+r'$^{+'+str(np.round(eln_gold[0],2))+r'}_{-'+str(np.round(eln_gold[1],2))+'}$')
            
    ax.fill_between(meanmag, 
                    linear(meanmag,m_gold+max(em_gold),ln_gold+max(eln_gold)), 
                    linear(meanmag,m_gold-max(em_gold),ln_gold-max(eln_gold)), 
                    interpolate=True, 
                    color='C1',alpha=0.3)

    # ax.plot(meanmag,linear(meanmag,m_all,ln_all),'C3',
            # label=r'All pair sub-samples, $\alpha = $'+str(np.round(m_all,2))+'$\pm$'+str(np.round(em_all,2)))
            
    # ax.fill_between(meanmag, 
                    # linear(meanmag,m_all+em_all,ln_all+eln_all), 
                    # linear(meanmag,m_all-em_all,ln_all-eln_all), 
                    # interpolate=True, 
                    # color='C3',alpha=0.3)

    # ax.plot(meanmag,lM200,'C7',label='MICE - True pairs')
    # ax.plot(0,0,'w,',label=' ')
    
    ax.legend(loc=4,ncol=3,fontsize=11)
    ax.set_ylabel(r'$\log (M_{200}/M_\odot h^{-1})$')
    ax.axis([10.3,11.2,10.4,13.4])
    
    # axr.plot(meanmag,np.ones(len(meanmag)),'k--')
    
    # axr.grid()
    # axr.legend(loc=2,ncol=5,fontsize=11)
    ax.set_xlabel(r'$ \log \langle L_{TOT}/(h^{-2} L_\odot) \rangle$')
    # axr.set_ylabel(r'$M^{Viola}_{200}/M_{200}$')
    # axr.set_ylim(0,7.)
    fig.savefig('../final_plots/mass_lum_iso.pdf',bbox_inches='tight')   
    
    # INNER PLOT
    # axin = axall[0].inset_axes([0.6,0.1,0.37,0.45]) 
 
    fig, axin = plt.subplots(1,1, figsize=(6,6))
 
    axin.plot(LGirardi,MGirardi,'o',color='olive')    
    axin.errorbar(LGirardi,MGirardi,markersize=4,
                     yerr=eMGirardi,fmt='none',
                     color='olive',marker='o')
                     
    axin.plot(MM_gold,LM_gold,'C1o')    
    axin.errorbar(MM_gold,LM_gold,markersize=4,
                     yerr=(eLM_gold).T,fmt='none',
                     color='C1',marker='o')
    axin.plot(MM,LM,'C4o')    
    axin.errorbar(MM,LM,markersize=4,
                     yerr=eLM.T,fmt='none',
                     color='C4',marker='o')
 
    
    meanmag = np.arange(9.2,14,0.5)
        
    axin.plot(meanmag,linear_viola(meanmag,m_v,ln_v),'C8--')
    axin.plot(meanmag[meanmag >= 10.5],linear_viola(meanmag,m_v,ln_v)[meanmag >= 10.5],'C8')
            
    axin.fill_between(meanmag, 
                    linear_viola(meanmag,m_v+em_v,ln_v+eln_v), 
                    linear_viola(meanmag,m_v-em_v,ln_v-eln_v), 
                    interpolate=True, 
                    color='C8',alpha=0.3)

    
    axin.plot(meanmag[(meanmag < 10.5)*(meanmag > 10.)],linear(meanmag,m,ln)[(meanmag < 10.5)*(meanmag > 10.)],'C4')
    axin.plot(meanmag,linear(meanmag,m,ln),'C4--')
            
    axin.fill_between(meanmag, 
                    linear(meanmag,m+em[1],ln+eln[1]), 
                    linear(meanmag,m-em[0],ln-eln[0]), 
                    interpolate=True, 
                    color='C4',alpha=0.3)

    axin.plot(meanmag,linear(meanmag,m_gold,ln_gold),'C1',
            label=r'Gold sample, $\alpha = $'+str(np.round(m_gold,2))+'$\pm$'+str(np.round(em_gold,2)))
            
    axin.fill_between(meanmag, 
                    linear(meanmag,m_gold+em_gold[1],ln_gold+eln_gold[1]), 
                    linear(meanmag,m_gold-em_gold[0],ln_gold-eln_gold[0]), 
                    interpolate=True, 
                    color='C1',alpha=0.3)

    # axin.plot(meanmag[(meanmag <= 10.5)*(meanmag >= 10.)],linear(meanmag,m_all,ln_all)[(meanmag <= 10.5)*(meanmag >= 10.)],'C3')
    # axin.plot(meanmag,linear(meanmag,m_all,ln_all),'C3--')
            # label=r'All pair sub-samples, $\alpha = $'+str(np.round(m_all,2))+'$\pm$'+str(np.round(em_all,2)))
            
    # axin.fill_between(meanmag, 
                    # linear(meanmag,m_all+em_all,ln_all+eln_all), 
                    # linear(meanmag,m_all-em_all,ln_all-eln_all), 
                    # interpolate=True, 
                    # color='C3',alpha=0.3)   
                    
    axin.plot(np.log10(0.991e10),np.log10(166.*0.991e10*0.6),'ko')
    axin.errorbar(np.log10(0.991e10),np.log10(166.*0.991e10*0.6),markersize=4,
                     yerr=32/(np.log(10)*166.),fmt='none',
                     color='k',marker='o')
                    
    axin.axis([9.8,12.5,11,15])
    
    axin.axvspan(9.8,10.1,color='C7',alpha=0.2)
    axin.axvspan(10.3,11.05,color='C3',alpha=0.2)
    axin.axvspan(10.45,12.5,color='C8',alpha=0.2)
    
    axin.text(9.85,13.7,'Galaxies',rotation='vertical')
    axin.text(10.31,13.7,'Pairs',rotation='vertical')
    axin.text(10.55,13.7,'Groups/Clusters',rotation='vertical')

    axin.set_xlabel(r'$ \log \langle L_{TOT}/(h^{-2} L_\odot) \rangle$')
    axin.set_ylabel(r'$\log (M_{200}/M_\odot h^{-1})$')
    
    fig.savefig('../final_plots/mass_lum_ext_iso.pdf',bbox_inches='tight')  
    
    # MRATIO PLOT 
    
    fig, axin = plt.subplots(1,1, figsize=(6,6))
 
    # axin.plot(LGirardi,10**(MGirardi-LGirardi),'o',color='olive')    
    # axin.errorbar(LGirardi,10**(MGirardi-LGirardi),markersize=4,
                     # yerr=eMGirardi*np.log(10.)*10**(MGirardi-LGirardi),
                     # fmt='none',
                     # color='olive',marker='o')
                     
    axin.plot(MM_gold[:3],10**(LM_gold-MM_gold)[:3],'C1o')    
    axin.errorbar(MM_gold[:3],10**(LM_gold-MM_gold)[:3],markersize=4,
                     yerr=np.array([(eLM_gold.T[0]*10**(LM_gold-MM_gold)*np.log(10.))[:3],(eLM_gold.T[1]*10**(LM_gold-MM_gold)*np.log(10.))[:3]]),
                     fmt='none',
                     color='C1',marker='o')
                     
    axin.plot(MM[:3],10**(LM-MM)[:3],'C4o')    
    axin.errorbar(MM[:3],10**(LM-MM)[:3],markersize=4,
                     yerr=np.array([(eLM.T[0]*10**(LM-MM)*np.log(10.))[:3],(eLM.T[1]*10**(LM-MM)*np.log(10.))[:3]]),
                     fmt='none',
                     color='C4',marker='o')
 
    
    meanmag = np.arange(9.,14,0.1)
        
    axin.plot(meanmag,10**(linear_viola(meanmag,m_v,ln_v)-meanmag),'C8--')
    axin.plot(meanmag[meanmag >= 10.5],10**(linear_viola(meanmag,m_v,ln_v)-meanmag)[meanmag >= 10.5],'C8')
            
    axin.fill_between(meanmag, 
                    10**(linear_viola(meanmag,m_v+em_v,ln_v+eln_v)-meanmag), 
                    10**(linear_viola(meanmag,m_v-em_v,ln_v-eln_v)-meanmag), 
                    interpolate=True, 
                    color='C8',alpha=0.3)

    

    # axin.plot(meanmag[(meanmag <= 10.5)*(meanmag >= 10.)],10**(linear(meanmag,m_all,ln_all)-meanmag)[(meanmag <= 10.5)*(meanmag >= 10.)],'C3')
    # axin.plot(meanmag,10**(linear(meanmag,m_all,ln_all)-meanmag),'C3--')
            
    # axin.fill_between(meanmag, 
                    # 10**(linear(meanmag,m_all+em_all,ln_all+eln_all)-meanmag), 
                    # 10**(linear(meanmag,m_all-em_all,ln_all-eln_all)-meanmag), 
                    # interpolate=True, 
                    # color='C3',alpha=0.3)   

    axin.plot(meanmag[(meanmag <= 11.1)*(meanmag >= 10.4)],10**(linear(meanmag,m_gold,ln_gold)-meanmag)[(meanmag <= 11.1)*(meanmag >= 10.4)],'C1')
    axin.plot(meanmag,10**(linear(meanmag,m_gold,ln_gold)-meanmag),'C1--')
            
    axin.fill_between(meanmag, 
                    10**(linear(meanmag,m_gold+em_gold[1],ln_gold+eln_gold[1])-meanmag), 
                    10**(linear(meanmag,m_gold-em_gold[0],ln_gold-eln_gold[0])-meanmag), 
                    interpolate=True, 
                    color='C1',alpha=0.3)   

    axin.plot(meanmag[(meanmag <= 11.1)*(meanmag >= 10.4)],10**(linear(meanmag,m,ln)-meanmag)[(meanmag <= 11.1)*(meanmag >= 10.4)],'C4')
    axin.plot(meanmag,10**(linear(meanmag,m,ln)-meanmag),'C4--')
            
    axin.fill_between(meanmag, 
                    10**(linear(meanmag,m+em[1],ln+eln[1])-meanmag), 
                    10**(linear(meanmag,m-em[0],ln-eln[0])-meanmag), 
                    interpolate=True, 
                    color='C4',alpha=0.3)   
                    
    # axin.plot(np.log10(0.991e10),166.*0.6,'ko')
    # axin.errorbar(np.log10(0.991e10),166.*0.6,markersize=4,
                     # yerr=32*0.6,fmt='none',
                     # color='k',marker='o')

    axin.plot(L_Mandel,ML_Mandel,'ko')
    axin.errorbar(L_Mandel,ML_Mandel,markersize=4,
                     yerr=np.array([eML_MandelM,eML_Mandelm]),fmt='none',
                     color='k',marker='o')

    axin.plot(L_Mandel,ML_Mandel2,'kx')
    axin.errorbar(L_Mandel,ML_Mandel2,markersize=4,
                     yerr=np.array([eML_MandelM2,eML_Mandelm2]),fmt='none',
                     color='k',marker='o')
                    
    axin.axis([9.2,12.5,0,800])
    
    # axin.axhspan((166-32)*0.6,(166+32)*0.6,color='C7',alpha=0.2)
    # axin.axvspan(9.2,10.1,color='C7',alpha=0.2)
    # axin.axvspan(10.3,11.05,color='C3',alpha=0.2)
    # axin.axvspan(10.45,12.5,color='C8',alpha=0.2)
    
    # axin.plot([9.2,10.1],[166*0.6]*2,'k')
    # axin.plot([9.2,12.5],[166*0.6]*2,'k--')
    axin.text(9.3,500,'Galaxies',rotation='vertical')
    axin.text(10.31,500,'Pairs',rotation='vertical')
    axin.text(10.8,500,'Groups/Clusters',rotation='vertical')

    axin.set_xlabel(r'$ \log \langle L_{TOT}/(h^{-2} L_\odot) \rangle$')
    axin.set_ylabel(r'$M_{200}/L_{TOT} h$')
    
    
    fig.savefig('../final_plots/mlratio_iso.pdf',bbox_inches='tight')   
    

'''
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

