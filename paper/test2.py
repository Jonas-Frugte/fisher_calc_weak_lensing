#!/usr/bin/env python3
import numpy as np
import camb
import matplotlib.pyplot as plt

def get_camb_spectra(lmax=3000):
    """
    Returns:
      ell: array of multipoles [0..lmax]
      unlensed: shape (lmax+1,4) for [TT,EE,BB,TE]
      lensed:   shape (lmax+1,4) same
    in muK^2 units.
    """
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, tau=0.06)
    pars.InitPower.set_params(ns=0.965)
    pars.set_for_lmax(lmax, lens_potential_accuracy=0)
    pars.WantCls = True
    
    # Unlensed
    pars.DoLensing = False
    r1 = camb.get_results(pars)
    unl = r1.get_cmb_power_spectra(pars, CMB_unit='muK', raw_cl=True)['unlensed_scalar']
    
    # Lensed
    pars.DoLensing = True
    r2 = camb.get_results(pars)
    lens = r2.get_cmb_power_spectra(pars, CMB_unit='muK', raw_cl=True)['total']
    
    ell = np.arange(lmax+1)
    return ell, unl[:lmax+1,:], lens[:lmax+1,:]

def add_noise(lensed_cls, beam_arcmin=4.0, Delta_T=1.0, Delta_P=np.sqrt(2.0)):
    """
    Add white noise to the lensed spectra in muK^2.
    beam_arcmin: FWHM in arcmin
    Delta_T,P: noise in muK·arcmin
    """
    obs = np.copy(lensed_cls)
    lmax = obs.shape[0] - 1
    from math import log, exp
    sigma = beam_arcmin*(np.pi/10800.)
    dT = Delta_T*(np.pi/10800.)
    dP = Delta_P*(np.pi/10800.)
    
    for L in range(2, lmax+1):
        factor = (L*(L+1)*sigma*sigma)/(8.*log(2.))
        nTT = dT*dT*np.exp(factor)
        nPP = dP*dP*np.exp(factor)
        obs[L,0] += nTT  # TT
        obs[L,1] += nPP  # EE
        obs[L,2] += nPP  # BB
        # TE cross-noise = 0
    return obs

def noise_TETE_polar(ell, unl_cls, obs_cls, Lvals=None, lmin=2, lmax=3000, n_l=40, n_phi=32):
    """
    A *toy* polar integration for TETE lensing noise:
      N(L)^{-1} ~ ∫ d^2 l1/(2π)^2 [ f^2 / Den ] 
    with f ~ C_unl^TE, Den ~ [C_obs^TT*C_obs^EE - (C_obs^TE)^2], etc.

    Returns nTE[L] from Lvals in [lmin..lmax].
    """
    if Lvals is None:
        Lvals = range(lmin, lmax+1)
    
    def Cunl(l, idx):
        # idx: 0->TT,1->EE,2->BB,3->TE
        if l<0 or l>lmax: return 0.
        return unl_cls[int(l), idx]
    
    def Cobs(l, idx):
        if l<0 or l>lmax: return 0.
        return obs_cls[int(l), idx]
    
    nTE = np.zeros(lmax+1)
    ls = np.linspace(lmin, lmax, n_l)
    d_l = (lmax - lmin)/(n_l - 1)
    d_phi = 2.*np.pi/n_phi
    
    for iL in Lvals:
        L = float(iL)
        # We'll accumulate integrand in 'sum_val'
        sum_val = 0.
        for l1 in ls:
            # We'll do a real polar measure factor inside the loop:
            # measure ~ l1 d_l d_phi, plus factor 1/(2π)^2
            for ip in range(n_phi):
                phi = ip*d_phi
                l1x = l1*np.cos(phi)
                l1y = l1*np.sin(phi)
                l2x = L - l1x
                l2y = -l1y
                l2 = np.sqrt(l2x*l2x + l2y*l2y)
                if l2<lmin or l2>lmax:
                    continue
                
                # Dot products
                dotL1 = L*l1x
                dotL2 = L*l2x
                
                # Numerator
                f_num = Cunl(l1,3)*dotL1 + Cunl(l2,3)*dotL2
                # Denominator piece
                # We'll do a symmetrical version:
                Den1 = Cobs(l1,0)*Cobs(l2,1) - (Cobs(l1,3)*Cobs(l2,3))
                Den2 = Cobs(l2,0)*Cobs(l1,1) - (Cobs(l2,3)*Cobs(l1,3))
                den = Den1 + Den2
                
                if den<=0.0:
                    # if negative or zero, skip to avoid blow-ups
                    continue
                integrand = (f_num*f_num)/(den)
                
                # Now multiply by measure: area = (1/(2π)^2)*(l1 d_l d_phi)
                measure = (1./(2.*np.pi)**2)* (l1*d_l*d_phi)
                sum_val += integrand*measure
        
        # Finally, for TETE, there's typically a factor 1/(L^2) 
        # in front of the integral to get N(L)^{-1}.
        # So N(L) = [ (L^2)* sum_val ]^{-1}
        if sum_val>1e-30:
            nTE[iL] = 1./( (L**2)*sum_val )
        else:
            nTE[iL] = 0.
    
    return nTE

def main():
    lmax = 3000
    ell, unl, lens = get_camb_spectra(lmax)
    obs = add_noise(lens, beam_arcmin=4, Delta_T=1.0, Delta_P=np.sqrt(2.0))
    
    # We'll do TETE lensing noise from 2..3000
    n_te = noise_TETE_polar(ell, unl, obs, lmin=2, lmax=lmax, n_l=40, n_phi=32)
    
    # Plot
    plt.figure()
    Lplot = np.arange(lmax+1)
    mask = (Lplot>=2)
    plt.loglog(Lplot[mask], Lplot[mask]*(Lplot[mask]+1)*n_te[mask], label=r'$N_\ell^{(TE\,TE)}$ (toy)')
    plt.ylim(1e-9, 1e19)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell(\ell+1)\,N_\ell$')
    plt.title('TETE Noise × L(L+1) (Toy Integration)')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
