# follows quadratic estimator of https://arxiv.org/pdf/astro-ph/0111606

print('boba')

import numpy as np
from scipy.integrate import quad, dblquad
from cosm_setup import lensing_spectra
from scipy.interpolate import interp1d

# Instantiate your spectra object (which has a camb.results instance in spectra.results)
spectra = lensing_spectra(fiducial_k_nls=True)

# Get C_l for TT, EE, BB, TE in DIMENSIONLESS UNITS!! up to some lmax
lmax = 5000
cls = spectra.results.get_lensed_scalar_cls(lmax=lmax, raw_cl=True)

# cls[:,0] = TT, cls[:,1] = EE, cls[:,2] = BB, cls[:,3] = TE
ells = np.arange(lmax+1)
fTT = interp1d(ells, cls[:,0], kind='cubic', bounds_error=False, fill_value=0.0)
fEE = interp1d(ells, cls[:,1], kind='cubic', bounds_error=False, fill_value=0.0)
fBB = interp1d(ells, cls[:,2], kind='cubic', bounds_error=False, fill_value=0.0)
fTE = interp1d(ells, cls[:,3], kind='cubic', bounds_error=False, fill_value=0.0)

def cmbps(l, type1, type2):
    """
    Returns C_l in μK^2 for the given combination (type1, type2).
    type1, type2 can be 'T', 'E', or 'B'.
    """
    if (type1, type2) == ('T', 'T'):
        return fTT(l)
    elif (type1, type2) == ('E', 'E'):
        return fEE(l)
    elif (type1, type2) == ('B', 'B'):
        return fBB(l)
    elif (type1, type2) in [('T', 'E'), ('E', 'T')]:
        return fTE(l)
    else:
        # For non-standard combos like TB or EB (normally zero in standard ΛCDM),
        # just return zero unless you have them computed.
        return 0.0

arcmintorad = np.pi / 10800 # converts arcminutes to radians

def cmbps_noise(l, type1, type2, sigma, Delta_T, Delta_P):
    noise = 0

    # units to input:
    # sigma: arcmin
    # Delta_T, Delta_P: microKelvin arcmin

    # constant time!
    Tcmb = 2.728e6 # in micro kelvin
    sigma *= arcmintorad # in radians
    Delta_T *= arcmintorad
    Delta_P *= arcmintorad
    
    if type1 == 'T' and type2 == 'T':
        noise = (Delta_T / Tcmb)**2 * np.exp(l * (l + 1) * sigma**2 / (8 * np.log(2)))
    if (type1 == 'E' and type2 == 'E') or (type1 == 'B' and type2 == 'B'):
        noise = (Delta_P / Tcmb)**2 * np.exp(l * (l + 1) * sigma**2 / (8 * np.log(2)))

    return noise

def cmbps_obs(l, type1, type2, sigma, Delta_T, Delta_P):
    # observed powerspectrum (with noise)
    return cmbps(l, type1, type2) + cmbps_noise(l, type1, type2, sigma, Delta_T, Delta_P)

def f(l1x, l1y, l2x, l2y, type1, type2):
    # angle is defined as: angle of l1 vec - angle of l2 vec
    # L = l1 + l2

    l1 = np.sqrt(l1x**2 + l1y**2)
    l2 = np.sqrt(l2x**2 + l2y**2)
    Ll1 = l1**2 + l2x*l1x + l2y*l1y
    Ll2 = l1x*l2x + l1y*l2y + l2**2
    angle12 = np.arctan2(l1y, l1x) - np.arctan2(l2y, l2x)

    if type1 == 'T' and type2 == 'T':
        return cmbps(l1, 'T', 'T') * Ll1 + cmbps(l2, 'T', 'T') * Ll2
    if (type1 == 'T' and type2 == 'E') or (type1 == 'E' and type2 == 'T'):
        return cmbps(l1, 'T', 'E') * np.cos(angle12) * Ll1 + cmbps(l2, 'T', 'E') * Ll2
    if (type1 == 'T' and type2 == 'B') or (type1 == 'B' and type2 == 'T'):
        return cmbps(l1, 'T', 'E') * np.sin(2 * angle12) * Ll1
    if type1 == 'E' and type2 == 'E':
        return (cmbps(l1, 'E', 'E') * Ll1 + cmbps(l2, 'E', 'E') * Ll2 ) * np.cos(2 * angle12)
    if (type1 == 'E' and type2 == 'B') or (type1 == 'B' and type2 == 'E'):
        return (cmbps(l1, 'E', 'E') * Ll1 - cmbps(l2, 'B', 'B') * Ll2 ) * np.sin(2 * angle12)
    if type1 == 'B' and type2 == 'B':
        return (cmbps(l1, 'B', 'B') * Ll1 + cmbps(l2, 'B', 'B') * Ll2 ) * np.cos(2 * angle12)
    else:
        print('fucked up types')
        pass

def F(l1x, l1y, l2x, l2y, type1, type2, sigma, Delta_T, Delta_P):
    l1 = np.sqrt(l1x**2 + l1y**2)
    l2 = np.sqrt(l2x**2 + l2y**2)
    # angle is defined as: angle of l1 vec - angle of l2 vec
    if type1 == type2:
        numerator = f(l1x, l1y, l2x, l2y, type1, type2)
        denominator = 2 * cmbps_obs(l1, type1, type1, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type1, sigma, Delta_T, Delta_P)
    elif (type1 == 'T' and type2 == 'B') or (type1 == 'E' and type2 == 'B'):
        numerator = f(l1x, l1y, l2x, l2y, type1, type2)
        denominator = cmbps_obs(l1, type1, type1, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type2, type2, sigma, Delta_T, Delta_P)
    else:
        numerator = cmbps_obs(l1, type2, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type1, sigma, Delta_T, Delta_P) * f(l1x, l1y, l2x, l2y, type1, type2) - cmbps_obs(l1, type1, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type2, sigma, Delta_T, Delta_P) * f(l2x, l2y, l1x, l1y, type1, type2)
        denominator = cmbps_obs(l1, type1, type1, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type2, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l1, type2, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type1, sigma, Delta_T, Delta_P) - (cmbps_obs(l1, type1, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type2, sigma, Delta_T, Delta_P))**2

    result = numerator * denominator**(-1)
            
    return result

def l2mag(L, l1, angle):
    return np.sqrt(L**2 - 2 * L * l1 * np.cos(angle) + l1**2)

# def l1l2ang(L, l1, angle):
#     return angle - np.arccos((L - l1 * np.cos(angle)) / np.sqrt(L**2 - 2 * L * l1 * np.cos(angle) + l1**2))

def l1l2ang(L, l1, angle):
    denom = l2mag(L, l1, angle)
    # Avoid division by zero:
    if denom == 0:
        return angle
    val = (L - l1 * np.cos(angle)) / denom
    clipped_val = np.clip(val, -1, 1)
    return angle - np.arccos(clipped_val)

# def A(L, type1, type2):
#     integrand = lambda angle, l1 : l1 * f(l1, l2mag(L, l1, angle), l1l2ang(L, l1, angle), type1, type2) * F(l1, l2mag(L, l1, angle), l1l2ang(L, l1, angle), type1, type2)
#     integrand2 = lambda l1 : quad(integrand, 0, 2 * np.pi, (l1,), epsrel = 1e-2, epsabs = 1e-2)[0]
#     integral = quad(integrand2, 0, 3000, epsrel = 1e-2, epsabs = 1e-2)[0]

#     return L**2 * (2 * np.pi)**2 / integral

# dblquad impl with cartesian coords

lmaxintdef = 5000

def A_integrand(l1x, l1y, L, type1, type2, sigma, Delta_T, Delta_P):
    # Avoid potential issues at the origin:
    if np.sqrt(l1x**2 + l1y**2) == 0:
        return 0
    # l2 = L - l1, we take L to be in x direction wlog
    l2x = L - l1x 
    l2y = 0. - l1y
    return f(l1x, l1y, l2x, l2y, type1, type2) * F(l1x, l1y, l2x, l2y, type1, type2, sigma, Delta_T, Delta_P)

def A(L, type1, type2, sigma, Delta_T, Delta_P, lmaxint = lmaxintdef):
    print(L, type1, type2)
    integrand = lambda lx, ly : A_integrand(lx, ly, L, type1, type2, sigma, Delta_T, Delta_P)
    
    # Integrate over the square domain
    integral, err = dblquad(integrand, -1*lmaxint, lmaxint,
                            lambda lx: -1*lmaxint, lambda lx: lmaxint,
                            epsrel=1e-1, epsabs=1e-1)
    
    # The estimator normalization A(L) is then given by:
    # A(L) = L^2 (2π)^2 / (∫ dℓ_x dℓ_y [f*F])
    return L**2 * (2 * np.pi)**2 / integral

# nested quad impl

# def N(L, type11, type12, type21, type22):
#     integrand = lambda angle, l1 : l1 * F(l1, l2mag(L, l1, angle), l1l2ang(L, l1, angle), type11, type12) * (
#         F(l1, l2mag(L, l1, angle), l1l2ang(L, l1, angle), type21, type22) * cmbps_obs(l1, type11, type21) * cmbps_obs(l2mag(L, l1, angle), type12, type22) + F(l2mag(L, l1, angle), l1, -1*l1l2ang(L, l1, angle), type21, type22) * cmbps_obs(l1, type11, type22) * cmbps_obs(l2mag(L, l1, angle), type22, type21)
#     )
#     integrand2 = lambda l1 : quad(integrand, 0, 2 * np.pi, (l1,), epsrel = 1e-2, epsabs = 1e-2)[0]
#     integral = quad(integrand2, 0, 3000, epsrel = 1e-2, epsabs = 1e-2)[0]

#     return L**(-2) * A(L, type11, type12) * A(L, type21, type22) * (2 * np.pi)**(-2) * integral

# dblquad impl

# def N(L, type11, type12, type21, type22):
#     integrand = lambda angle, l1: l1 * F(l1, l2mag(L, l1, angle), l1l2ang(L, l1, angle), type11, type12) * (
#         F(l1, l2mag(L, l1, angle), l1l2ang(L, l1, angle), type21, type22) * cmbps_obs(l1, type11, type21) * cmbps_obs(l2mag(L, l1, angle), type12, type22)
#         + F(l2mag(L, l1, angle), l1, -1 * l1l2ang(L, l1, angle), type21, type22) * cmbps_obs(l1, type11, type22) * cmbps_obs(l2mag(L, l1, angle), type22, type21)
#     )
#     # dblquad integrates: ∫[l1=0,3000] ∫[angle=0,2π] f(angle, l1) d(angle) d(l1)
#     integral, err = dblquad(integrand, 0, 3000, lambda l1: 0, lambda l1: 2 * np.pi,
#                             epsrel=1e-2, epsabs=1e-2)
    
#     return L**(-2) * A(L, type11, type12) * A(L, type21, type22) * (2 * np.pi)**(-2) * integral

# dblquad impl in cartesian coords

def N_integrand(l1x, l1y, L, type11, type12, type21, type22, sigma, Delta_T, Delta_P):
    l2x = L - l1x 
    l2y = 0. - l1y

    l1 = np.sqrt(l1x**2 + l1y**2)
    l2 = np.sqrt(l2x**2 + l2y**2)

    term1 = (F(l1x, l1y, l2x, l2y, type11, type12, sigma, Delta_T, Delta_P) *
            F(l1x, l1y, l2x, l2y, type21, type22, sigma, Delta_T, Delta_P) *
            cmbps_obs(l1, type11, type12, sigma, Delta_T, Delta_P) *
            cmbps_obs(l2, type21, type22, sigma, Delta_T, Delta_P))
    term2 = (F(l1x, l1y, l2x, l2y, type11, type12, sigma, Delta_T, Delta_P) * 
            F(l2x, l2y, l1x, l1y, type21, type22, sigma, Delta_T, Delta_P) *
            cmbps_obs(l1, type11, type22, sigma, Delta_T, Delta_P) *
            cmbps_obs(l2, type21, type12, sigma, Delta_T, Delta_P))
    
    return term1 + term2

def N(L, type11, type12, type21, type22, sigma, Delta_T, Delta_P, lmaxint = lmaxintdef):
    print(L, type11, type12, type21, type22)
    if type11 == type21 and type12 == type22:
        return A(L, type11, type12, sigma, Delta_T, Delta_P, lmaxint = lmaxint)
    
    else:
        integrand = lambda lx, ly : N_integrand(lx, ly, L, type11, type12, type21, type22, sigma, Delta_T, Delta_P)

        # Integrate over lx and ly in the square [-3000, 3000] x [-3000, 3000]
        integral, err = dblquad(integrand, -1*lmaxint, lmaxint,
                                lambda lx: -1*lmaxint, lambda lx: lmaxint,
                                epsrel=1e-1, epsabs=1e-1)
        
        return L**(-2) * A(L, type11, type12, sigma, Delta_T, Delta_P) * A(L, type21, type22, sigma, Delta_T, Delta_P) * (2 * np.pi)**(-2) * integral

#############################################################################################
configs = (('T', 'T'), ('T', 'E'), ('T', 'B'), ('E', 'E'), ('E', 'B'), ('B', 'B'))
# configs = (('E', 'E'), ('E', 'B'), ('B', 'B'))
#############################################################################################

num_configs = len(configs)

config_mat = lambda L, sigma, Delta_T, Delta_P, lmaxint : np.array([[N(L, *config1, *config2, sigma, Delta_T, Delta_P, lmaxint = lmaxint) for config1 in configs] for config2 in configs])

def lps_noise(L, sigma, Delta_T, Delta_P, lmaxint = lmaxintdef):
    return 1 / np.sum(np.linalg.inv(config_mat(L, sigma, Delta_T, Delta_P, lmaxint)))

if __name__ == '__main__':
    print('kiki')
    import numpy as np
    import os
    #import multiprocessing
    #from functools import partial
    # Create a logarithmically spaced array for L.
    Ls = np.logspace(np.log10(50), np.log10(2000), 16)
    np.savetxt("cmb_noise_files/ls.txt", Ls)
    
    def create_noise_vals(sigma, Delta_T, Delta_P):
        # Use multiprocessing to compute N for each L in parallel.
        # with multiprocessing.Pool() as pool:
        #     func = partial(lps_noise, sigma=sigma, Delta_T=Delta_T, Delta_P=Delta_P)
        #     noise_vals = pool.map(func, Ls)
        noise_vals = [lps_noise(L, sigma, Delta_T, Delta_P) for L in Ls]
        
        # Ensure the output directory exists.
        out_dir = "cmb_noise_files"
        os.makedirs(out_dir, exist_ok=True)
        
        # Filename includes sigma, Delta_T, and Delta_P.
        filename = f"{out_dir}/Ns_sigma{sigma}_DeltaT{Delta_T}_DeltaP{Delta_P}.txt"
        np.savetxt(filename, noise_vals)
        
        print(f'\nParameters: sigma={sigma}, Delta_T={Delta_T}, Delta_P={Delta_P}\nComputed Ns:\n', noise_vals)
    
    # Noise configurations: (sigma, Delta_T, Delta_P)
    # noise_configs = [
    #     (1, 0, 6),  # fsky = 0.5
    #     (1, 0, 3),  # fsky = 0.05
    #     (3, 0, 1),  # fsky = 0.5
    # ]

    noise_configs = [
        (4, 1, np.sqrt(2)),  # fsky = 0.5
    ]

    
    # Loop over the noise configurations.
    for sigma, Delta_T, Delta_P in noise_configs:
        create_noise_vals(sigma, Delta_T, Delta_P)
