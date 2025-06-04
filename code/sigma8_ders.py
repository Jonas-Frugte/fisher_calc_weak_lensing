# first arguments of spectra class: H0=67.4, ombh2=0.0224, omch2=0.120, ns=0.965, As=2.1e-9, mnu=0.06, w0 = -1., wa = 0.
import numpy as np
import cosm_setup as cs
# this needs to be in the same order as used for fisher matrix calculations, otherwise you will die later trying to fix param constraints :)
par_names = ('H', 'ombh2', 'omch2', 'ns', 'mnu', 'tau', 'As', 'w0')
fiducial_cosm_par = (67.4, 0.0224, 0.120, 0.965, 2.1e-9, 0.06, 0.06, -1)
num_pars = len(par_names)

dx = 0.05 

dx_coeffs = [2, 1, -1, -2]

def dx_coeffs_to_str(coeff):
    if coeff == 2:
        return '2p'
    if coeff == 1:
        return '1p'
    if coeff == -1:
        return '1m'
    if coeff == -2:
        return '2m'
    else:
        raise Exception

sigma8s = np.zeros((num_pars, len(dx_coeffs)))

for j in range(len(dx_coeffs)):
    for i in range(num_pars):
        perturbed_params = [fiducial_cosm_par[k] * (1 + dx_coeffs[j] * dx * int(k == i)) for k in range(num_pars)]
        spectra = cs.lensing_spectra(*perturbed_params)
        sigma8s[i, j] = spectra.results.get_sigma8_0()

print(sigma8s)

ders = np.zeros(num_pars)
for i in range(num_pars):
    ders[i] = (-1 * sigma8s[i, 0] + 8 * sigma8s[i, 1] - 8 * sigma8s[i, 2] + sigma8s[i, 3]) / (12 * dx * fiducial_cosm_par[i])

print(ders)

# IMPORTANT!
# the derivatives are gonna be in order: (H0, ombh2, omch2, ns, As, tau, mnu, w0),
# but for the fisher matrices they need to be switched around to be in order: ('H', 'ombh2', 'omch2', 'ns', 'mnu', 'tau', 'As', 'w0')
# mnu and As need to be switched basically

