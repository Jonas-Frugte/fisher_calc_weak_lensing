import data_importer as di
import cosm_setup as cs
import numpy as np

ls = np.arange(2, 2001, 2)

lbs_eq_int = np.abs(np.array([
    di.lbs_der_py(l, l, l, b'c', b'c', b'c', 200, b'snr') * float(l)**10 / (2 * np.pi)**2 / 8 for l in ls
]))

lbs_fold_int = np.abs(np.array([
    di.lbs_der_py(l, 2*l, 2*l, b'c', b'c', b'c', 200, b'snr') * float(l)**10 / (2 * np.pi)**2 / 8 for l in ls
]))

spectra = cs.lensing_spectra()

ls_nz_samp_eq = ls[np.nonzero(lbs_eq_int)][0::5]
ls_nz_samp_fold = ls[np.nonzero(lbs_fold_int)][0::5]

lbs_eq_ex = np.array([
    spectra.lbs(l, l, l, ('c', 'c', 'c')) * float(l)**10 / (2 * np.pi)**2 / 8 for l in ls_nz_samp_fold
])

lbs_fold_ex = np.array([
    spectra.lbs(l, 2*l, 2*l, ('c', 'c', 'c')) * float(l)**10 / (2 * np.pi)**2 / 8 for l in ls_nz_samp_fold
])

import matplotlib.pyplot as plt

plt.loglog(ls, lbs_eq_int)
plt.loglog(ls, lbs_fold_int)
plt.loglog(ls_nz_samp_eq, lbs_eq_ex, linestyle = '--')
plt.loglog(ls_nz_samp_fold, lbs_fold_ex, linestyle = '--')

plt.savefig('plots/bispec_int_check.png')