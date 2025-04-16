import cosm_setup as cs
spectra = cs.lensing_spectra()

import data_importer as di


import numpy as np
import matplotlib.pyplot as plt

ls = np.arange(10, 2001, 20)

# di.lps_f_obs_test gives lensing powerspectrum plus noise spectrum of lensing potential
# 'c' is for cmb lensing

bisp_cov_eq = np.array([
    float(l)**4 * di.lps_f_obs_test(l, b'c', b'c')**3 * float(l)**6 / (2 * np.pi)**2 / 8 for l in ls
])

bisp_cov_fold = np.array([
    float(l)**4 * di.lps_f_obs_test(l, b'c', b'c') * di.lps_f_obs_test(l/2, b'c', b'c')**2 * float(l)**6 / 2**4 / (2 * np.pi)**2 / 8 for l in ls
])

plt.loglog(ls, bisp_cov_eq, label = 'bisp cov, eq')
plt.loglog(ls, bisp_cov_fold, label = 'bisp cov, fold')

plt.ylabel(r'$l^4(C^\kappa_{l_1} + N_{l_1})(C^\kappa_{l_2} + N_{l_2})(C^\kappa_{l_3} + N_{l_3})$')
plt.xlabel(r'$l$')

plt.legend()
plt.savefig('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/bisp_cov.png', dpi = 300)