import cosm_setup as cs

spectra = cs.lensing_spectra(
    ombh2=0.0223,
    omch2=0.119 - 0.0223,
    As=2.13e-9,
    ns=0.965,
    fiducial_k_nls = True)

import numpy as np
import matplotlib.pyplot as plt

ls = np.arange(10, 2001, 2)

lbs_eq_ex = np.array([
    spectra.lbs_tree(l, l, l, ('c', 'c', 'c')) * float(l)**10 / (2 * np.pi)**2 / 8 for l in ls
])

lbs_fold_ex = np.array([
    spectra.lbs_tree(l, l/2, l/2, ('c', 'c', 'c')) * float(l)**10 / 2**4 / (2 * np.pi)**2 / 8 for l in ls
])

bisp_eq_tosh = np.loadtxt('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/lbs_eq_tree.txt')

bisp_fold_tosh = np.loadtxt('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/lbs_fold_tree.txt')

plt.loglog(bisp_eq_tosh[:,0], bisp_eq_tosh[:,1], label = 'tosh, eq', color = 'blue')
plt.loglog(bisp_fold_tosh[:,0], bisp_fold_tosh[:,1], label = 'tosh, fold', color = 'red')

plt.loglog(ls, lbs_eq_ex, label = 'eq', linestyle = '--', color = 'blue')
plt.loglog(ls, lbs_fold_ex, label = 'fold', linestyle = '--', color = 'red')

plt.legend()
plt.savefig('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/plot.png', dpi = 300)