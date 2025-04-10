import cosm_setup as cs

spectra = cs.lensing_spectra()
# spectra = cs.lensing_spectra(
#     # ombh2=0.0223,
#     # omch2=0.119,
#     # As=2.13e-9,
#     # ns=0.965,
#     )

import numpy as np
import matplotlib.pyplot as plt

ls = np.arange(10, 2001, 20)

lbs_eq_ex_tree = np.array([
    spectra.lbs_flat(l, l, l, ('c', 'c', 'c'), model = 'tree') * float(l)**10 / (2 * np.pi)**2 / 8 for l in ls
])

lbs_fold_ex_tree = np.array([
    spectra.lbs_flat(l, l/2, l/2, ('c', 'c', 'c'), model = 'tree') * float(l)**10 / 2**4 / (2 * np.pi)**2 / 8 for l in ls
])

lbs_eq_ex_nonl = np.array([
    spectra.lbs_flat(l, l, l, ('c', 'c', 'c')) * float(l)**10 / (2 * np.pi)**2 / 8 for l in ls
])

lbs_fold_ex_nonl = np.array([
    spectra.lbs_flat(l, l/2, l/2, ('c', 'c', 'c')) * float(l)**10 / 2**4 / (2 * np.pi)**2 / 8 for l in ls
])

bisp_eq_tosh_nonl = np.loadtxt('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/lbs_eq_nonl.txt')

bisp_fold_tosh_nonl = np.loadtxt('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/lbs_fold_nonl.txt')

bisp_eq_tosh_tree = np.loadtxt('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/lbs_eq_tree.txt')

bisp_fold_tosh_tree = np.loadtxt('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/lbs_fold_tree.txt')

plt.loglog(bisp_eq_tosh_tree[:,0], bisp_eq_tosh_tree[:,1], label = 'tosh, eq, tree', color = 'blue', linestyle = '--')
plt.loglog(bisp_fold_tosh_tree[:,0], bisp_fold_tosh_tree[:,1], label = 'tosh, fold, tree', color = 'red', linestyle = '--')

plt.loglog(bisp_eq_tosh_nonl[:,0], bisp_eq_tosh_nonl[:,1], label = 'tosh, eq, nonl', color = 'blue')
plt.loglog(bisp_fold_tosh_nonl[:,0], bisp_fold_tosh_nonl[:,1], label = 'tosh, fold, nonl', color = 'red')

plt.loglog(ls, lbs_eq_ex_tree, linestyle = 'dotted', color = 'black')
plt.loglog(ls, lbs_fold_ex_tree, linestyle = 'dotted', color = 'black')

plt.loglog(ls, lbs_eq_ex_nonl, linestyle = 'dotted', color = 'black')
plt.loglog(ls, lbs_fold_ex_nonl, linestyle = 'dotted', color = 'black')

plt.legend()
plt.savefig('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/plot.png', dpi = 300)