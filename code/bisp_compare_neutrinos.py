import cosm_setup as cs



import numpy as np
import matplotlib.pyplot as plt

ls = np.arange(10, 2001, 20)

for mnu in (0, 0.75*0.06, 0.9*0.06, 0.95 * 0.06, 0.06, 1.05*0.06, 1.1*0.06, 1.25*0.06, 0.2):
    spectra = cs.lensing_spectra(mnu = mnu)

    lbs_eq_ex_nonl = np.array([
        spectra.lbs_flat(l, l, l, ('c', 'c', 'c')) * float(l)**10 / (2 * np.pi)**2 / 8 for l in ls
    ])

    lbs_fold_ex_nonl = np.array([
        spectra.lbs_flat(l, l/2, l/2, ('c', 'c', 'c')) * float(l)**10 / 2**4 / (2 * np.pi)**2 / 8 for l in ls
    ])

    plt.loglog(ls, lbs_eq_ex_nonl, label=f'eq, mnu={mnu}')
    plt.loglog(ls, lbs_fold_ex_nonl, linestyle = 'dotted', label=f'fold, mnu={mnu}')

plt.legend()
plt.savefig('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/toshiya_bisp_curves/plot.png', dpi = 300)