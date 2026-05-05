import numpy as np
import matplotlib.pyplot as plt
from cosm_setup import lensing_spectra

cosmo = lensing_spectra()

k_h = np.logspace(np.log10(0.01), np.log10(0.5), 300)  # h/Mpc
k   = k_h * cosmo.h                                    # Mpc^-1

n        = np.array([cosmo.ns_eff(ki, 0.0)          for ki in k])
n_smooth = np.array([cosmo.smoothed_ns_eff(ki, 0.0) for ki in k])

plt.figure(figsize=(6, 4))
plt.plot(k_h, n,        'k--', lw=1,   label='unsmoothed')
plt.plot(k_h, n_smooth, 'r-',  lw=1.5, label='smoothed')
plt.xscale('log')
plt.xlim(0.01, 0.5)
plt.ylim(-3, 0.5)
plt.xlabel(r'$k\ [h/\mathrm{Mpc}]$')
plt.ylabel(r'$n(k)$')
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('n_k_smoothed.png', dpi=150)
plt.show()