import cosm_setup as cs
import numpy as np
import matplotlib.pyplot as plt

spectra = cs.lensing_spectra()

ls = np.arange(10, 2001, 100)
lbss = np.array([spectra.lbs(l, l, l, ('c', 'c', 'c')) for l in ls])
print(lbss)

plt.loglog(ls, lbss * ls**4 / (2 * np.pi)**2 * (ls * (ls + 1) / 2) ** 3)
plt.savefig('test.png')