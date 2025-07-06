import data_importer_new as di
import numpy as np
import matplotlib.pyplot as plt

ls = np.arange(5, 3000)

sigma = 3
Delta_T = 1

cmbps_TT = np.array([di.lps_der_test(l, b't', b't', b'snr') for l in ls])
cmbps_TT_noise_namikawa = np.array([di.cmbps_noise_test(l, sigma, Delta_T) for l in ls])
cmbps_TT_noise_hu = np.array([di.cmbps_noise_test(l, 4, 1) for l in ls])


plt.loglog(ls, ls**2 * cmbps_TT, label = 'C^TT')
plt.loglog(ls, ls**2 * cmbps_TT_noise_namikawa, label = f'noise, sigma = {sigma}, Delta_T = {Delta_T}')
plt.loglog(ls, ls**2 * cmbps_TT_noise_hu, label = f'noise, sigma = {4}, Delta_T = {1}')

plt.legend()
plt.savefig('cmbps_vs_noise.png')
