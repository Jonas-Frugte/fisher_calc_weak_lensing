import matplotlib.pyplot as plt
import numpy as np
import cosm_setup
import data_importer

spectra = cosm_setup.lensing_spectra()

fig, ax = plt.subplots(2, figsize=(12, 14))

ls = np.arange(50, 491, 10)

types = (b's', b's', b's')

# Equilateral configuration
exact_data_1 = [spectra.lbs(l, l, l, ('s', 's', 's')) for l in ls]
interp_data_1 = [data_importer.lbs_f(l, l, l, *types, 5000) for l in ls]

exact_data_2 = [spectra.lbs(2*l, l, l, ('s', 's', 's')) for l in ls]
interp_data_2 = [data_importer.lbs_f(2*l, l, l, *types, 500) for l in ls]

# Plotting the data for each configuration in the current subplot
ax[0].loglog(ls, exact_data_1, linestyle='--', label='Exact 1', marker = 'x')
# ax[0].loglog(ls, interp_data_1, label='Interpolated 1')

ax[0].loglog(ls, exact_data_2, linestyle='--', label='Exact 2', marker = 'x')
# ax[0].loglog(ls, interp_data_2, label='Interpolated 2')

# ax[1].semilogx(ls,
#             [interp_data_1[i]/exact_data_1[i] for i in range(len(ls))], label = '1')
# ax[1].semilogx(ls,
#             [interp_data_2[i]/exact_data_2[i] for i in range(len(ls))], label = '2')

ax[0].legend()
ax[1].legend()

plt.savefig('/home3/p319950/ResearchProject/plots_and_figures/interp_accuracy.png')
print('done plotting')
#fig.savefig('/Users/jonasfrugte/Desktop/ResearchProject/SNR_calculation/plots_and_figures/poopie.png')