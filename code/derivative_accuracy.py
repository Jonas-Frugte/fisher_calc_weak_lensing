import data_importer_new as di
import numpy as np
import matplotlib.pyplot as plt

ls = np.arange(2, 2000, 1)
ls_even = np.arange(2, 2000, 2)

fig, axs = plt.subplots(7, 2, figsize=(8.3 * 2, 11.7 * 2))
pars = [b'H', b'ombh2', b'omch2', b'As', b'mnu', b'ns', b'w0']
pars_latex = [r'$H_0$', r'$\Omega_b h^2$', r'$\Omega_c h^2$', r'$A_s$', r'$m_{\nu}$', r'$n_s$', r'$w_0$']


dds = [-2, -1, 0, 1, 2]
colors = ['darkred', 'red', 'black', 'blue', 'darkblue']

for i in range(len(pars)):
    for dd_index in range(len(dds)):
        axs[i,0].semilogx(ls_even, [di.lbs_der_test(l, l, l, b'c', b'c', b'c', 200, pars[i], dds[dd_index]) / di.lbs_der_test(l, l, l, b'c', b'c', b'c', 200, pars[i], 0) for l in ls_even], label = f'{dds[dd_index]}', lw = 1, color = colors[dd_index])
        axs[i,0].set_ylabel(r'$\partial B^{\psi_{\text{CMB}}}_{lll} / \partial$' + pars_latex[i])

        axs[i,1].semilogx(ls, [di.lps_der_test(l, b'c', b'c', pars[i], dds[dd_index]) / di.lps_der_test(l, b'c', b'c', pars[i], 0) for l in ls], label = f'{pars[i]}, {dds[dd_index]}', lw = 1, color = colors[dd_index])
        axs[i,1].set_ylabel(r'$\partial C^{\psi_{\text{CMB}}}_{l} / \partial$' + pars_latex[i])

axs[-1, 0].set_xlabel(r'$l$')
axs[-1, 1].set_xlabel(r'$l$')

axs[5, 0].set_ylim(0.998, 1.002)
axs[5, 1].set_ylim(0.998, 1.002)
axs[6, 1].set_ylim(0.995, 1.005)
axs[2, 1].set_ylim(0.998, 1.002)
axs[1, 1].set_ylim(0.99, 1.01)


handles, labels = axs[0,0].get_legend_handles_labels()
fig.legend(handles, [r'-10%', r'-5%', r'0%', r'5%', r'10%'], loc = 'right')

fig.savefig('plots/derivative_accuracy.png')