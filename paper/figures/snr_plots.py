import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Define data
ell_max = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]

f_sky = 0.5

snrb_c_s3 = f_sky * np.array([18.14225774292767, 104.80963360462812, 246.33294403758597, 406.55575384332985, 565.893477899891, 723.8984338793172, 883.4983606420614, 1050.0227116349033])

snrb_s_s3 = f_sky * np.array([463.85614017322854, 1522.8225121434366, 2593.8277015312924, 3530.456621967403, 4305.844778597493, 4941.576011585971, 5457.833730995741, 5876.8344414029])

snrb_c_s4 = f_sky * np.array([21.634443589995616, 166.04498908145288, 510.24879700176723, 1052.5932963778246, 1725.518886984883, 2466.6245011624437, 3231.621124084956, 4010.7797969493395])

snrb_s_s4 = f_sky * np.array([1613.0939430053656, 11693.189598379018, 31673.989578902532, 58205.021807229314, 87412.42982093079, 116807.14803651144, 144721.82601691395, 170335.10666953446])

ls = np.linspace(2, 2000, 201)
plot_datas = np.zeros((4, len(ls)))
datas = (snrb_c_s3, snrb_s_s3, snrb_c_s4, snrb_s_s4)
for i in range(4):
    data = datas[i]
    interp_func = interp1d([0] + ell_max, np.sqrt([0] + list(data)), kind='cubic')
    plot_datas[i, :] = [interp_func(l) for l in ls]

fig, axs = plt.subplots(1, 2, figsize = (8.3, 3))

# Plot each array
axs[1].plot(ls, plot_datas[0], label='CMB, stage 3', linestyle='--', color='black')
axs[1].plot([ls[-1]], [plot_datas[0][-1]], marker='o', color='black', label='_nolegend_')

axs[0].plot(ls, plot_datas[1], label='galaxy, stage 3', linestyle='--', color='black')
axs[0].plot([ls[-1]], [plot_datas[1][-1]], marker='o', color='black', label='_nolegend_')

axs[1].plot(ls, plot_datas[2], label='CMB, stage 4', linestyle='-', color='black')
axs[1].plot([ls[-1]], [plot_datas[2][-1]], marker='o', color='black', label='_nolegend_')

axs[0].plot(ls, plot_datas[3], label='galaxy, stage 4', linestyle='-', color='black')
axs[0].plot([ls[-1]], [plot_datas[3][-1]], marker='o', color='black', label='_nolegend_')


axs[0].minorticks_on()
axs[1].minorticks_on()

axs[0].set_xlabel(r'$l_{\max}$', fontsize=12)
axs[0].set_ylabel(r'$S/N$', fontsize=12)
axs[1].set_xlabel(r'$l_{\max}$', fontsize=12)
# axs[1].set_ylabel(r'$S/N$', fontsize=12)

# Increase the x-axis limit to add extra space on the right
axs[0].set_ylim(0, np.sqrt(snrb_s_s4[-1]) * 1.125)
axs[1].set_ylim(0, np.sqrt(snrb_c_s4[-1]) * 1.125)

x_max = ell_max[-1]
axs[0].set_xlim(2, x_max + 250)
axs[1].set_xlim(2, x_max + 250)

# Annotate each line with the y-value at the right end
offset = 50  # horizontal offset for the text
axs[1].text(x_max + offset, np.sqrt(snrb_c_s3)[-1], f'{np.sqrt(snrb_c_s3)[-1]:.0f}', color='black', ha='left', va='center')
axs[0].text(x_max + offset, np.sqrt(snrb_s_s3)[-1], f'{np.sqrt(snrb_s_s3)[-1]:.0f}', color='black', ha='left', va='center')
axs[1].text(x_max + offset, np.sqrt(snrb_c_s4)[-1], f'{np.sqrt(snrb_c_s4)[-1]:.0f}', color='black', ha='left', va='center')
axs[0].text(x_max + offset, np.sqrt(snrb_s_s4)[-1], f'{np.sqrt(snrb_s_s4)[-1]:.0f}', color='black', ha='left', va='center')

axs[0].legend(fontsize=10, loc='upper left')
axs[1].legend(fontsize=10, loc='upper left')

plt.tight_layout()

# Save the figure
plt.savefig('/Users/jonasfrugte/Desktop/fisher_calc_weak_lensing/paper/figures/snrplots.pdf', dpi=300)
# plt.show()
