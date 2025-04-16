import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Define data
ell_max = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]

snrb_c_s3 = [15.07615989323126, 60.14746748173047, 103.12636014205, 129.35899320252307, 142.68006253688625, 149.21360488222513, 152.46213942042544, 154.16826858416584]

snrb_s_s3 = [531.1622436479552, 1725.2752469768607, 2949.2997415134773, 4030.1790464195456, 4935.08531114439, 5676.463829544848, 6279.344472866087, 6769.365663353458]

snrb_c_s4 = [23.083226199124194, 160.95563855174456, 453.86854903039114, 846.2154979047796, 1236.5912661335751, 1556.0849450590474, 1785.694771513511, 1939.702697620804]

snrb_s_s4 = [1841.1019217665996, 13139.433138837652, 35830.00873208229, 66306.94638041034, 100272.43117960829, 134537.47452003675, 167174.8737905056, 197196.08660044952]

ls = np.linspace(2, 2000, 201)
plot_datas = np.zeros((4, len(ls)))
datas = (snrb_c_s3, snrb_s_s3, snrb_c_s4, snrb_s_s4)
for i in range(4):
    data = datas[i]
    interp_func = interp1d([0] + ell_max, np.sqrt([0] + data), kind='cubic')
    plot_datas[i, :] = [interp_func(l) for l in ls]

fig, axs = plt.subplots(1, 2, figsize = (8.3, 3))

# Plot each array
axs[1].plot(ls, plot_datas[0], label='CMB, stage 3', linestyle = '--', color='red')
axs[0].plot(ls, plot_datas[1], label='galaxy, stage 3', linestyle = '--', color='blue')
axs[1].plot(ls, plot_datas[2], label='CMB, stage 4', color='red')
axs[0].plot(ls, plot_datas[3], label='galaxy, stage 4', color='blue')

axs[0].minorticks_on()
axs[1].minorticks_on()

axs[0].set_xlabel(r'$l_{\max}$', fontsize=12)
axs[0].set_ylabel(r'$S/N$', fontsize=12)
axs[1].set_xlabel(r'$l_{\max}$', fontsize=12)
# axs[1].set_ylabel(r'$S/N$', fontsize=12)

# Increase the x-axis limit to add extra space on the right
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
plt.savefig('/Users/jonasfrugte/Desktop/Research_Project/fisher_calc_weak_lensing/paper/figures/snrplots.pdf', dpi=300)
# plt.show()
