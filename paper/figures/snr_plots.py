import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Define data
ell_max = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]

snrb_c_s3 = [0.4606180302137393, 2.428277643136825, 4.556412710716516, 5.924272288700067, 6.620955446735307, 6.958473094521213, 7.12003774130501, 7.201717334076159]

snrb_s_s3 = [14.937583099693954, 44.263581102561155, 70.4841150784228, 91.60588005582235, 108.3741775299598, 121.44133496518877, 131.79346525781503, 139.98497349255402]

snrb_c_s4 = [23.083226199124194, 160.95563855174456, 453.86854903039114, 846.2154979047796, 1236.5912661335751, 1556.0849450590474, 1785.694771513511, 1939.702697620804]

snrb_s_s4 = [1841.1019217665996, 13139.433138837652, 35830.00873208229, 66306.94638041034, 100272.43117960829, 134537.47452003675, 167174.8737905056, 197196.08660044952]

ls = np.linspace(2, 2000, 201)
plot_datas = np.zeros((4, len(ls)))
datas = (snrb_c_s3, snrb_s_s3, snrb_c_s4, snrb_s_s4)
for i in range(4):
    data = datas[i]
    interp_func = interp1d([0] + ell_max, np.sqrt([0] + data), type='cubic')
    plot_datas[i, :] = [interp_func(l) for l in ls]


# Set figure size
plt.figure(figsize=(7,5))

# Plot each array
plt.plot(ls, plot_datas[0], marker='x', label='CMB, stage 3', linestyle = '--', color='red')
plt.plot(ls, plot_datas[1], marker='x', label='galaxy, stage 3', linestyle = '--', color='blue')
plt.plot(ls, plot_datas[2], marker='x', label='CMB, stage 4', color='red')
plt.plot(ls, plot_datas[3], marker='x', label='galaxy, stage 4', color='blue')

plt.minorticks_on()
plt.xlabel(r'$l_{\max}$', fontsize=12)
plt.ylabel(r'$S/N$', fontsize=12)

# Increase the x-axis limit to add extra space on the right
x_max = ell_max[-1]
plt.xlim(ell_max[0] - 50, x_max + 250)

# Annotate each line with the y-value at the right end
offset = 50  # horizontal offset for the text
plt.text(x_max + offset, np.sqrt(snrb_c_s3)[-1], f'{np.sqrt(snrb_c_s3)[-1]:.2f}', color='black', ha='left', va='center')
plt.text(x_max + offset, np.sqrt(snrb_s_s3)[-1], f'{np.sqrt(snrb_s_s3)[-1]:.2f}', color='black', ha='left', va='center')
plt.text(x_max + offset, np.sqrt(snrb_c_s4)[-1], f'{np.sqrt(snrb_c_s4)[-1]:.2f}', color='black', ha='left', va='center')
plt.text(x_max + offset, np.sqrt(snrb_s_s4)[-1], f'{np.sqrt(snrb_s_s4)[-1]:.2f}', color='black', ha='left', va='center')

plt.legend(fontsize=10, loc='upper left')
plt.tight_layout()

# Save the figure
plt.savefig('/Users/jonasfrugte/Desktop/Research_Project/paper/figures/snrplots.pdf', dpi=300)
# plt.show()
