import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Define data
ell_max = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]

f_sky = 0.5

snrb_c_s3 = f_sky * np.array([12.672454984617664, 83.7377129725408, 206.95352488728747, 350.9902759829114, 498.35684081032775, 645.7108221860103, 796.0813174568198, 953.7584589205301])

snrb_s_s3 = f_sky * np.array([477.5259136914146, 2533.295459864549, 5493.551176152028, 8644.382159891089, 11636.918745628755, 14316.53223335389, 16656.336652636706, 18662.00612425752])

snrb_c_s3_pb = f_sky * np.array([7.825690349725596, 46.151237529667654, 107.18588258187486, 177.07130085865384, 247.7876624988021, 319.56645397880607, 394.29703253136023, 474.84034728265567])

snrb_s_s3_pb = f_sky * np.array([464.03174572011704, 2474.5313315586527, 5380.676014059787, 8479.734824268631, 11425.708268033992, 14064.629823584384, 16368.973212497156, 18344.038048655635])

snrb_c_s4 = f_sky * np.array([15.286201971070938, 134.55984159950066, 435.4551401432371, 923.6584823791366, 1546.5050107983147, 2239.308423902829, 2962.6579436700445, 3703.6935030765526])

snrb_s_s4 = f_sky * np.array([1037.0653293826213, 10727.141251236242, 36738.626526038235, 80019.31027627057, 137303.10313148203, 203689.07845263026, 274894.5746942937, 347161.44922220515])

snrb_c_s4_pb = f_sky * np.array([9.42808790586853, 76.53291859391692, 237.99881712534457, 499.29733575525745, 828.0210947262487, 1192.5574598676799, 1572.9647473021628, 1966.6182767589248])

snrb_s_s4_pb = f_sky * np.array([1008.5197756342941, 10508.392680065226, 36116.6775596652, 78809.14679562827, 135363.4396650899, 200927.61688691552, 271252.22563812265, 342616.16867419565])

ls = np.linspace(2, 2000, 201)
plot_datas = np.zeros((8, len(ls)))
datas = (snrb_c_s3, snrb_s_s3, snrb_c_s4, snrb_s_s4, snrb_c_s3_pb, snrb_s_s3_pb, snrb_c_s4_pb, snrb_s_s4_pb)
for i in range(len(datas)):
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

# post born curves
axs[1].plot(ls, plot_datas[4], label='CMB, stage 3, PB corr.', linestyle='--', color='blue')
axs[1].plot([ls[-1]], [plot_datas[4][-1]], marker='o', color='blue', label='_nolegend_')

# axs[0].plot(ls, plot_datas[5], label='galaxy, stage 3, PB corr.', linestyle='--', color='blue')
# axs[0].plot([ls[-1]], [plot_datas[5][-1]], marker='o', color='blue', label='_nolegend_')

axs[1].plot(ls, plot_datas[6], label='CMB, stage 4, PB corr.', linestyle='-', color='blue')
axs[1].plot([ls[-1]], [plot_datas[6][-1]], marker='o', color='blue', label='_nolegend_')

# axs[0].plot(ls, plot_datas[7], label='galaxy, stage 4, PB corr.', linestyle='-', color='blue')
# axs[0].plot([ls[-1]], [plot_datas[7][-1]], marker='o', color='blue', label='_nolegend_')


axs[0].minorticks_on()
axs[1].minorticks_on()

axs[0].set_xlabel(r'$l_{\max}$', fontsize=10)
axs[0].set_ylabel(r'$S/N$', fontsize=10)
axs[1].set_xlabel(r'$l_{\max}$', fontsize=10)
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

axs[1].text(x_max + offset, np.sqrt(snrb_c_s3_pb)[-1], f'{np.sqrt(snrb_c_s3_pb)[-1]:.0f}', color='black', ha='left', va='center')
# axs[0].text(x_max + offset, np.sqrt(snrb_s_s3_pb)[-1], f'{np.sqrt(snrb_s_s3_pb)[-1]:.0f}', color='black', ha='left', va='center')
axs[1].text(x_max + offset, np.sqrt(snrb_c_s4_pb)[-1], f'{np.sqrt(snrb_c_s4_pb)[-1]:.0f}', color='black', ha='left', va='center')
# axs[0].text(x_max + offset, np.sqrt(snrb_s_s4_pb)[-1], f'{np.sqrt(snrb_s_s4_pb)[-1]:.0f}', color='black', ha='left', va='center')

axs[0].legend(fontsize=10, loc='upper left')
axs[1].legend(fontsize=10, loc='upper left')

axs[0].grid(True, which='both', linestyle='--', linewidth=0.5)
axs[1].grid(True, which='both', linestyle='--', linewidth=0.5)

plt.tight_layout()

# Save the figure
plt.savefig('../paper/figures/snrplots.pdf')
plt.savefig('../paper/figures/snrplots.png')