import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Define data
ell_max = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]

f_sky = 0.5

snrb_c_s3 = f_sky * np.array([4.95584804530861, 57.754982819045495, 154.05291503946256, 276.6628467874657, 397.4145481594399, 523.771418830852, 647.0739636714212, 778.3940264209907])

snrb_s1_s3 = f_sky * np.array([27.801064868683856, 55.123772937516094, 75.81743181970205, 87.82112428098876, 96.4168771831648, 101.8658328529791, 105.92308290115169, 108.64057393249895])  # galaxy bin 1 (closest)

snrb_s4_s3 = f_sky * np.array([55.4390243747189, 363.9721039655066, 735.5774592936141, 1106.6013974001628, 1433.4976445062948, 1728.9939419158582, 1978.6783733051618, 2199.778813954112])  # galaxy bin 4 (furthest)

snrb_c_s4 = f_sky * np.array([5.955877546302534, 93.94237095629144, 328.09090975260284, 745.871071776269, 1262.1890643224106, 1862.9246282461932, 2464.324622132588, 3089.590206305989])

snrb_s1_s4 = f_sky * np.array([485.7564893208991, 2041.2067046903053, 3649.7777417494276, 4899.77759278681, 5828.983891164869, 6508.791709394715, 7005.24789655677, 7370.969222935724])  # galaxy bin 1 (closest)

snrb_s4_s4 = f_sky * np.array([158.40880895872243, 2301.3245662989275, 7591.28696254538, 16066.547921560106, 26332.20898433008, 38017.25724836697, 49751.827456985215, 61650.04365787259])  # galaxy bin 4 (furthest)

snrb_c_s3_pb = f_sky * np.array([4.10101637964598, 31.0042831168067, 81.02172018011971, 136.23185686655094, 195.76751930922507, 253.40093374243162, 315.6601621979454, 379.37620141325164])

snrb_c_s4_pb = f_sky * np.array([5.039317946715429, 52.155545330656814, 185.30449559273902, 394.70817783613717, 674.9949288792279, 973.4328294931372, 1294.17170077849, 1611.4446809847316])

ls = np.linspace(2, 2000, 201)
plot_datas = np.zeros((10, len(ls)))  # increased to 10 for new curves
datas = (snrb_c_s3, snrb_s1_s3, snrb_s4_s3, snrb_c_s4, snrb_s1_s4, snrb_s4_s4, snrb_c_s3_pb, snrb_c_s4_pb)
for i in range(len(datas)):
    data = datas[i]
    interp_func = interp1d([0] + ell_max, np.sqrt([0] + list(data)), kind='cubic')
    plot_datas[i, :] = [interp_func(l) for l in ls]

fig, axs = plt.subplots(1, 2, figsize = (8.3, 3))

# Plot each array
axs[1].plot(ls, plot_datas[0], label='CMB, stage 3', linestyle='--', color='black')
axs[1].plot([ls[-1]], [plot_datas[0][-1]], marker='o', color='black', label='_nolegend_')

axs[0].plot(ls, plot_datas[1], label='galaxy bin 1, stage 3', linestyle='--', color='black')
axs[0].plot([ls[-1]], [plot_datas[1][-1]], marker='o', color='black', label='_nolegend_')

axs[0].plot(ls, plot_datas[2], label='galaxy bin 4, stage 3', linestyle='--', color='blue')
axs[0].plot([ls[-1]], [plot_datas[2][-1]], marker='o', color='blue', label='_nolegend_')

axs[1].plot(ls, plot_datas[3], label='CMB, stage 4', linestyle='-', color='black')
axs[1].plot([ls[-1]], [plot_datas[3][-1]], marker='o', color='black', label='_nolegend_')

axs[0].plot(ls, plot_datas[4], label='galaxy bin 1, stage 4', linestyle='-', color='black')
axs[0].plot([ls[-1]], [plot_datas[4][-1]], marker='o', color='black', label='_nolegend_')

axs[0].plot(ls, plot_datas[5], label='galaxy bin 4, stage 4', linestyle='-', color='blue')
axs[0].plot([ls[-1]], [plot_datas[5][-1]], marker='o', color='blue', label='_nolegend_')

# post born curves
axs[1].plot(ls, plot_datas[6], label='CMB, stage 3, \n post-Born corr.', linestyle='--', color='blue')
axs[1].plot([ls[-1]], [plot_datas[6][-1]], marker='o', color='blue', label='_nolegend_')

axs[1].plot(ls, plot_datas[7], label='CMB, stage 4, \n post-Born corr.', linestyle='-', color='blue')
axs[1].plot([ls[-1]], [plot_datas[7][-1]], marker='o', color='blue', label='_nolegend_')


axs[0].minorticks_on()
axs[1].minorticks_on()

axs[0].set_xlabel(r'$l_{\max}$', fontsize=10)
axs[0].set_ylabel(r'$S/N$', fontsize=10)
axs[1].set_xlabel(r'$l_{\max}$', fontsize=10)
# axs[1].set_ylabel(r'$S/N$', fontsize=12)

# Increase the x-axis limit to add extra space on the right
axs[0].set_ylim(0, np.sqrt(snrb_s4_s4[-1]) * 1.125)
axs[1].set_ylim(0, np.sqrt(snrb_c_s4[-1]) * 1.125)

x_max = ell_max[-1]
axs[0].set_xlim(2, x_max + 250)
axs[1].set_xlim(2, x_max + 250)

# Annotate each line with the y-value at the right end
offset = 50  # horizontal offset for the text
axs[1].text(x_max + offset, np.sqrt(snrb_c_s3)[-1], f'{np.sqrt(snrb_c_s3)[-1]:.0f}', color='black', ha='left', va='center')
axs[0].text(x_max + offset, np.sqrt(snrb_s1_s3)[-1], f'{np.sqrt(snrb_s1_s3)[-1]:.0f}', color='black', ha='left', va='center')
axs[0].text(x_max + offset, np.sqrt(snrb_s4_s3)[-1], f'{np.sqrt(snrb_s4_s3)[-1]:.0f}', color='blue', ha='left', va='center')
axs[1].text(x_max + offset, np.sqrt(snrb_c_s4)[-1], f'{np.sqrt(snrb_c_s4)[-1]:.0f}', color='black', ha='left', va='center')
axs[0].text(x_max + offset, np.sqrt(snrb_s1_s4)[-1], f'{np.sqrt(snrb_s1_s4)[-1]:.0f}', color='black', ha='left', va='center')
axs[0].text(x_max + offset, np.sqrt(snrb_s4_s4)[-1], f'{np.sqrt(snrb_s4_s4)[-1]:.0f}', color='blue', ha='left', va='center')

axs[1].text(x_max + offset, np.sqrt(snrb_c_s3_pb)[-1], f'{np.sqrt(snrb_c_s3_pb)[-1]:.0f}', color='black', ha='left', va='center')
axs[1].text(x_max + offset, np.sqrt(snrb_c_s4_pb)[-1], f'{np.sqrt(snrb_c_s4_pb)[-1]:.0f}', color='black', ha='left', va='center')

axs[0].legend(fontsize=10, loc='upper left')
axs[1].legend(fontsize=10, loc='upper left')

axs[0].grid(True, which='both', linestyle='--', linewidth=0.5)
axs[1].grid(True, which='both', linestyle='--', linewidth=0.5)

plt.tight_layout()

# Save the figure
plt.savefig('snrplots.pdf')
plt.savefig('snrplots.png')