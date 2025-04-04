import matplotlib.pyplot as plt
import numpy as np

# Define data
ell_max = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]

# snrp_c = [28177.879112088674, 96206.55433780684, 175853.35058235622, 243509.61866955983, 289830.6864329934, 317440.70691750693, 332633.4570946259, 340837.0890864015]

# snrp_s = [23441.8759542758, 67396.65825513331, 113536.84597591546, 156335.21499038063, 194043.68425517925, 226409.47923775774, 253826.42455016664, 276933.44511557865]

# new
# Convergence bispectrum SNR^2:
# [0.4606180302137393, 2.428277643136825, 4.556412710716516, 5.924272288700067, 6.620955446735307, 6.958473094521213, 7.12003774130501, 7.201717334076159]
# Shear bispectrum SNR^2:
# [14.937583099693954, 44.263581102561155, 70.4841150784228, 91.60588005582235, 108.3741775299598, 121.44133496518877, 131.79346525781503, 139.98497349255402]

snrb_c_s3 = [0.4606180302137393, 2.428277643136825, 4.556412710716516, 5.924272288700067, 6.620955446735307, 6.958473094521213, 7.12003774130501, 7.201717334076159]

snrb_s_s3 = [14.937583099693954, 44.263581102561155, 70.4841150784228, 91.60588005582235, 108.3741775299598, 121.44133496518877, 131.79346525781503, 139.98497349255402]

snrb_c_s4 = [0.7032144111436794, 6.542889778749563, 20.500697932075894, 39.87427410664932, 59.254915437097765, 74.95646527206883, 86.07133128729137, 93.3095234367019]

snrb_s_s4 = [50.95719570421603, 307.30117549383294, 748.571492858229, 1289.9230963915802, 1863.479395132727, 2419.0805380400534, 2928.1361348954356, 3387.6648308426447]

# Set figure size
plt.figure(figsize=(7,5))

# Plot each array
plt.plot(ell_max, np.sqrt(snrb_c_s3), marker='x', label='CMB, stage 3', linestyle = '--', color='red')
plt.plot(ell_max, np.sqrt(snrb_s_s3), marker='x', label='galaxy, stage 3', linestyle = '--', color='blue')
plt.plot(ell_max, np.sqrt(snrb_c_s4), marker='x', label='CMB, stage 4', color='red')
plt.plot(ell_max, np.sqrt(snrb_s_s4), marker='x', label='galaxy, stage 4', color='blue')

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
