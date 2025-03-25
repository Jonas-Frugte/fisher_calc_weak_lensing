import matplotlib.pyplot as plt
import numpy as np

# Define data
ell_max = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]

snrp_c = np.array([15989.588551558505, 31618.168645038517, 39928.86302292515, 
                     44034.90338600421, 46105.011942758254, 47230.495325004005, 
                     47860.80819355356, 48209.52024209044])

snrp_s = np.array([22253.427444706474, 66208.20974556397, 112348.39746634613, 
                     155146.7664808113, 192855.2357456099, 225221.0307281884, 
                     252637.9760405973, 275744.9966060093])

snrb_c = np.array([0.4662713058660203, 1.7371421603331405, 2.871500243034231, 
                     3.65500342260037, 4.15133536761963, 4.4659414552800625, 
                     4.664800647728199, 4.786812407063142])

snrb_s = np.array([34.868772005476124, 263.4320514442226, 685.031116048985, 
                     1215.6310184788954, 1779.353150346276, 2329.5943353466655, 
                     2840.639407847829, 3300.6119934544467])

# Set figure size
plt.figure(figsize=(7,5))

# Plot each array
plt.plot(ell_max, np.sqrt(snrp_c), marker='x', label='CMB, powerspectrum', linestyle='--', color='blue')
plt.plot(ell_max, np.sqrt(snrp_s), marker='x', label='galaxy, powerspectrum', color='blue')
plt.plot(ell_max, np.sqrt(snrb_c), marker='x', label='CMB, bispectrum', linestyle='--', color='red')
plt.plot(ell_max, np.sqrt(snrb_s), marker='x', label='galaxy, bispectrum', color='red')

plt.minorticks_on()
plt.xlabel(r'$l_{\max}$', fontsize=12)
plt.ylabel(r'$S/N$', fontsize=12)

# Increase the x-axis limit to add extra space on the right
x_max = ell_max[-1]
plt.xlim(ell_max[0] - 50, x_max + 250)

# Annotate each line with the y-value at the right end
offset = 50  # horizontal offset for the text
plt.text(x_max + offset, np.sqrt(snrp_c)[-1], f'{np.sqrt(snrp_c)[-1]:.2f}', color='black', ha='left', va='center')
plt.text(x_max + offset, np.sqrt(snrp_s)[-1], f'{np.sqrt(snrp_s)[-1]:.2f}', color='black', ha='left', va='center')
plt.text(x_max + offset, np.sqrt(snrb_c)[-1], f'{np.sqrt(snrb_c)[-1]:.2f}', color='black', ha='left', va='center')
plt.text(x_max + offset, np.sqrt(snrb_s)[-1], f'{np.sqrt(snrb_s)[-1]:.2f}', color='black', ha='left', va='center')

plt.legend(fontsize=10, loc='upper left')
plt.tight_layout()

# Save the figure
plt.savefig('/Users/jonasfrugte/Desktop/Research_Project/paper/figures/snrplots.pdf', dpi=300)
# plt.show()
