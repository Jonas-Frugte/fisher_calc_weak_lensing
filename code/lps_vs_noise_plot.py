import cosm_setup as cs
import matplotlib.pyplot as plt
import numpy as np

fig, axs = plt.subplots(1, 2, figsize=(8.3, 4))

# LENSING SPECTRA
l_values = np.logspace(np.log10(2), np.log10(3000), 100)
spectra = cs.lensing_spectra()
gal_lps = np.array([spectra.lps(l, ('s', 's')) for l in l_values])
cmb_lps = np.array([spectra.lps(l, ('c', 'c')) for l in l_values])

# CMB NOISE QUAD ESTIMATOR
ls_cmbn = np.loadtxt('cmb_noise_files/ls_1_3000_64.txt')
cmbn_301 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma3_DeltaT0.71_DeltaP1.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_106 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma1_DeltaT4.2_DeltaP6.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_TT = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_TT.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_TE = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_TE.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_TB = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_TB.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_EE = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_EE.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_EB = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_EB.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_test_vals = [cmbn_TT, cmbn_TE, cmbn_TB, cmbn_EE, cmbn_EB]
names = ['TT', 'TE', 'TB', 'EE', 'EB']

# SIMONS OBSERVATORY NOISE LEVELS
conv_noise_goal_data_array = np.loadtxt('./cmb_noise_files/nlkk_v3_1_0_deproj0_SENS2_fsky0p4_it_lT30-3000_lP30-5000.dat')
SO_noise_goal = lambda l: 0.4 * conv_noise_goal_data_array[l-2, 7] * 4.
SO_noise_goal_values = np.array([SO_noise_goal(int(l)) for l in l_values])

# GALAXY LENSING NOISE
sigma_rms_S3 = 0.3
sigma_rms_S4 = 0.3
arcmin2_to_steradian = (np.pi / (180 * 60)) ** 2
n_g_S3 = 5 / arcmin2_to_steradian
n_g_S4 = 30 / arcmin2_to_steradian
conversion_factor = 4 / ((l_values - 1) * l_values * (l_values + 1) * (l_values + 2))
gal_nps_S3 = np.array([sigma_rms_S3 ** 2 / n_g_S3 for l in l_values]) * conversion_factor
gal_nps_S4 = np.array([sigma_rms_S4 ** 2 / n_g_S4 for l in l_values]) * conversion_factor

# ------- Plot 1: Galaxy lensing -------
axs[0].loglog(l_values, l_values**4 * gal_lps, label='Galaxy lens. pot. powerspec.', color='black')
axs[0].loglog(l_values, l_values**4 * gal_nps_S3, label=r'S3 Noise, $n_g = 5\,\mathrm{arcmin}^{-2}$', linestyle='dashed')
axs[0].loglog(l_values, l_values**4 * gal_nps_S4, label=r'S4 Noise, $n_g = 30\,\mathrm{arcmin}^{-2}$', linestyle='dashed')
axs[0].set_xlabel('$l$')
axs[0].set_ylabel(r'$l^4 C_l^{\psi_g\psi_g}$')
axs[0].set_xlim(2, 2000)
axs[0].grid(True, which='both', linestyle='--', linewidth=0.5)
axs[0].legend(
    loc='upper center',
    bbox_to_anchor=(0.5, -0.15),
    ncol=3,
    fontsize='small'
)

# ------- Plot 2: CMB lensing -------
axs[1].loglog(l_values, l_values**4 * cmb_lps, label='CMB lens. pot. powerspec.', color='black')
axs[1].loglog(l_values, SO_noise_goal_values, label='SO Noise, goal', linestyle='dashed', color='green')
axs[1].loglog(ls_cmbn, ls_cmbn**4 * cmbn_106, label=r'S3 noise, $\sigma=1,\ \Delta_T=4.2,\ \Delta_P=6$', linestyle='dashed')
axs[1].loglog(ls_cmbn, ls_cmbn**4 * cmbn_301, label=r'S4 noise, $\sigma=3,\ \Delta_T=0.71,\ \Delta_P=1$', linestyle='dashed')
axs[1].set_xlabel('$l$')
axs[1].set_ylabel(r'$l^4 C_l^{\psi_c\psi_c}$')
axs[1].set_xlim(2, 2000)
axs[1].grid(True, which='both', linestyle='--', linewidth=0.5)
axs[1].legend(
    loc='upper center',
    bbox_to_anchor=(0.5, -0.15),
    ncol=2,
    fontsize='small'
)

# space for legends
fig.subplots_adjust(bottom=0.25, wspace=0.3)

# save
plt.savefig("plots/spectraplusnoise.png", format="png", dpi=300)
plt.savefig("plots/spectraplusnoise.pdf", format="pdf", dpi=300)
