import cosm_setup as cs
import matplotlib.pyplot as plt
import numpy as np

fig, axs = plt.subplots(1, 2, figsize = (14, 6))

# original noise is of deflection signal, so divide by l(l+1) to convert to lens potential


ls_cmbn = np.loadtxt('cmb_noise_files/ls_1_3000_64.txt')
cmbn_41sqrt2 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_301 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma3_DeltaT0_DeltaP1.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_106 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma1_DeltaT0_DeltaP6.txt')) / (ls_cmbn * (ls_cmbn + 1))
# cmbn_410 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP0.txt')) / (ls_cmbn * (ls_cmbn + 1))

# !!!!!! still with lmin = 10 !!!!!!!!
cmbn_TT = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_TT.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_TE = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_TE.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_TB = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_TB.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_EE = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_EE.txt')) / (ls_cmbn * (ls_cmbn + 1))
cmbn_EB = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma4_DeltaT1_DeltaP1.41_EB.txt')) / (ls_cmbn * (ls_cmbn + 1))

cmbn_test_vals = [cmbn_TT, cmbn_TE, cmbn_TB, cmbn_EE, cmbn_EB]
names = ['TT', 'TE', 'TB', 'EE', 'EB']

# Compute lensing power spectrum
l_values = np.logspace(np.log10(2), np.log10(3000), 100)
spectra = cs.lensing_spectra(fiducial_k_nls=True) # don't need matter bispectrum
cl_ss = np.array([spectra.lps(l, ('s', 's')) for l in l_values])  # Shear-shear power spectrum

folder_file_path = '/scratch/p319950/data/'
filepath_convergence_noise_file_path = folder_file_path + 'conv_noise.dat'

conv_noise_data_array = np.loadtxt(filepath_convergence_noise_file_path)
conv_pot_noise = lambda l : conv_noise_data_array[l-2, 7] * 4. # * (l * 1.0)**(-2) * (l + 1.0)**(-2)
cl_gg = np.array([spectra.lps(l, ('c', 'c')) for l in l_values])
conv_pot_noise_values = np.array([conv_pot_noise(int(l)) for l in l_values])

filepath_convergence_noise_goal_file_path = folder_file_path + 'conv_noise_goal.dat'

conv_noise_goal_data_array = np.loadtxt(filepath_convergence_noise_goal_file_path)
conv_pot_noise_goal = lambda l : conv_noise_goal_data_array[l-2, 7] * 4. # * (l * 1.0)**(-2) * (l + 1.0)**(-2)
conv_pot_noise_goal_values = np.array([conv_pot_noise_goal(int(l)) for l in l_values])

# Constants for noise
sigma_rms_LSST = 0.26  # LSST ellipticity dispersion
sigma_rms_Euclid = 0.3  # Euclid ellipticity dispersion

# Number density per steradian (converted from per arcmin^2)
arcmin2_to_steradian = (np.pi / (180 * 60)) ** 2
n_g_LSST = 26 # / arcmin2_to_steradian
n_g_Euclid = 30 # / arcmin2_to_steradian

# Compute shear noise power spectra
# do I need the factor of 4 here??
N_l_shear_LSST = sigma_rms_LSST ** 2 / (n_g_LSST)
N_l_shear_Euclid = sigma_rms_Euclid ** 2 / (n_g_Euclid)

# Compute lensing potential noise properly
N_l_phi_LSST = np.array([N_l_shear_LSST for l in l_values]) # * (l_values * (l_values + 1))**2
N_l_phi_Euclid = np.array([N_l_shear_Euclid for l in l_values]) # * (l_values * (l_values + 1))**2

# Convert to lensing potential using the right factor
conversion_factor = 4 / ((l_values - 1) * l_values * (l_values + 1) * (l_values + 2))
N_l_phi_LSST *= conversion_factor
N_l_phi_Euclid *= conversion_factor

# Plot
axs[0].loglog(l_values, l_values**2 * (l_values + 1)**2 * cl_ss, label='Galaxy lensing potential powerspectrum', color='black')
axs[0].loglog(l_values, l_values**2 * (l_values + 1)**2 * N_l_phi_LSST, label='LSST Noise', linestyle='dashed', color='purple')
axs[0].loglog(l_values, l_values**2 * (l_values + 1)**2 * N_l_phi_Euclid, label='Euclid Noise', linestyle='dashed', color='red')
axs[0].set_xlabel('$l$')
axs[0].set_ylabel('$l^4 C_l$')
axs[0].legend()
axs[0].grid(True, which='both', linestyle='--', linewidth=0.5)

axs[1].loglog(l_values, l_values**2 * (l_values + 1)**2 * cl_gg, label='CMB lensing potential powerspectrum', color='black')
axs[1].loglog(l_values, conv_pot_noise_values, label='SO Noise, baseline', linestyle='dashed', color='blue')
axs[1].loglog(l_values, conv_pot_noise_goal_values, label='SO Noise, goal', linestyle='dashed', color='blue')

#axs[1].loglog(ls_cmbn, ls_cmbn**2 * (ls_cmbn + 1)**2 * cmbn_41sqrt2, label=r'quad est. noise, $\sigma = 4$, $\Delta_T = 1$, $\Delta_P = \sqrt{2}$', linestyle='dashed', color='red')

# for i in range(len(cmbn_test_vals)):
#     axs[1].loglog(ls_cmbn, ls_cmbn**2 * (ls_cmbn + 1)**2 * cmbn_test_vals[i], linestyle='dotted', label = names[i])
# axs[1].loglog(ls_cmbn, ls_cmbn**2 * (ls_cmbn + 1)**2 * cmbn_410, label=r'quad est. noise, $\sigma = 4$, $\Delta_T = 1$, $\Delta_P = \infty$', linestyle='dotted', color='red')
axs[1].loglog(ls_cmbn, ls_cmbn**2 * (ls_cmbn + 1)**2 * cmbn_301, label=r'quad est. noise, $\sigma = 3$, $\Delta_T = \infty$, $\Delta_P = 1$ (S4 Noise Toshiya (2016))', linestyle='dashed', color='darkviolet')
axs[1].loglog(ls_cmbn, ls_cmbn**2 * (ls_cmbn + 1)**2 * cmbn_106, label=r'quad est. noise, $\sigma = 1$, $\Delta_T = \infty$, $\Delta_P = 6$ (S3-wide Noise Toshiya (2016))', linestyle='dashed', color='green')

axs[1].set_xlabel('$l$')
axs[1].set_ylabel('$l^2 (l+1)^2 C_l$')
axs[1].legend()
axs[1].grid(True, which='both', linestyle='--', linewidth=0.5)
axs[1].set_xlim(2, 2000)
fig.suptitle('CMB and galaxy lensing potential powerspectra with noise')
fig.tight_layout()

# Save the figure as a vector-based format for inclusion in papers
#plt.savefig("spectraplusnoise.pdf", format="pdf", dpi=300)
plt.savefig("spectraplusnoise.png", format="png", dpi=300)