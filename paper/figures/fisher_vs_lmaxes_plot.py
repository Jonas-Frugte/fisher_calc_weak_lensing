import numpy as np
import matplotlib.pyplot as plt

SNR2_c = np.array([21752.37109375, 42467.37109375, 61431.392578125, 80951.72265625, 103566.880859375, 127987.7421875])

SNR2_s = np.array([58.471466064453125, 209.129638671875, 472.9197998046875, 857.3757019042969, 1365.5498962402344, 1985.4261779785156])

SNR2_full = np.array([15973925.273764089, 21700691.52666128, 25006818.371009313, 27387274.95562951, 29266159.18653201, 31273214.900456935])

# Given Data
l_values = np.array([500, 1000, 1500, 2000, 2500, 3000])

Fisher_c = {
    "H0": np.array([14.077168464660645, 31.897847175598145, 50.746644020080566, 72.83032321929932, 99.00564861297607, 127.41051578521729]),
    "ombh2": np.array([4112962.0, 9336375.0, 20095801.0, 26203080.5, 38208011.5, 53317144.5]),
    "omch2": np.array([13600145.0, 32862975.0, 52410357.0, 77972975.0, 110291517.0, 147623545.0]),
    "As": np.array([1.9673580907200617e+22, 4.143315935600502e+22, 6.3830039936278815e+22, 8.84181108864527e+22, 1.1733016440005867e+23, 1.5001670052396254e+23]),
    "ns": np.array([12025.2470703125, 28862.3662109375, 55617.6474609375, 118703.4716796875, 219686.6748046875, 355107.7529296875]),
    "mnu": np.array([20274.103515625, 44399.357421875, 67587.205078125, 92397.80859375, 121286.11328125, 154003.509765625]),
}

Fisher_s = {
    "H0": np.array([0.27174660563468933, 0.980530709028244, 2.1956381499767303, 3.912103980779648, 6.173850864171982, 8.896407455205917]),
    "ombh2": np.array([43218.79296875, 174591.40234375, 407098.21484375, 743700.71484375, 1195916.71484375, 1761281.83984375]),
    "omch2": np.array([62721.1875, 276152.34375, 634106.46875, 1145792.4375, 1821605.5625, 2635287.3125]),
    "As": np.array([6.681161735362262e+19, 2.6235997145981675e+20, 5.815499914061661e+20, 1.0208579799201601e+21, 1.59092147450602e+21, 2.2675664190014164e+21]),
    "ns": np.array([95.61677551269531, 763.2185821533203, 2276.7918243408203, 4943.978103637695, 8984.006912231445, 14529.572830200195]),
    "mnu": np.array([85.50411987304688, 320.78245544433594, 708.7311248779297, 1249.472396850586, 1961.3191375732422, 2815.858139038086]),
}

Fisher_full = {
    "H0": np.array([4266.01521830935, 7111.420249883808, 8937.133902515721, 10232.61698307895, 11299.719275861975, 12436.480527096597]),
    "ombh2": np.array([2077065333.8645613, 3529874774.623591, 8230689040.315075, 12454985671.09074, 18136645958.554726, 20693198967.705315]),
    "omch2": np.array([5209024703.990475, 8666390335.668026, 10952289312.67326, 13034649615.6269, 15058592704.440039, 17100423717.402327]),
    "As": np.array([1.5315134367951586e+25, 2.1273318211413612e+25, 2.4815247880764978e+25, 2.745803747189947e+25, 2.9652898322314457e+25, 3.190373220272011e+25]),
    "ns": np.array([15060116.363071065, 17496111.440253276, 19897615.562558837, 23380381.22884735, 28916259.137280077, 34575529.96079264]),
    "mnu": np.array([14300425.868088072, 19508700.507455535, 22862434.155531663, 25690325.463063225, 27770792.768567003, 29864098.513198208]),
}

fiducial_values = {
    "H0": 67.4,
    "ombh2": 0.0224,
    "omch2": 0.120,
    "As": 2.1e-9,
    "ns": 0.965,
    "mnu": 0.06,
}

updated_titles = [
    r"H_0",
    r"\Omega_b h^2",
    r"\Omega_c h^2",
    r"A_s",
    r"n_s",
    r"m_{\nu}",
]

# Generate the plot
height = 9
fig, axs = plt.subplots(2, 3, figsize=(1.5 * height, height))

parameters = list(Fisher_c.keys())

for idx, param in enumerate(parameters):
    ax = axs[idx // 3, idx % 3]
    param_fiducial = fiducial_values[param]
    c_values = 1 / (np.sqrt(Fisher_c[param]) * param_fiducial)
    s_values = 1 / (np.sqrt(Fisher_s[param]) * param_fiducial)
    full_values = 1 / (np.sqrt(Fisher_full[param]) * param_fiducial)
    
    ax.semilogy(l_values, c_values, label='convergence only', marker = '+', color = 'forestgreen')
    ax.semilogy(l_values, s_values, label='shear only', marker = '+', color = 'darkcyan')
    ax.semilogy(l_values, full_values, label='all lensing bispectra', marker = '+', color = 'royalblue')
    
    ax.set_title(rf"${updated_titles[idx]}$", fontsize=14)
    ax.set_xlabel(r"$l_{max}$", fontsize=12)
    y_label = rf"$(\sqrt{{F_{{{updated_titles[idx]},{updated_titles[idx]}}}}} {updated_titles[idx]})^{{-1}}$"
    ax.set_ylabel(y_label, fontsize=12)
    ax.legend()
    ax.minorticks_on()
    ax.grid(True, which='major', linestyle='-', linewidth=0.8)  # Major grid
    ax.grid(True, which='minor', linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.savefig('/Users/jonasfrugte/Desktop/ResearchProject/Paper/figures/fisher_results_theoretical.png', dpi = 300)


SNR_labels = [
    'convergence only',
    'shear only',
    'all lensing bispectra'
]

SNR_data = [
    SNR2_c,
    SNR2_s,
    SNR2_full
]

fig, ax = plt.subplots(figsize = (6 * 1.5, 5))
for i in range(len(SNR_labels)):
    ax.semilogy(l_values, 1/np.sqrt(SNR_data[i]), label = SNR_labels[i], marker = '+')

ax.set_xlabel(r"$l_{max}$", fontsize=12)
ax.set_ylabel("relative uncertainty of bispectra amplitude (SNR)", fontsize = 12)
ax.legend()
plt.tight_layout()
ax.minorticks_on()
ax.grid(True, which='major', linestyle='-', linewidth=0.8)  # Major grid
ax.grid(True, which='minor', linestyle='--', linewidth=0.5)
plt.savefig("/Users/jonasfrugte/Desktop/ResearchProject/Paper/figures/SNR_results_theoretical.png", dpi = 300)

