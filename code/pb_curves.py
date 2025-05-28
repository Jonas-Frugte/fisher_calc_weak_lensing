import numpy as np
import matplotlib.pyplot as plt
import os

import cosm_setup as cs
spectra = cs.lensing_spectra()

ls = np.arange(100, 2000, 2)

lbs_eq_ex_tree = np.array([
    spectra.lbs_flat(l, l, l, ('c', 'c', 'c'), model = 'tree') * float(l)**6 / 8 for l in ls
])

lbs_fold_ex_tree = np.array([
    spectra.lbs_flat(l, l/2, l/2, ('c', 'c', 'c'), model = 'tree') * float(l)**6 / 2**4 / 8 for l in ls
])

# Define the directory and file information (filename: label)
data_directory = 'pb_curves/'
file_info = {
    'lss_lbs_eq_nonl.txt': 'LSS Lensing Bispectrum Nonlinear',
    'lss_lbs_eq_tree.txt': 'LSS Lensing Bispectrum Tree-Level',
    'lss_lbs_fold_nonl.txt': 'LSS Lensing Bispectrum Nonlinear',
    'lss_lbs_fold_tree.txt': 'LSS Lensing Bispectrum Tree-Level',
    'pb_lbs_eq.txt': 'Post-Born Lensing Bispectrum',
    'pb_lbs_fold.txt': 'Post-Born Lensing Bispectrum'
}

# Create the figure and subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# --- Plotting for 'eq' files ---
ax1.set_title('Equilateral Lensing Bispectra')
ax1.set_xlabel('$\ell$') # Assuming x-axis is multipole l
ax1.set_ylabel('Bispectrum $B(\ell, \ell, \ell)$') # Adjust ylabel as needed

for filename, label in file_info.items():
    if 'eq' in filename:
        path = os.path.join(data_directory, filename)
        # Data has 2 columns: 0 for x-values, 1 for y-values
        # print(path)
        data = np.loadtxt(path)
        # print(data)
        ax1.loglog(data[:, 0], data[:, 1], label=label)



# --- Plotting for 'fold' files ---
ax2.set_title('Folded Lensing Bispectra')
ax2.set_xlabel('$\ell$') # Assuming x-axis is multipole l
ax2.set_ylabel('Bispectrum $B(\ell, \ell, \ell)$') # Adjust ylabel as needed

for filename, label in file_info.items():
    if 'fold' in filename:
        path = os.path.join(data_directory, filename)
        # Data has 2 columns: 0 for x-values, 1 for y-values
        data = np.loadtxt(path)
        ax2.loglog(data[:, 0], data[:, 1], label=label)


import data_importer_new as di

# lbs_pb_flat_f(double l1, double l2, double l3, int num_samples)

pb_lbs_eq_vals = [di.lbs_pb_flat_f(l, l, l, 200) for l in ls]
pb_lbs_fold_vals = [np.abs(di.lbs_pb_flat_f(l, l/2, l/2, 200)) for l in ls]
lbs_eq_tree = [float(l)**6 / 8 * di.lbs_flat_f(l, l, l, b'c', b'c', b'c', 200) for l in ls]
lbs_fold_tree = [float(l)**6 / 8 / 2**4 * di.lbs_flat_f(l, l/2, l/2, b'c', b'c', b'c', 200) for l in ls]



ax1.loglog(ls, pb_lbs_eq_vals, linestyle = '--', label = 'pb')
ax1.loglog(ls, lbs_eq_tree, linestyle = '--', label = 'nonl lbs')
ax2.loglog(ls, pb_lbs_fold_vals, linestyle = '--', label = 'pb')
ax2.loglog(ls, lbs_fold_tree, linestyle = '--', label = 'nonl lbs')

ax1.loglog(ls, lbs_eq_ex_tree, linestyle = 'dotted', color = 'black', label = 'lin lbs')
ax2.loglog(ls, lbs_fold_ex_tree, linestyle = 'dotted', color = 'black', label = 'lin lbs')

ax1.legend()
ax1.grid(True, which="both", ls="--", alpha=0.5)

ax2.legend()
ax2.grid(True, which="both", ls="--", alpha=0.5)


# Adjust layout to prevent overlapping titles/labels
plt.tight_layout()

# Show the plot
plt.savefig('pb_curves.png')