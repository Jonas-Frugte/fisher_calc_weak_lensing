import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import norm
from numpy.linalg import inv

np.set_printoptions(
    precision=1,
    suppress=False, 
    formatter={'float_kind': lambda x: f"{x:.2e}"}
)

######################
# Bispectra Matrices #
######################
visb_c = np.array([
    [ 3.66134328e-03,  5.64294103e-04, -5.59024163e-03, -8.28197198e-03,  9.21839119e-05, -3.43845273e-03, -1.21009701e-03],
    [ 5.64294103e-04,  9.39182080e-05, -8.78913206e-04, -1.41245463e-03,  1.42075865e-05, -5.24072512e-04, -1.83365409e-04],
    [-5.59024163e-03, -8.78913206e-04,  8.58128614e-03,  1.29940930e-02, -1.40723722e-04,  5.23439712e-03,  1.84009658e-03],
    [-8.28197198e-03, -1.41245463e-03,  1.29940930e-02,  2.15235269e-02, -2.08249724e-04,  7.65644963e-03,  2.67469517e-03],
    [ 9.21839119e-05,  1.42075865e-05, -1.40723722e-04, -2.08249724e-04,  2.32194050e-06, -8.65918443e-05, -3.04632746e-05],
    [-3.43845273e-03, -5.24072512e-04,  5.23439712e-03,  7.65644963e-03, -8.65918443e-05,  3.23545639e-03,  1.13828815e-03],
    [-1.21009701e-03, -1.83365409e-04,  1.84009658e-03,  2.67469517e-03, -3.04632746e-05,  1.13828815e-03,  4.01520985e-04]
])

visb_s = np.array([
    [ 1.43008841e+01,  1.58650740e+00, -1.14635424e+01, -8.04509013e+00,  1.95515860e-01, -5.98457532e+00, -5.21205911e+00],
    [ 1.58650740e+00,  1.78726181e-01, -1.28311279e+00, -9.14981146e-01,  2.18306721e-02, -6.68745848e-01, -5.78208226e-01],
    [-1.14635424e+01, -1.28311279e+00,  9.23957956e+00,  6.54726869e+00, -1.57349525e-01,  4.81985548e+00,  4.17900899e+00],
    [-8.04509013e+00, -9.14981146e-01,  6.54726869e+00,  4.73292722e+00, -1.11131095e-01,  3.40652074e+00,  2.93071743e+00],
    [ 1.95515860e-01,  2.18306721e-02, -1.57349525e-01, -1.11131095e-01,  2.68124786e-03, -8.20991755e-02, -7.12976000e-02],
    [-5.98457532e+00, -6.68745848e-01,  4.81985548e+00,  3.40652074e+00, -8.20991755e-02,  2.51368749e+00,  2.18218608e+00],
    [-5.21205911e+00, -5.78208226e-01,  4.17900899e+00,  2.93071743e+00, -7.12976000e-02,  2.18218608e+00,  1.89999791e+00]
])

visb_f = np.array([
    [ 1.44863736e+01,  1.60644872e+00, -1.15924280e+01, -8.14549164e+00,  1.97919543e-01, -6.05641501e+00, -5.28028125e+00],
    [ 1.60644872e+00,  1.80999610e-01, -1.29933301e+00, -9.27103116e-01,  2.20864292e-02, -6.77304682e-01, -5.85552679e-01],
    [-1.15924280e+01, -1.29933301e+00,  9.37147269e+00,  6.64810299e+00, -1.59417592e-01,  4.88501137e+00,  4.22882478e+00],
    [-8.14549164e+00, -9.27103116e-01,  6.64810299e+00,  4.81725236e+00, -1.12756306e-01,  3.45659641e+00,  2.96580280e+00],
    [ 1.97919543e-01,  2.20864292e-02, -1.59417592e-01, -1.12756306e-01,  2.71512907e-03, -8.31587413e-02, -7.21540131e-02],
    [-6.05641501e+00, -6.77304682e-01,  4.88501137e+00,  3.45659641e+00, -8.31587413e-02,  2.54832853e+00,  2.20862447e+00],
    [-5.28028125e+00, -5.85552679e-01,  4.22882478e+00,  2.96580280e+00, -7.21540131e-02,  2.20862447e+00,  1.92565489e+00]
])

########################
# Powerspectra Matrices
########################
visp_c = np.array([
    [7.21859499e+00,  1.83473419e+00, -1.45468172e+01, -1.93821124e+00,  2.77227313e-01, -9.02220580e+00, -3.51298653e+00],
    [1.83473419e+00,  5.11468909e-01, -3.82290639e+00, -8.81690525e-01,  7.04920642e-02, -2.24589063e+00, -8.84381549e-01],
    [-1.45468172e+01, -3.82290639e+00,  2.96664255e+01,  5.01581744e+00, -5.58578281e-01,  1.80461349e+01,  7.05560491e+00],
    [-1.93821124e+00, -8.81690525e-01,  5.01581744e+00,  4.25336634e+00, -7.20688928e-02,  1.95900114e+00,  8.69447068e-01],
    [2.77227313e-01,  7.04920642e-02, -5.58578281e-01, -7.20688928e-02,  1.06669982e-02, -3.46858189e-01, -1.34905800e-01],
    [-9.02220580e+00, -2.24589063e+00,  1.80461349e+01,  1.95900114e+00, -3.46858189e-01,  1.13345007e+01,  4.39975793e+00],
    [-3.51298653e+00, -8.84381549e-01,  7.05560491e+00,  8.69447068e-01, -1.34905800e-01,  4.39975793e+00,  1.71129547e+00]
])

visp_s = np.array([
    [5.78350209e+02,  6.84835881e+01, -4.71217702e+02, -2.52313437e+02,  8.36536571e+00, -2.49174166e+02, -2.46662642e+02],
    [6.84835881e+01,  8.24006706e+00, -5.63418508e+01, -3.09043701e+01,  9.97208461e-01, -2.97548996e+01, -2.92084350e+01],
    [-4.71217702e+02, -5.63418508e+01,  3.86240859e+02,  2.09778080e+02, -6.84380558e+00,  2.04076045e+02,  2.00984091e+02],
    [-2.52313437e+02, -3.09043701e+01,  2.09778080e+02,  1.18355893e+02, -3.70083229e+00,  1.10603502e+02,  1.07579402e+02],
    [8.36536571e+00,  9.97208461e-01, -6.84380558e+00, -3.70083229e+00,  1.21341173e-01, -3.61711370e+00, -3.56799205e+00],
    [-2.49174166e+02, -2.97548996e+01,  2.04076045e+02,  1.10603502e+02, -3.61711370e+00,  1.07855534e+02,  1.06284242e+02],
    [-2.46662642e+02, -2.92084350e+01,  2.00984091e+02,  1.07579402e+02, -3.56799205e+00,  1.06284242e+02,  1.05211495e+02]
])

visp_f = np.array([
    [5.88405126e+02,  7.05983979e+01, -4.88362694e+02, -2.57408810e+02,  8.66423625e+00, -2.58762845e+02, -2.51316256e+02],
    [7.05983979e+01,  8.79146209e+00, -6.05110535e+01, -3.21085595e+01,  1.07162945e+00, -3.21251937e+01, -3.02111649e+01],
    [-4.88362694e+02, -6.05110535e+01,  4.19640517e+02,  2.17989008e+02, -7.44094073e+00,  2.23344637e+02,  2.09126551e+02],
    [-2.57408810e+02, -3.21085595e+01,  2.17989008e+02,  1.24804734e+02, -3.81280992e+00,  1.13695789e+02,  1.09756543e+02],
    [8.66423625e+00,  1.07162945e+00, -7.44094073e+00, -3.81280992e+00,  1.32376315e-01, -3.97703005e+00, -3.71222356e+00],
    [-2.58762845e+02, -3.21251937e+01,  2.23344637e+02,  1.13695789e+02, -3.97703005e+00,  1.19679129e+02,  1.10940553e+02],
    [-2.51316256e+02, -3.02111649e+01,  2.09126551e+02,  1.09756543e+02, -3.71222356e+00,  1.10940553e+02,  1.07385733e+02]
])

########### NORMALIZATION ###########

param_values_normalization = np.array([
    67.4,         # H_0
    0.0224,       # Omega_b h^2
    0.120,        # Omega_c h^2
    0.965,        # n_s
    0.06,         # m_\nu
    2.1e-9,       # A_s
    -1.0          # w_0
])

normalization_matrix = np.outer(param_values_normalization, param_values_normalization)
# If needed, multiply each matrix by this normalization_matrix.

########### SIGMA8 ################

def append_row_column(mats, vec):
    """
    For each matrix in 'mats':
      - Enlarge the matrix by 1 row + 1 column
      - New column entries: old_row(mat) · vec
      - New row same as new column
      - New bottom-right corner = sum( mat * outer(vec, vec) )
      - Returns an updated list (same references replaced).
    """
    matsnew = mats
    for i in range(len(mats)):
        mat = mats[i]
        N = mat.shape[0]
        
        # Compute new col/row
        new_col = mat @ vec  # shape (N,)
        new_row = new_col.copy()
        
        # Diagonal entry
        diag_entry = np.sum(np.outer(vec, vec) * mat)
        
        # Create an enlarged (N+1)x(N+1) array
        enlarged = np.zeros((N+1, N+1), dtype=mat.dtype)
        
        # Copy old block
        enlarged[:N, :N] = mat
        
        # Fill new column and row
        enlarged[:N, N] = new_col
        enlarged[N, :N] = new_row
        
        # Fill bottom-right entry
        enlarged[N, N] = diag_entry
        
        # Reassign the enlarged matrix back to the list
        matsnew[i] = enlarged
    
    return matsnew

s8_ders = np.array([2.76945054e-03, -6.50570238e+00,  4.30234556e+00,  3.06272435e-01,
                    1.93167857e+08, -2.15681667e-01, -2.01175433e-01])

visb_c, visb_s, visb_f, visp_c, visp_s, visp_f = append_row_column(
    [visb_c, visb_s, visb_f, visp_c, visp_s, visp_f],
    s8_ders
)

y = (0.0224 + 0.120) / (0.3 * 0.674**2) # used to calculate S8
S8_ders = [0.811 * y**(-1/2) * 0.5 * y * (0.674)**(-1) * (-2) * 0.01, 
           0.811 * y**(-1/2) * (0.3 * 0.674**2)**(-1) * 0.5,
           0.811 * y**(-1/2) * (0.3 * 0.674**2)**(-1) * 0.5,
           0,
           0,
           0,
           0,
           y**(1/2)]

visb_c, visb_s, visb_f, visp_c, visp_s, visp_f = append_row_column(
    [visb_c, visb_s, visb_f, visp_c, visp_s, visp_f],
    S8_ders
)

omm_ders = [0.811 * y**(-1/2) * 0.5 * y * (0.674)**(-1) * (-2) * 0.01, 
           1/0.674**2,
           1/0.674**2,
           0,
           0,
           0,
           0,
           0,
           0]

visb_c, visb_s, visb_f, visp_c, visp_s, visp_f = append_row_column(
    [visb_c, visb_s, visb_f, visp_c, visp_s, visp_f],
    omm_ders
)
######### PRINT EIGENVALUES ############

def print_sorted_eigens(matrix_name, mat):
    """Compute eigen-decomposition of a matrix, sort by descending eigenvalues, and print."""
    vals, vecs = np.linalg.eig(mat)
    idx = np.argsort(vals)[::-1]
    vals = vals[idx]
    vecs = vecs[:, idx]

    print(f"Eigenvalues ({matrix_name}), sorted:")
    print(vals)
    print(f"Eigenvectors ({matrix_name}) in columns, matching eigenvalues above:")
    for i, val in enumerate(vals):
        print(f"Eigenvector for eigenvalue {val:.6e}:")
        print(vecs[:, i])

# Indices to remove: 6,1,4,5 (0-based). This means we keep [0,2,3].
keep_indices = [0, 7, 9]  # We'll keep H_0 and sigma_8

# Reduced copies for bispectra
visb_c_reduced = visb_c[np.ix_(keep_indices, keep_indices)]
visb_s_reduced = visb_s[np.ix_(keep_indices, keep_indices)]
visb_f_reduced = visb_f[np.ix_(keep_indices, keep_indices)]

# Reduced copies for powerspectra
visp_c_reduced = visp_c[np.ix_(keep_indices, keep_indices)]
visp_s_reduced = visp_s[np.ix_(keep_indices, keep_indices)]
visp_f_reduced = visp_f[np.ix_(keep_indices, keep_indices)]

######### PLOTS #########

# def confidence_ellipse(ax, mean, cov, color='blue', nstd=1.0, **kwargs):
#     """
#     Plot the 2D confidence ellipse of a Gaussian distribution given by mean & covariance
#     onto the Matplotlib axes 'ax'.
#     - mean: (2,) array of the parameter means
#     - cov: (2x2) covariance matrix
#     - color: color for the ellipse
#     - nstd: number of standard deviations (1 => ~68% region, 2 => ~95%, etc.)
#     """
#     vals, vecs = np.linalg.eigh(cov)
#     # Angle to rotate the ellipse
#     angle = np.degrees(np.arctan2(*vecs[:, 1][::-1]))
#     # Major/minor axes
#     width, height = 2 * nstd * np.sqrt(vals)

#     ellip = Ellipse(
#         xy=mean, width=width, height=height,
#         angle=angle, edgecolor=color, facecolor='none', lw=2, **kwargs
#     )
#     ax.add_patch(ellip)

# new version that returns flipped version which I think is correct
def confidence_ellipse(ax, mean, cov, color='blue', nstd=1.0, **kwargs):
    """
    Plot the 2D confidence ellipse of a Gaussian distribution given by mean & covariance
    onto the Matplotlib axes 'ax'.
    - mean: (2,) array of the parameter means
    - cov: (2x2) covariance matrix
    - color: color for the ellipse
    - nstd: number of standard deviations (1 => ~68% region, 2 => ~95%, etc.)
    """
    vals, vecs = np.linalg.eigh(cov)
    # Angle to rotate the ellipse
    angle = 45 - (np.degrees(np.arctan2(*vecs[:, 1][::-1])) - 45)
    # Major/minor axes
    width, height = 2 * nstd * np.sqrt(vals)

    ellip = Ellipse(
        xy=mean, width=width, height=height,
        angle=angle, edgecolor=color, facecolor='none', lw=2, **kwargs
    )
    ax.add_patch(ellip)

def plot_corner(
    cov_matrices,
    param_names,
    param_values,
    labels,
    colors=None,
    nstd=1.0,
    figsize=(6,6),
    dpi=100
):
    """
    Generate a corner-plot style figure showing 1D Gaussians on the diagonal
    and 2D 1σ confidence ellipses off-diagonal.

    Parameters
    ----------
    cov_matrices : list of np.ndarray
        Each element is an NxN covariance matrix for the same N parameters.
    param_names : list of str
        Names of the parameters (length N).
    param_values : list or np.ndarray of shape (N,)
        Fiducial (mean) values of each parameter, used as the center for each distribution.
    labels : list of str
        A label for each covariance matrix (must match length of cov_matrices).
    colors : list of str, optional
        Colors to use for each covariance matrix. If None, defaults are used.
    nstd : float
        Number of standard deviations (1.0 => ~68%) for the 2D ellipse.
    figsize : tuple
        Figure size.
    dpi : int
        Figure DPI.

    Returns
    -------
    fig, axes : Matplotlib Figure and array of Axes
    """

    # Optionally define your own manual ranges for certain parameters:
    custom_ranges = {
        #r"$H_0$":      (67.4 - 50, 67.4 + 50),   # Example
        #r"$\sigma_8$": (0.811 - 1e-5, 0.811 + 1e-5), # Example
        #r"$S_8$": (-50, 50)
        # You can add more if you want:
        # r"$\Omega_b h^2$": (0.015, 0.03),
        # r"$\Omega_c h^2$": (0.09, 0.14),
        # r"$n_s$":       (0.90, 1.0),
        # r"$w_0$":       (-1.2, -0.8)
    }

    N = len(param_names)
    n_matrices = len(cov_matrices)

    # Default color cycle if none given
    if colors is None:
        color_cycle = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
        colors = [color_cycle[i % len(color_cycle)] for i in range(n_matrices)]

    fig, axes = plt.subplots(N, N, figsize=figsize, dpi=dpi)
    axes = np.atleast_2d(axes)  # ensure it's 2D

    # Diagonal subplots: 1D Gaussians + print sigma^2
    for i in range(N):
        ax = axes[i, i]

        for k, (cov, color) in enumerate(zip(cov_matrices, colors)):
            sigma = np.sqrt(cov[i, i])
            x_min = param_values[i] - 40*sigma
            x_max = param_values[i] + 40*sigma
            x_plot = np.linspace(x_min, x_max, 1000)
            pdf = norm.pdf(x_plot, loc=param_values[i], scale=sigma)

            ax.plot(x_plot, pdf, color=color, lw=2)
            
            # Print sigma^2 in matching color, offset each line a bit
            sigma_sq = sigma**2
            ax.text(
                0.05,
                0.90 - 0.08*k,
                rf"$\sigma^2 = {sigma_sq:.2e}$",
                color=color,
                transform=ax.transAxes
            )

        # X-limits from the first matrix or whichever you want
        sigma_ref = np.sqrt(cov_matrices[0][i, i])
        ax.set_xlim(param_values[i] - 2*sigma_ref, param_values[i] + 2*sigma_ref)
        
        # # If you have a custom range for this diagonal parameter, override it:
        # if param_names[i] in custom_ranges:
        #     ax.set_xlim(*custom_ranges[param_names[i]])

        ax.set_xlabel(param_names[i])
        ax.set_yticks([])

    # Off-diagonal subplots: 2D ellipses
    for i in range(N):
        for j in range(i+1, N):
            ax_low = axes[j, i]  # lower triangle
            ax_high = axes[i, j] # upper triangle
            ax_high.set_visible(False)  # hide mirrored subplot if you want

            # Draw each ellipse
            for (cov, color) in zip(cov_matrices, colors):
                subcov = cov[[i, j]][:, [i, j]]
                mean_ij = [param_values[i], param_values[j]]
                confidence_ellipse(ax_low, mean_ij, subcov, color=color, nstd=nstd)

            # Default: ±5σ around param_values
            sigma_i = np.sqrt(cov_matrices[0][i, i])
            sigma_j = np.sqrt(cov_matrices[0][j, j])
            x_min = param_values[i] - 1*sigma_i
            x_max = param_values[i] + 1*sigma_i
            y_min = param_values[j] - 1*sigma_j
            y_max = param_values[j] + 1*sigma_j

            ax_low.set_xlim(x_min, x_max)
            ax_low.set_ylim(y_min, y_max)

            # If these parameters have manual ranges, override them
            if param_names[i] in custom_ranges:
                ax_low.set_xlim(*custom_ranges[param_names[i]])
            if param_names[j] in custom_ranges:
                ax_low.set_ylim(*custom_ranges[param_names[j]])

            if j == N - 1:
                ax_low.set_xlabel(param_names[i])
            if i == 0:
                ax_low.set_ylabel(param_names[j])

    # Make a legend in the top-right corner
    fig.legend(labels, loc='upper right', bbox_to_anchor=(0.98, 0.98))

    plt.tight_layout()
    return fig, axes


if __name__ == "__main__":
    # Parameter info
    param_names_latex = [
        r"$H_0$",         
        r"$\Omega_b h^2$",
        r"$\Omega_c h^2$",
        r"$n_s$",         
        r"$m_\nu$",       
        r"$A_s$",         
        r"$w_0$",         
        r"$\sigma_8$",
        r"$S_8$",
        r"$\Omega_m$"
    ]

    param_values = [
        67.4,        # H_0
        0.0224,      # Omega_b h^2
        0.120,       # Omega_c h^2
        0.965,       # n_s
        0.06,        # m_nu
        2.1e-9,      # A_s
        -1.0,        # w_0
        0.811,       # sigma8
        0.811 * y**(1/2), # S8
        (0.120 + 0.0224) / 0.674**2 # Omega_m
    ]

    # Keep only the chosen indices for the final covariance
    param_names_latex_kept = [param_names_latex[i] for i in keep_indices]
    param_values_kept      = [param_values[i] for i in keep_indices]

    # Inverted sums of the reduced covariances
    cov_matrices = [
        inv(visb_c_reduced + visp_c_reduced),
        inv(visb_s_reduced + visp_s_reduced),
        inv(visb_f_reduced + visp_f_reduced)
    ]

    labels = ["CMB", "Galaxy", "CMB + Galaxy"]

    cov_matrices = [
        inv(visp_f_reduced),
        # inv(visb_f_reduced),
        inv(visb_f_reduced + visp_f_reduced)
    ]

    labels = ['powerspectra', 'power- $+$ bispectra']

    fig, axes = plot_corner(
        cov_matrices=cov_matrices,
        param_names=param_names_latex_kept,
        param_values=param_values_kept,
        labels=labels,
        colors=["tab:blue", "tab:orange", 'black'],
        nstd=1.0,
        figsize=(8,8),
        dpi=100
    )

    #print(np.linalg.inv(visp_f_reduced + visb_f_reduced))
    # plt.show()
    plt.savefig('paramconstraints_sigma8_diffspectra.pdf', dpi = 300)
