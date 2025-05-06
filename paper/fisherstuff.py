import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import norm
from numpy.linalg import inv
import matplotlib.lines as mlines

np.set_printoptions(
    precision=1,
    suppress=False, 
    formatter={'float_kind': lambda x: f"{x:.2e}"}
)

fisher_matrices_dir = '/Users/jonasfrugte/Desktop/Research_Project/fisher_calc_weak_lensing/code/fisher_matrices'

######################
# Bispectra Matrices #
######################
visb_c = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_approx_c.txt')

visb_s = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_approx_s.txt')

visb_f = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_approx_both.txt')

########################
# Powerspectra Matrices
########################
visp_c = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_c.txt')

visp_s = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_s.txt')

visp_f = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_both.txt')

######################
# Prior Matrices #
######################

planck = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_approx_c_planck.txt')

planck_prior = np.diag([0.0033, 1/(0.0005)**2, 0, 1/(0.02)**2, 0, 0, 0])

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

visb_c, visb_s, visb_f, visp_c, visp_s, visp_f, planck, planck_prior = append_row_column(
    [visb_c, visb_s, visb_f, visp_c, visp_s, visp_f, planck, planck_prior],
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

visb_c, visb_s, visb_f, visp_c, visp_s, visp_f, planck, planck_prior = append_row_column(
    [visb_c, visb_s, visb_f, visp_c, visp_s, visp_f, planck, planck_prior],
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

visb_c, visb_s, visb_f, visp_c, visp_s, visp_f, planck, planck_prior = append_row_column(
    [visb_c, visb_s, visb_f, visp_c, visp_s, visp_f, planck, planck_prior],
    omm_ders
)

omm = (0.0224 + 0.120) / (0.674**2)
sigma8 = 0.81

sigmaomm_ders = [0, 
           0,
           0,
           0,
           0,
           0,
           0,
           omm**0.25,
           0,
           0.25 * sigma8 * omm**(-0.75)]

visb_c, visb_s, visb_f, visp_c, visp_s, visp_f, planck, planck_prior = append_row_column(
    [visb_c, visb_s, visb_f, visp_c, visp_s, visp_f, planck, planck_prior],
    sigmaomm_ders
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
#keep_indices = [0, 7, 9]
keep_indices=[10]
# keep_indices = [0, 1, 2, 3, 4, 5, 6]
#keep_indices = [4, 6]

# Reduced copies for bispectra
visb_c_reduced = visb_c[np.ix_(keep_indices, keep_indices)]
visb_s_reduced = visb_s[np.ix_(keep_indices, keep_indices)]
visb_f_reduced = visb_f[np.ix_(keep_indices, keep_indices)]

# Reduced copies for powerspectra
visp_c_reduced = visp_c[np.ix_(keep_indices, keep_indices)]
visp_s_reduced = visp_s[np.ix_(keep_indices, keep_indices)]
visp_f_reduced = visp_f[np.ix_(keep_indices, keep_indices)]

planck_reduced = planck[np.ix_(keep_indices, keep_indices)]
planck_prior_reduced = planck_prior[np.ix_(keep_indices, keep_indices)]

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
        angle=angle, edgecolor=color, facecolor=color, lw=1, alpha=0.2, **kwargs
    )
    ax.add_patch(ellip)
    ellip = Ellipse(
        xy=mean, width=width, height=height,
        angle=angle, edgecolor=color, facecolor='none', lw=1, **kwargs
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
            x_min = param_values[i] - 100*sigma
            x_max = param_values[i] + 100*sigma
            x_plot = np.linspace(x_min, x_max, 5000)
            pdf = norm.pdf(x_plot, loc=param_values[i], scale=sigma)

            ax.fill_between(x_plot, pdf, color=color, lw=1, alpha = 0.2)
            ax.plot(x_plot, pdf, color=color, lw=1, alpha = 1)
            
            # Print sigma^2 in matching color, offset each line a bit
            ax.text(
                0.05,
                0.90 - 0.08*k,
                rf"${sigma:.1e}$",
                color=color,
                transform=ax.transAxes
            )

        # X-limits from the first matrix or whichever you want
        sigma_ref = np.sqrt(cov_matrices[0][i, i])
        ax.set_xlim(param_values[i] - 1.5*sigma_ref, param_values[i] + 1.5*sigma_ref)
        
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
    lines = [mlines.Line2D([], [], color = colors[i], label = labels[i]) for i in range(min(len(colors), len(labels)))]
    fig.legend(handles=lines, loc='upper right', bbox_to_anchor=(0.98, 0.98))

    plt.tight_layout()
    return fig, axes

from tabulate import tabulate

def save_table(param_names, param_vals, which_pars, constraints, labels):
    table = [[0 for i in range(len(constraints) + 1)] for j in range(len(which_pars))]
    for i in range(len(which_pars)):
        # param names and fiducial values
        table[i][0] = param_names[i]
        #table[i][1] = param_vals[i]
        #constraints
        for j in range(len(constraints)):
            table[i][j+1] = np.abs(np.sqrt(constraints[j][i][i]) / param_vals[i]) * 100
    labels = ['Par'] + labels
    print(tabulate(table, headers=labels, tablefmt='latex_raw', floatfmt='.10f'))
        
    

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
        r"$\Omega_m$",
        r"$\sigma_8\Omega_m^{0.25}$"
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
        (0.120 + 0.0224) / 0.674**2, # Omega_m
        0.811 * ((0.120 + 0.0224) / 0.674**2)**0.25
    ]

    # Keep only the chosen indices for the final covariance
    param_names_latex_kept = [param_names_latex[i] for i in keep_indices]
    param_values_kept      = [param_values[i] for i in keep_indices]

    # cov_matrices = [
    #     inv(visb_c_reduced + visp_c_reduced),
    #     inv(visb_s_reduced + visp_s_reduced),
    #     inv(visb_f_reduced + visp_f_reduced)
    # ]

    # labels = ["CMB", "Galaxy", "CMB + Galaxy"]

    # cov_matrices = [
    #     inv(visp_c_reduced),
    #     inv(visb_c_reduced),
    #     inv(visb_c_reduced + visp_c_reduced)
    # ]

    # labels = ['cmb powerp', 'cmb bisp', 'cmb power + bisp']

    # cov_matrices = [
    #     inv(visp_c_reduced + visb_c_reduced),
    #     inv(visp_s_reduced + visb_s_reduced),
    #     inv(visp_f_reduced),
    #     inv(visb_f_reduced),
    #     inv(visb_f_reduced + visp_f_reduced)
    # ]

    # labels = ['CMB + Gal Powersp', 'CMB + Gal Bisp', 'CMB Power- + Bisp', 'Gal Power- + Bisp', 'CMB + Gal Power- + Bisp']

    cov_matrices = [
        inv(planck_reduced),
        inv(planck_prior_reduced + planck_reduced)
    ]

    labels = ['Planck, no priors', 'Planck']

    fig, axes = plot_corner(
        cov_matrices=cov_matrices,
        param_names=param_names_latex_kept,
        param_values=param_values_kept,
        labels=labels,
        colors=["tab:blue", "tab:orange", 'tab:green', 'tab:red','black'],
        nstd=1.0,
        figsize=(16,16),
        dpi=100
    )

    #print(np.linalg.inv(visp_f_reduced + visb_f_reduced))
    #plt.show()
    #plt.savefig('/Users/jonasfrugte/Desktop/Research_Project/fisher_calc_weak_lensing/paper/figures/param_constraints_all.pdf', dpi = 300)

    save_table(param_names_latex_kept, param_values_kept, which_pars=keep_indices,constraints=cov_matrices, labels=labels)
