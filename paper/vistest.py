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
# Planck #
######################
planck = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_c_planck.txt')

planck_prior = np.diag([12 / 60**2, 1/(0.0005)**2, 12 / 1**2, 1/(0.02)**2, 1e6, 1e19, 1e6])

planck_p_prior = planck + planck_prior

viscmb_f = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_both_cmb.txt')

######################
# Takada & Jain #
######################

visp_s_takada_jain = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_s_takada_jain.txt')

visb_s_takada_jain = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_approx_s_takada_jain.txt')

prior_takada_jain = np.diag([1 / 13**2, 1/(0.0002)**2, 0, 1/(0.008)**2, 1e8, 0, 0])

# pars = [b'H', b'ombh2', b'omch2', b'ns', b'mnu', b'As', b'w0']


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


keep_indices = [0, 1]
#keep_indices=[10]
#keep_indices = [0, 1, 2, 3, 4, 5, 6]
#keep_indices = [4, 6]

def process_fishes(fishes, keep_indices, derived_param_derivss = []):
    fishes_inv = [inv(fish) for fish in fishes]
    
    # marginalize over unwanted params
    for i in range(len(fishes)):
        fishes_inv[i] = fishes_inv[i][np.ix_(keep_indices, keep_indices)]
    
    return fishes_inv

# print(0.1**(-0.5) * np.sqrt(process_fishes(
#     [prior_takada_jain + visp_s_takada_jain, prior_takada_jain + visp_s_takada_jain, prior_takada_jain + visb_s_takada_jain + visb_s_takada_jain],
#     [6],   
#     derived_param_derivss = [s8_ders, S8_ders, omm_ders, sigmaomm_ders]
#     )))

######### PLOTS #########

plt_scale = 1.5

# new version that returns flipped version which I think is correct
def confidence_ellipse(ax, mean, cov, color='blue', nstd=1.0, isbig=False, **kwargs):
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
    if isbig == True:
        width_small, height_small = (0.5 * width, 0.5 * height)
        ellip = Ellipse(
        xy=mean, width=width_small, height=height_small,
        angle=angle, edgecolor=color, facecolor='none', lw=1, linestyle = '--', **kwargs
        )
        ax.add_patch(ellip)

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
        ax.set_xlim(param_values[i] - plt_scale*sigma_ref, param_values[i] + plt_scale*sigma_ref)
        
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
                if color == "tab:blue":
                    print('yay')
                    confidence_ellipse(ax_low, mean_ij, subcov, color=color, nstd=nstd, isbig=True)
                else:
                    confidence_ellipse(ax_low, mean_ij, subcov, color=color, nstd=nstd)

            # Default: ±σ around param_values
            sigma_i = np.sqrt(cov_matrices[0][i, i])
            sigma_j = np.sqrt(cov_matrices[0][j, j])
            x_min = param_values[i] - plt_scale*sigma_i
            x_max = param_values[i] + plt_scale*sigma_i
            y_min = param_values[j] - plt_scale*sigma_j
            y_max = param_values[j] + plt_scale*sigma_j

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
    print(tabulate(table, headers=labels, tablefmt='latex_raw', floatfmt='.2f'))
        
    

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
        0.0224      # Omega_b h^2
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

    # fish_matrices = [
    #     0.5 * viscmb_f,
    #     0.5 * (viscmb_f + 1.5 * visp_c),
    #     0.5 * (viscmb_f + 1.5 * visb_c),
    #     0.5 * (viscmb_f + 1.5 * visb_c + 1.5 * visp_c)
    # ]

    # labels = ['cmb T+E', 'cmb powersp', 'cmb bisp', 'cmb powersp + bisp']


    # fish_matrices = [
    #     visp_c + visb_c,
    #     visp_s + visb_s,
    #     visp_f,
    #     visb_f,
    #     visb_f + visp_f
    # ]

    # labels = ['CMB + Gal Powersp', 'CMB + Gal Bisp', 'CMB Power- + Bisp', 'Gal Power- + Bisp', 'CMB + Gal Power- + Bisp']
    vis1 = np.array([[1,-0.49],[-0.49,0.5]])
    vis2 = np.array([[1,-0.99],[-0.99,1]])

    fish_matrices = [
        vis1,
        vis2,
        vis1 + vis2
    ]

    labels = ['CMB + Gal Powersp', 'CMB + Gal Bisp', 'CMB Power- + Bisp']

    # fish_matrices = [
    #     planck,
    #     planck_prior,
    #     planck_prior + planck
    # ]

    # labels = ['Planck, no priors', 'only prior', 'Planck']

    cov_matrices = process_fishes(fish_matrices, keep_indices)

    fig, axes = plot_corner(
        cov_matrices=cov_matrices,
        param_names=param_names_latex_kept,
        param_values=param_values_kept,
        labels=labels,
        colors=["tab:blue", "tab:orange", 'tab:green', 'tab:red','black'],
        nstd=1.0,
        figsize=(6,6),
        dpi=100
    )



    #print(np.linalg.inv(visp_f_reduced + visb_f_reduced))
    plt.show()
    #plt.savefig('/Users/jonasfrugte/Desktop/Research_Project/fisher_calc_weak_lensing/paper/figures/param_constraints_tight.pdf', dpi = 300)

    save_table(param_names_latex_kept, param_values_kept, which_pars=keep_indices,constraints=cov_matrices, labels=labels)
