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
# Convergence (bisp)
visb_c = np.array([
    [ 1.03291163e+00,  2.06865762e-01, -1.73748438e+00, -3.36334319e+00,  2.57507000e-02, -9.02272360e-01, -3.31060410e-01],
    [ 2.06865762e-01,  4.21529585e-02, -3.47999633e-01, -6.88880813e-01,  5.15774102e-03, -1.79857024e-01, -6.60428685e-02],
    [-1.73748438e+00, -3.47999633e-01,  2.92890430e+00,  5.69889944e+00, -4.32792348e-02,  1.51535768e+00,  5.56382986e-01],
    [-3.36334319e+00, -6.88880813e-01,  5.69889944e+00,  1.12884346e+01, -8.38662127e-02,  2.92025774e+00,  1.07246662e+00],
    [ 2.57507000e-02,  5.15774102e-03, -4.32792348e-02, -8.38662127e-02,  6.41959141e-04, -2.24938576e-02, -8.25204814e-03],
    [-9.02272360e-01, -1.79857024e-01,  1.51535768e+00,  2.92025774e+00, -2.24938576e-02,  7.89137341e-01,  2.89453449e-01],
    [-3.31060410e-01, -6.60428685e-02,  5.56382986e-01,  1.07246662e+00, -8.25204814e-03,  2.89453449e-01,  1.06244427e-01]
])

# Shear (bisp)
visb_s = np.array([
    [ 1.88500828e+02,  2.23183949e+01, -1.56671938e+02, -1.21036153e+02,  2.64358967e+00, -8.11044392e+01, -6.82077843e+01],
    [ 2.23183949e+01,  2.67274426e+00, -1.86760429e+01, -1.45708725e+01,  3.14533131e-01, -9.65684766e+00, -8.07652633e+00],
    [-1.56671938e+02, -1.86760429e+01,  1.30780412e+02,  1.01553923e+02, -2.20358748e+00,  6.76397718e+01,  5.66880738e+01],
    [-1.21036153e+02, -1.45708725e+01,  1.01553923e+02,  7.98965790e+01, -1.70943608e+00,  5.24849058e+01,  4.37788446e+01],
    [ 2.64358967e+00,  3.14533131e-01, -2.20358748e+00, -1.70943608e+00,  3.71345276e-02, -1.14005023e+00, -9.56561957e-01],
    [-8.11044392e+01, -9.65684766e+00,  6.76397718e+01,  5.24849058e+01, -1.14005023e+00,  3.49837143e+01,  2.93454941e+01],
    [-6.82077843e+01, -8.07652633e+00,  5.66880738e+01,  4.37788446e+01, -9.56561957e-01,  2.93454941e+01,  2.46840703e+01]
])

# Full (bisp)
visb_f = np.array([
    [ 2.01152576e+02,  2.39436681e+01, -1.65245682e+02, -1.29401157e+02,  2.80885024e+00, -8.61094114e+01, -7.34006370e+01],
    [ 2.39436681e+01,  2.87573019e+00, -1.97888660e+01, -1.58821652e+01,  3.35826030e-01, -1.02954789e+01, -8.67647734e+00],
    [-1.65245682e+02, -1.97888660e+01,  1.37095994e+02,  1.09969714e+02, -2.32086537e+00,  7.12415855e+01,  5.98399489e+01],
    [-1.29401157e+02, -1.58821652e+01,  1.09969714e+02,  9.54901163e+01, -1.85248433e+00,  5.71387125e+01,  4.67103573e+01],
    [ 2.80885024e+00,  3.35826030e-01, -2.32086537e+00, -1.85248433e+00,  3.94024363e-02, -1.20855446e+00, -1.01924135e+00],
    [-8.61094114e+01, -1.02954789e+01,  7.12415855e+01,  5.71387125e+01, -1.20855446e+00,  3.70900980e+01,  3.11819115e+01],
    [-7.34006370e+01, -8.67647734e+00,  5.98399489e+01,  4.67103573e+01, -1.01924135e+00,  3.11819115e+01,  2.66147129e+01]
])


########################
# Powerspectra Matrices #
########################
# Convergence (powersp)
visp_c = np.array([
    [ 2.32028112e-02,  4.27071021e-03, -4.22204262e-02,  7.86457698e-03,  8.88930828e-04, -3.07148538e-02, -1.16040148e-02],
    [ 4.27071021e-03,  1.06670342e-03, -8.52286520e-03, -5.63129298e-04,  1.67250593e-04, -5.41950885e-03, -2.08184572e-03],
    [-4.22204262e-02, -8.52286520e-03,  7.88470797e-02, -8.85210795e-03, -1.62686201e-03,  5.52535403e-02,  2.09699803e-02],
    [ 7.86457698e-03, -5.63129298e-04, -8.85210795e-03,  1.79329956e-02,  2.81091927e-04, -1.22130705e-02, -4.32433568e-03],
    [ 8.88930828e-04,  1.67250593e-04, -1.62686201e-03,  2.81091927e-04,  3.41574492e-05, -1.17458092e-03, -4.43854877e-04],
    [-3.07148538e-02, -5.41950885e-03,  5.52535403e-02, -1.22130705e-02, -1.17458092e-03,  4.08735504e-02,  1.54068164e-02],
    [-1.16040148e-02, -2.08184572e-03,  2.09699803e-02, -4.32433568e-03, -4.43854877e-04,  1.54068164e-02,  5.81406135e-03]
])

# Shear (powersp)
visp_s = np.array([
    [ 1.55460255e-01,  1.30661865e-02, -1.02459167e-01, -2.72638126e-02,  1.97256033e-03, -5.70104649e-02, -6.59993897e-02],
    [ 1.30661865e-02,  1.16301229e-03, -8.92513358e-03, -2.78826083e-03,  1.69083125e-04, -4.89871624e-03, -5.55076897e-03],
    [-1.02459167e-01, -8.92513358e-03,  6.90791194e-02,  2.03701806e-02, -1.31609020e-03,  3.80909614e-02,  4.35104143e-02],
    [-2.72638126e-02, -2.78826083e-03,  2.03701806e-02,  8.64095287e-03, -3.70941937e-04,  1.07866253e-02,  1.15964010e-02],
    [ 1.97256033e-03,  1.69083125e-04, -1.31609020e-03, -3.70941937e-04,  2.51981784e-05, -7.28989310e-04, -8.37684916e-04],
    [-5.70104649e-02, -4.89871624e-03,  3.80909614e-02,  1.07866253e-02, -7.28989310e-04,  2.11095788e-02,  2.42139490e-02],
    [-6.59993897e-02, -5.55076897e-03,  4.35104143e-02,  1.15964010e-02, -8.37684916e-04,  2.42139490e-02,  2.80329426e-02]
])

# Full (powersp)
visp_f = np.array([
    [ 9.07460464e-02,  8.81447255e-03, -7.26122479e-02, -1.35071037e-02,  1.40358341e-03, -4.23964124e-02, -3.91825896e-02],
    [ 8.81447255e-03,  1.13253693e-03, -8.79681857e-03, -1.95488516e-03,  1.67301251e-04, -5.09140107e-03, -3.86507714e-03],
    [-7.26122479e-02, -8.79681857e-03,  7.42059503e-02,  8.11622779e-03, -1.45442629e-03,  4.58687115e-02,  3.22468453e-02],
    [-1.35071037e-02, -1.95488516e-03,  8.11622779e-03,  1.53624615e-02, -7.11737173e-05, -3.25625345e-04,  5.13964021e-03],
    [ 1.40358341e-03,  1.67301251e-04, -1.45442629e-03, -7.11737173e-05,  2.91364749e-05, -9.33617357e-04, -6.28252866e-04],
    [-4.23964124e-02, -5.09140107e-03,  4.58687115e-02, -3.25625345e-04, -9.33617357e-04,  3.05023641e-02,  1.91893364e-02],
    [-3.91825896e-02, -3.86507714e-03,  3.22468453e-02,  5.13964021e-03, -6.28252866e-04,  1.91893364e-02,  1.70019014e-02]
])

# normalization
param_values_normalization = np.array([
        67.4,         # H_0
        0.0224,       # Omega_b h^2
        0.120,        # Omega_c h^2
        0.965,        # n_s
        0.06,         # m_nu
        2.1e-9,       # A_s
        -1.0          # w_0
    ])
normalization_matrix = np.outer(param_values_normalization, param_values_normalization)
for mat in (visb_c, visb_s, visb_f, visp_c, visp_s, visp_f):
    mat *= normalization_matrix


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


#!!!#
keep_indices = [0, 3]
#!!!#



# Reduced copies for bispectra
visb_c_reduced = visb_c[np.ix_(keep_indices, keep_indices)]
visb_s_reduced = visb_s[np.ix_(keep_indices, keep_indices)]
visb_f_reduced = visb_f[np.ix_(keep_indices, keep_indices)]

# Reduced copies for powerspectra
visp_c_reduced = visp_c[np.ix_(keep_indices, keep_indices)]
visp_s_reduced = visp_s[np.ix_(keep_indices, keep_indices)]
visp_f_reduced = visp_f[np.ix_(keep_indices, keep_indices)]

# Optionally run the same eigen analysis on the reduced versions:
# print_sorted_eigens('b_c_reduced', visb_c_reduced)
# print_sorted_eigens('b_s_reduced', visb_s_reduced)
# print_sorted_eigens('b_f_reduced', visb_f_reduced)
# print_sorted_eigens('p_c_reduced', visp_c_reduced)
# print_sorted_eigens('p_s_reduced', visp_s_reduced)
# print_sorted_eigens('p_f_reduced', visp_f_reduced)
print_sorted_eigens('p+b_f', visb_f + visp_f)
print_sorted_eigens('p+b_f_reduced', visb_f_reduced + visp_f_reduced)

######### PLOTS #########


def confidence_ellipse(ax, mean, cov, color='blue', nstd=1.0, **kwargs):
    """
    Plot the 2D confidence ellipse of a Gaussian distribution given by mean & covariance
    onto the Matplotlib axes 'ax'.
    - mean: (2,) array of the parameter means
    - cov: (2x2) covariance matrix
    - color: color for the ellipse
    - nstd: number of standard deviations (1 => ~68% region, 2 => ~95%, etc.)
    """

    # Eigen-decompose the 2x2 covariance
    vals, vecs = np.linalg.eigh(cov)
    # The angle to rotate the ellipse
    angle = np.degrees(np.arctan2(*vecs[:, 1][::-1]))
    # Width and height of ellipse (major/minor axis)
    width, height = 2 * nstd * np.sqrt(vals)

    # Draw the ellipse
    ellip = Ellipse(xy=mean, width=width, height=height,
                    angle=angle, edgecolor=color, facecolor='none', lw=2, **kwargs)
    ax.add_patch(ellip)
    pass


def plot_corner(
    cov_matrices,
    param_names,
    param_values,
    labels,
    colors=None,
    nstd=1.0,
    figsize=(9,9),
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
        Number of standard deviations to show for the 2D ellipse (1.0 => ~68%).
    figsize : tuple
        Figure size.
    dpi : int
        Figure DPI.

    Returns
    -------
    fig, axes : Matplotlib Figure and array of Axes
    """

    # Number of parameters
    N = len(param_names)
    n_matrices = len(cov_matrices)

    if colors is None:
        # Some default color cycle (can customize)
        color_cycle = color_cycle = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
        # If fewer colors than covariance matrices, cycle repeats
        colors = [color_cycle[i % len(color_cycle)] for i in range(n_matrices)]

    fig, axes = plt.subplots(N, N, figsize=figsize, dpi=dpi)
    # Ensure axes is 2D array
    axes = np.atleast_2d(axes)

    # For each diagonal axis: plot 1D Gaussians
    for i in range(N):
        ax = axes[i, i]
        # We'll plot param range around the fiducial value
        # by estimating e.g. ±5 sqrt(cov) or so
        for cov, color in zip(cov_matrices, colors):
            sigma = np.sqrt(cov[i, i])  # std dev for parameter i
            x_min = param_values[i] - 5*sigma
            x_max = param_values[i] + 5*sigma
            x_plot = np.linspace(x_min, x_max, 300)
            # Normal distribution centered on param_values[i], scale = sigma
            pdf = norm.pdf(x_plot, loc=param_values[i], scale=sigma)
            # We'll scale the PDF so different lines are visible
            ax.plot(x_plot, pdf, color=color, lw=2)
        ax.set_xlim(
            param_values[i] - 5*np.sqrt(cov_matrices[0][i, i]),
            param_values[i] + 5*np.sqrt(cov_matrices[0][i, i])
        )
        ax.set_xlabel(param_names[i])
        ax.set_yticks([])  # no y-ticks on the diagonal

    # For each off diagonal: plot 2D ellipses
    for i in range(N):
        for j in range(i+1, N):
            ax = axes[j, i]  # lower triangle
            ax_top = axes[i, j]  # upper triangle (we can mirror or hide)
            # We'll only draw on the lower triangle (and mirror or hide the upper)
            ax_top.set_visible(False)  # hide the upper triangle if you like

            # For each covariance, add ellipse
            for cov, color in zip(cov_matrices, colors):
                # Slice out the 2x2 block
                subcov = cov[[i, j]][:, [i, j]]
                # Means
                mean_ij = [param_values[i], param_values[j]]
                confidence_ellipse(ax, mean_ij, subcov, color=color, nstd=nstd)

            # We pick a range for x / y
            x_sigma = np.sqrt(cov_matrices[0][i, i])
            y_sigma = np.sqrt(cov_matrices[0][j, j])
            ax.set_xlim(param_values[i] - 5*x_sigma, param_values[i] + 5*x_sigma)
            ax.set_ylim(param_values[j] - 5*y_sigma, param_values[j] + 5*y_sigma)
            if j == N - 1:
                ax.set_xlabel(param_names[i])
            if i == 0:
                ax.set_ylabel(param_names[j])

            # Turn on minor tick etc. if you want
            # ax.minorticks_on()

    # Legend: top right
    # We'll create a separate axis on top of figure for the legend
    fig.legend(labels, loc='upper right', bbox_to_anchor=(0.98, 0.98))

    plt.tight_layout()
    return fig, axes


if __name__ == "__main__":
    # setting up parameters
    param_names_latex = [
        r"$H_0$",         # Hubble parameter
        r"$\Omega_b h^2$",# baryon density
        r"$\Omega_c h^2$",# cold dark matter density
        r"$n_s$",         # spectral index
        r"$m_\nu$",       # sum of neutrino masses
        r"$A_s$",         # amplitude of scalar fluctuations
        r"$w_0$"          # dark energy equation of state
    ]

    param_values = [
        67.4,         # H_0
        0.0224,       # Omega_b h^2
        0.120,        # Omega_c h^2
        0.965,        # n_s
        0.06,         # m_nu
        2.1e-9,       # A_s
        -1.0          # w_0
    ]
    param_values = [0, 0, 0, 0, 0, 0, 0]

    param_names_latex_kept = [param_names_latex[i] for i in keep_indices]
    param_values_kept      = [param_values[i] for i in keep_indices]

    # Covariance matrices
    cov_matrices = [inv(visb_c_reduced + visp_c_reduced),
                    inv(visb_s_reduced + visp_s_reduced),
                    inv(visb_f_reduced + visp_f_reduced)]

    labels = ["CMB", "Galaxy", "CMB + Galaxy"]

    fig, axes = plot_corner(
        cov_matrices=cov_matrices,
        param_names=param_names_latex_kept,
        param_values=param_values_kept,
        labels=labels,
        colors=["tab:blue", "tab:orange", 'black'],  # optional
        nstd=1.0,
        figsize=(8,8),
        dpi=100
    )

    plt.show()