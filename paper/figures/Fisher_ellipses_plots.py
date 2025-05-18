import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import norm
from numpy.linalg import inv

# # Define fiducial parameter values
# fiducial_values = {
#     r"$H_0$": 67.4,
#     r"$\Omega_b h^2$": 0.0224,
#     r"$\Omega_c h^2$": 0.120,
#     r"$A_s$": 2.1e-9,
#     r"$n_s$": 0.965,
#     r"$\Sigma m_\nu$": 0.06,
# }

fiducial_values = {
    r"$H_0$": 67.4,
    r"$\Omega_b h^2$": 0.0224,
    r"$\Omega_c h^2$": 0.120,
    r"$n_s$": 0.965
}

fiducial_array = np.array([
    fiducial_values[param] for param in fiducial_values.keys()
])

# Function to normalize covariance values
def normalize_covariance(covariance_matrix, fiducial_array):
    normalization_matrix = np.outer(fiducial_array, fiducial_array)
    return (covariance_matrix / normalization_matrix) * 100

# Function to plot normalized confidence ellipses
def plot_multiple_fisher_ellipses_with_percentages_no_text(fisher_matrices, parameter_names, labels, confidence_level=0.68):
    num_params = len(parameter_names)
    chi2_val = {0.68: 2.3, 0.95: 5.99}[confidence_level]
    colors = ['blue', 'red', 'green', 'purple', 'orange']

    fig, axes = plt.subplots(num_params, num_params, figsize=(12, 12))
    for i in range(num_params):
        for j in range(num_params):
            ax = axes[i, j]
            if i == j:
                # Plot 1D Gaussians for each Fisher matrix on the diagonal
                for k, fisher_matrix in enumerate(fisher_matrices):
                    covariance_matrix = np.linalg.inv(fisher_matrix)
                    print([covariance_matrix[i, i] for i in range(num_params)])
                    covariance_normalized = normalize_covariance(covariance_matrix, fiducial_array)
                    sigma = np.sqrt(covariance_normalized[i, i]) #!!!
                    x = np.linspace(-3 * sigma, 3 * sigma, 1000)
                    y = np.exp(-0.5 * (x / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
                    ax.plot(x, y, label=labels[k], color=colors[k % len(colors)], lw=2)
                ax.set_xlim(-3 * sigma, 3 * sigma)
                ax.set_yticks([])
                if i == num_params - 1:
                    ax.set_xlabel(parameter_names[i])
                if i == 0:
                    ax.legend()
            elif i < j:
                ax.axis('off')
            else:
                # Plot 2D confidence ellipses for each Fisher matrix
                for k, fisher_matrix in enumerate(fisher_matrices):
                    covariance_matrix = np.linalg.inv(fisher_matrix)
                    covariance_normalized = normalize_covariance(covariance_matrix, fiducial_array)
                    sub_cov = covariance_normalized[np.ix_([i, j], [i, j])]
                    eigenvalues, eigenvectors = np.linalg.eig(sub_cov)
                    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
                    width, height = 2 * np.sqrt(chi2_val * eigenvalues) #!!!
                    ellipse = Ellipse(xy=(0, 0), width=width, height=height, angle=angle,
                                      edgecolor=colors[k % len(colors)], fill=False, lw=2, label=labels[k])
                    ax.add_patch(ellipse)
                ax.set_xlim(-15, 15)
                ax.set_ylim(-15, 15)
                ax.axhline(0, color='gray', linestyle='--', lw=0.5)
                ax.axvline(0, color='gray', linestyle='--', lw=0.5)
                if j == 0:
                    ax.set_ylabel(parameter_names[i])
                if i == num_params - 1:
                    ax.set_xlabel(parameter_names[j])

    plt.tight_layout()
    plt.show()

Fisher_conv = np.array([[ 2.27552541e-02,  9.52373123e+00, -2.34613171e+01,
         7.85964847e-01,  8.94470394e-01],
       [ 9.52373123e+00,  4.39944531e+03, -1.01984287e+04,
         2.68961182e+02,  3.76008972e+02],
       [-2.34613171e+01, -1.01984287e+04,  2.47482363e+04,
        -7.77074829e+02, -9.21555786e+02],
       [ 7.85964847e-01,  2.68961182e+02, -7.77074829e+02,
         3.83012772e+01,  3.08532257e+01],
       [ 8.94470394e-01,  3.76008972e+02, -9.21555786e+02,
         3.08532257e+01,  3.51004562e+01]])[:-1, :-1]
Fisher_shear = np.array([[ 3.65161806e-01,  1.37609497e+02, -2.23485840e+02,
        -6.03931236e+00,  7.29286051e+00],
       [ 1.37609497e+02,  5.27060742e+04, -8.44470078e+04,
        -2.35060327e+03,  2.74217090e+03],
       [-2.23485840e+02, -8.44470078e+04,  1.37204047e+05,
         3.72692871e+03, -4.39039111e+03],
       [-6.03931236e+00, -2.35060327e+03,  3.72692871e+03,
         1.19891571e+02, -1.21028679e+02],
       [ 7.29286051e+00,  2.74217090e+03, -4.39039111e+03,
        -1.21028679e+02,  1.45192963e+02]])[:-1, :-1]
Fisher_full = np.array([[ 7.83514130e-01,  2.93683159e+02, -5.03949975e+02,
        -1.11990863e+01,  1.64104009e+01],
       [ 2.93683159e+02,  1.13348283e+05, -1.89850897e+05,
        -4.37220182e+03,  6.23645730e+03],
       [-5.03949975e+02, -1.89850897e+05,  3.36733074e+05,
         6.48228083e+03, -1.10115871e+04],
       [-1.11990863e+01, -4.37220182e+03,  6.48228083e+03,
         3.59674275e+02, -1.85644498e+02],
       [ 1.64104009e+01,  6.23645730e+03, -1.10115871e+04,
        -1.85644498e+02,  3.62955147e+02]])[:-1, :-1]
# Parameter names
# param_names = [
#     r"$H_0$", r"$\Omega_b h^2$", r"$\Omega_c h^2$", r"$A_s$", r"$n_s$", r"$\Sigma m_\nu$"
# ]
# param_names = [
#     r"$H_0$", r"$\Omega_b h^2$", r"$\Omega_c h^2$", r"$n_s$", r"$\Sigma m_\nu$"
# ]

param_names = [
    r"$H_0$", r"$\Omega_b h^2$", r"$\Omega_c h^2$", r"$n_s$"
]

# Experiment labels
labels = ["CMB Lensing", "Galaxy Lensing", "CMB + Galaxy Lensing"]

# Plotting
plot_multiple_fisher_ellipses_with_percentages_no_text(
    [Fisher_conv, Fisher_shear, Fisher_full],
    param_names,
    labels
)

########### PLOTTING ############

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
        color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
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
    # Example usage:
    # Suppose we have 3 parameters: [H0, omega_b, omega_c]
    param_names = ["H0", "ωb", "ωc"]
    param_values = [70.0, 0.022, 0.12]  # some fiducials

    # Two (3x3) covariance matrices
    cov_matrices = [inv(cov1), cov2]

    labels = ["Experiment A", "Experiment B"]

    fig, axes = plot_corner(
        cov_matrices=cov_matrices,
        param_names=param_names,
        param_values=param_values,
        labels=labels,
        colors=["tab:blue", "tab:orange"],  # optional
        nstd=1.0,
        figsize=(8,8),
        dpi=100
    )

    plt.show()


