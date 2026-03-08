import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import norm
from numpy.linalg import inv
import matplotlib.lines as mlines
import getdist
from getdist import plots, MCSamples

# Enable LaTeX rendering for proper display of parameter names
matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['font.family'] = 'serif'

np.set_printoptions(
    precision=1,
    suppress=False, 
    formatter={'float_kind': lambda x: f"{x:.2e}"}
)

fisher_matrices_dir = '/home3/p319950/ResearchProject/fisher_calc_weak_lensing/code/fisher_matrices'

######################
# Bispectra Matrices
######################

visb_c = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_c.txt')

visb_s = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_s.txt')

visb_f = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_both.txt')

visb_c_pb = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_c_pb.txt')

visb_s_pb = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_s_pb.txt')

visb_f_pb = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_both_pb.txt')

########################
# Powerspectra Matrices
########################

# CMB anisotropy powerspectra computed with lmin = 30

visp_c = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_c.txt')

visp_s = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_s.txt')

visp_f = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_both.txt')

######################
# Planck
######################

#planck = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_c_planck.txt')

planck_prior = np.diag([
    12 / 60**2, 
    1/(0.0005)**2, 
    12 / 1**2, 
    1/(0.02)**2, 
    1e18, 
    0.063**(-2), 
    0.01, 
    0.01, 
    0.01
    ])

#planck_p_prior = planck + planck_prior

######################
# CMB Stage 4
######################

viscmb_f = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_bothCMB.txt')

viscmb_t = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_t.txt')

viscmb_e = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_e.txt')

######################
# Takada & Jain
######################

# visp_s_takada_jain = np.loadtxt(fisher_matrices_dir + '/fish_mat_powersp_s_takada_jain.txt')

# visb_s_takada_jain = np.loadtxt(fisher_matrices_dir + '/fish_mat_bisp_approx_s_takada_jain.txt')

# prior_takada_jain = np.diag([1 / 13**2, 1/(0.0002)**2, 0, 1/(0.008)**2, 1e8, 0, 0])

########### Fiducial Vals and latex names ###########

# fiducial values for each cosmological parameter (used to undo derivative
# normalization applied during Fisher computation)
fiducial_param_vals = np.array([
    67.4,         # H_0
    0.0224,       # Omega_b h^2
    0.120,        # Omega_c h^2
    0.965,        # n_s
    2.1e-9,       # A_s
    0.063,        # tau
    0.06,         # m_\nu
    -1.0,         # w_0
    7.8           # logT_AGN
])

param_names_latex = [
        r"$H_0$",         
        r"$\Omega_{b} h^2$",
        r"$\Omega_{c} h^2$",
        r"$n_s$",         
        r"$A_s$",         
        r"$\tau$",      
        r"$m_\nu$", 
        r"$w_0$",
        r"$\log T_{\mathrm{AGN}}$"
    ]

# after loading the Fisher matrices we will want to rescale them by the
# fiducial parameter values (undo the normalization applied in
# data_importer_new).  The matrices are loaded above; we perform the
# rescaling just below once the fiducial vector is defined.

def _unnormalize_matrix(mat):
    """Return mat_{ij} = mat_{ij}^{normalized} / (fid[i]*fid[j])."""
    return mat / np.outer(fiducial_param_vals, fiducial_param_vals)

########### SIGMA8 ################

def add_parameter(mats, vec):
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

# These were last updated 15 feb 2026

# Apply undo-scaling to all imported Fisher matrices so that further
# manipulations operate on un-normalized objects.
for _name in [
        'visb_c', 'visb_s', 'visb_f',
        'visb_c_pb', 'visb_s_pb', 'visb_f_pb',
        'visp_c', 'visp_s', 'visp_f',
        'viscmb_f', 'viscmb_t', 'viscmb_e'
    ]:
    if _name in globals():
        globals()[_name] = _unnormalize_matrix(globals()[_name])

s8_ders = np.array([
    2.80598548e-03, # H_0
    -6.56385299e+00, # Omega_b h^2 
    4.32895218e+00,  # Omega_c h^2
    3.09054725e-01,# n_s
    1.91861249e+08, # A_s
    -6.92875026e-04, # tau
    -6.51946137e-02, # m_\nu
    -2.02959792e-01, # w_0
    0               # logT_AGN # is 0 because 8Mpc h^-1 is much larger than the scales affected by AGN feedback, so we expect no effect on sigma8
    ])

y = (0.0224 + 0.120) / (0.3 * 0.674**2) # used to calculate S8
S8_ders = [
    (0.811 * y**(-1/2) * 0.5 * y * (0.674)**(-1) * (-2) * 0.01), 
    (0.811 * y**(-1/2) * (0.3 * 0.674**2)**(-1) * 0.5),
    (0.811 * y**(-1/2) * (0.3 * 0.674**2)**(-1) * 0.5),
    0,
    0,
    0,
    0,
    0,
    0,
    y**(1/2)]

omm_ders = [
    (0.811 * y**(-1/2) * 0.5 * y * (0.674)**(-1) * (-2) * 0.01), 
    1 / 0.674**2,
    1 / 0.674**2,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0]

omm = (0.0224 + 0.120) / (0.674**2)
sigma8 = 0.811

sigmaomm_ders = [
    0, 
    0,
    0,
    0,
    0,
    0,
    0,
    (omm**0.25),
    0,
    0,
    0,
    (0.25 * sigma8 * omm**(-0.75))]


def process_fishes(
        fishes, 
        which_indices_to_keep, 
        derived_param_derivss = [s8_ders, S8_ders, omm_ders, sigmaomm_ders], 
        derived_param_fid_vals = [0.811, 0.811 * (0.0224 + 0.120) / (0.674**2), (0.0224 + 0.120) / (0.674**2), 0.811 * ((0.0224 + 0.120) / (0.674**2))**(1/4)], 
        derived_param_names = [r"$\sigma_8$", r"$S_8$", r"$\Omega_m$", r"$\sigma_8 \Omega_m^{0.25}$"]
        ):
    
    fishes_inv = [inv(fish) for fish in fishes]
    fiducial_param_vals_kept = fiducial_param_vals
    param_names_latex_kept = param_names_latex

    # extend to include derived params
    for i in range(len(derived_param_derivss)):
        fishes_inv = add_parameter(fishes_inv, derived_param_derivss[i])
        fiducial_param_vals_kept = np.append(fiducial_param_vals_kept, derived_param_fid_vals[i])
        param_names_latex_kept.append(derived_param_names[i])

    if which_indices_to_keep == 1:
        keep_indices = [6, 7, 9, 11] # "tight"
    if which_indices_to_keep == 2:
        keep_indices = [0, 1, 2, 3, 5, 6, 7, 8] # "lcdm"
    if which_indices_to_keep == 3:
        keep_indices = [4, 6]

    # marginalize over unwanted params
    for i in range(len(fishes)):
        fishes_inv[i] = fishes_inv[i][np.ix_(keep_indices, keep_indices)]
    fiducial_param_vals_kept = [fiducial_param_vals_kept[i] for i in keep_indices]
    param_names_latex_kept = [param_names_latex_kept[i] for i in keep_indices]

    return fishes_inv, fiducial_param_vals_kept, param_names_latex_kept

def plot_confidence_ellipses(
        fishes,
        fish_names, 
        which_indices_to_keep, 
        derived_param_derivss = [s8_ders, S8_ders, omm_ders, sigmaomm_ders], 
        derived_param_fid_vals = [0.811, 0.811 * (0.0224 + 0.120) / (0.674**2), (0.0224 + 0.120) / (0.674**2), 0.811 * ((0.0224 + 0.120) / (0.674**2))**(1/4)], 
        derived_param_names = [r"$\sigma_8$", r"$S_8$", r"$\Omega_m$", r"$\sigma_8 \Omega_m^{0.25}$"]
):
    fishes_inv, fiducial_param_vals_kept, param_names_latex_kept = process_fishes(fishes, which_indices_to_keep, derived_param_derivss, derived_param_fid_vals, derived_param_names)

    random_state = np.random.default_rng(seed=42)
    
    sampss = []
    for i in range(len(fishes_inv)):
        # Generate samples from the covariance matrix
        try:            
            samples = random_state.multivariate_normal(
                mean=fiducial_param_vals_kept, 
                cov=fishes_inv[i], 
                size=100000
            )
        except RuntimeWarning("covariance is not symmetric positive-semidefinite"):
            print(fishes_inv[i])
                  
        # Create MCSamples object
        samps = MCSamples(
            samples=samples, 
            names=param_names_latex_kept, 
            label=fish_names[i],
            ignore_rows=0.0
        )
        sampss.append(samps)
    
    # Create the triangle plot once with all samples
    g = plots.get_subplot_plotter()
    g.triangle_plot(sampss, filled=True)
    
    return g

def select_plot_type(fish_pond_number):
    '''
    fish pond numbers correspond to:
    0: for the table, strong prior
    1: for the table, weak prior
    2, 3, 4, 5: for plots, tight parameters
    6, 7, 8, 9: for plots, lcdm parameters
    0.5: for the table, strong prior, with pb corrections
    1.5: for the table, weak prior, with pb corrections
    '''

    plt_name = ''
    which_indices_to_keep = 1
    f_sky = 0.5

    if fish_pond_number == 0:
        fish_matrices = [
            f_sky * viscmb_f,
            f_sky * (viscmb_f + visp_c),
            f_sky * (viscmb_f + visb_c),
            f_sky * (viscmb_f + visp_c + visb_c),
            f_sky * (viscmb_f + visp_s),
            f_sky * (viscmb_f + visb_s),
            f_sky * (viscmb_f + visp_s + visb_s),
            f_sky * (viscmb_f + visp_f),
            f_sky * (viscmb_f + visb_f),
            f_sky * (viscmb_f + visp_f + visb_f)
        ]

        labels = ['prior', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$']


    if fish_pond_number == 1:
        fish_matrices = [
            planck_prior,
            planck_prior + f_sky * (visp_c),
            planck_prior + f_sky * (visb_c),
            planck_prior + f_sky * (visp_c + visb_c),
            planck_prior + f_sky * (visp_s),
            planck_prior + f_sky * (visb_s),
            planck_prior + f_sky * (visp_s + visb_s),
            planck_prior + f_sky * (visp_f),
            planck_prior + f_sky * (visb_f),
            planck_prior + f_sky * (visp_f + visb_f),
        ]

        labels = ['prior', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$']

    if fish_pond_number == 0.5:
        fish_matrices = [
            f_sky * viscmb_f,
            f_sky * (viscmb_f + visp_c),
            f_sky * (viscmb_f + visb_c_pb),
            f_sky * (viscmb_f + visp_c + visb_c_pb),
            f_sky * (viscmb_f + visp_s),
            f_sky * (viscmb_f + visb_s_pb),
            f_sky * (viscmb_f + visp_s + visb_s_pb),
            f_sky * (viscmb_f + visp_f),
            f_sky * (viscmb_f + visb_f_pb),
            f_sky * (viscmb_f + visp_f + visb_f_pb)
        ]

        labels = ['prior', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$']


    if fish_pond_number == 1.5:
        fish_matrices = [
            planck_prior,
            planck_prior + f_sky * (visp_c),
            planck_prior + f_sky * (visb_c_pb),
            planck_prior + f_sky * (visp_c + visb_c_pb),
            planck_prior + f_sky * (visp_s),
            planck_prior + f_sky * (visb_s_pb),
            planck_prior + f_sky * (visp_s + visb_s_pb),
            planck_prior + f_sky * (visp_f),
            planck_prior + f_sky * (visb_f_pb),
            planck_prior + f_sky * (visp_f + visb_f_pb),
        ]

        labels = ['prior', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$', r'$C_\ell$', r'$B_{\ell_1\ell_2\ell_3}$', r'$C_\ell + B_{\ell_1\ell_2\ell_3}$']

    if fish_pond_number == 2:
        fish_matrices = [
            planck_prior + f_sky * visp_c,
            planck_prior + f_sky * visb_c,
            planck_prior + f_sky * (visp_c + visb_c),
            planck_prior + f_sky * (visp_f + visb_f)
        ]

        labels = [r'weak prior + $C_\ell^{\psi_{\text{CMB}}}$', r'weak prior + $B_{\ell_1\ell_2\ell_3}^{\psi_{\text{CMB}}}$', r'weak prior + $C_\ell^{\psi_{\text{CMB}}} + B_{\ell_1\ell_2\ell_3}^{\psi_{\text{CMB}}}$', 'weak prior + All lensing spec.']

        plt_name = 'param_constraints_tight_cmb_weak_prior.pdf'
        which_indices_to_keep = 1


    if fish_pond_number == 3:
        fish_matrices = [
            planck_prior + f_sky * visp_s,
            planck_prior + f_sky * visb_s,
            planck_prior + f_sky * (visp_s + visb_s),
            planck_prior + f_sky * (visp_f + visb_f),
        ]

        labels = [r'weak prior + $C_\ell^{\psi_{\text{gal}}}$', r'weak prior + $B_{\ell_1\ell_2\ell_3}^{\psi_{\text{gal}}}$', r'weak prior + $C_\ell^{\psi_{\text{gal}}} + B_{\ell_1\ell_2\ell_3}^{\psi_{\text{gal}}}$', 'weak prior + All lensing spec.']

        plt_name = 'param_constraints_tight_gal_weak_prior.pdf'
        which_indices_to_keep = 1

    if fish_pond_number == 4:
        fish_matrices = [
            f_sky * (viscmb_f + visp_c),
            f_sky * (viscmb_f + visb_c),
            f_sky * (viscmb_f + visp_c + visb_c),
            f_sky * (viscmb_f + visp_f + visb_f)
        ]

        labels = [r'CMB $T+E$ prior + $C_\ell^{\psi_{\text{CMB}}}$', r'CMB $T+E$ prior + $B_{\ell_1\ell_2\ell_3}^{\psi_{\text{CMB}}}$', r'CMB $T+E$ prior + $C_\ell^{\psi_{\text{CMB}}} + B_{\ell_1\ell_2\ell_3}^{\psi_{\text{CMB}}}$', 'CMB $T+E$ prior + All lensing spec.']

        plt_name = 'param_constraints_tight_cmb_cmb_prior.pdf'
        which_indices_to_keep = 1


    if fish_pond_number == 5:
        fish_matrices = [
            f_sky * (viscmb_f + visp_s),
            f_sky * (viscmb_f + visb_s),
            f_sky * (viscmb_f + visp_s + visb_s),
            f_sky * (viscmb_f + visp_f + visb_f)
        ]

        labels = [r'CMB $T+E$ prior + $C_\ell^{\psi_{\text{gal}}}$', r'CMB $T+E$ prior + $B_{\ell_1\ell_2\ell_3}^{\psi_{\text{gal}}}$', r'CMB $T+E$ prior + $C_\ell^{\psi_{\text{gal}}} + B_{\ell_1\ell_2\ell_3}^{\psi_{\text{gal}}}$', 'CMB $T+E$ prior + All lensing spec.']

        plt_name = 'param_constraints_tight_gal_cmb_prior.pdf'
        which_indices_to_keep = 1


    if fish_pond_number == 6:
        fish_matrices = [
            planck_prior + f_sky * visp_c,
            planck_prior + f_sky * visb_c,
            planck_prior + f_sky * (visp_c + visb_c),
            planck_prior + f_sky * (visp_f + visb_f)
        ]

        labels = [r'weak prior + $C_\ell^{\psi_{\text{CMB}}}$', r'weak prior + $B_{\ell_1\ell_2\ell_3}^{\psi_{\text{CMB}}}$', r'weak prior + $C_\ell^{\psi_{\text{CMB}}} + B_{\ell_1\ell_2\ell_3}^{\psi_{\text{CMB}}}$', 'weak prior + All lensing spec.']

        plt_name = 'param_constraints_lcdm_cmb_weak_prior.pdf'
        which_indices_to_keep = 2


    if fish_pond_number == 7:
        fish_matrices = [
            planck_prior + f_sky * visp_s,
            planck_prior + f_sky * visb_s,
            planck_prior + f_sky * (visp_s + visb_s),
            planck_prior + f_sky * (visp_f + visb_f),
        ]

        labels = [r'weak prior + $C_\ell^{\psi_{\text{gal}}}$', r'weak prior + $B_{\ell_1\ell_2\ell_3}^{\psi_{\text{gal}}}$', r'weak prior + $C_\ell^{\psi_{\text{gal}}} + B_{\ell_1\ell_2\ell_3}^{\psi_{\text{gal}}}$', 'weak prior + All lensing spec.']

        plt_name = 'param_constraints_lcdm_gal_weak_prior.pdf'
        which_indices_to_keep = 2

    if fish_pond_number == 8:
        fish_matrices = [
            f_sky * (viscmb_f + visp_c),
            f_sky * (viscmb_f + visb_c),
            f_sky * (viscmb_f + visp_c + visb_c),
            f_sky * (viscmb_f + visp_f + visb_f)
        ]

        labels = [r'CMB $T+E$ prior + $C_\ell^{\psi_{\text{CMB}}}$', r'CMB $T+E$ prior + $B_{\ell_1\ell_2\ell_3}^{\psi_{\text{CMB}}}$', r'CMB $T+E$ prior + $C_\ell^{\psi_{\text{CMB}}} + B_{\ell_1\ell_2\ell_3}^{\psi_{\text{CMB}}}$', 'CMB $T+E$ prior + All lensing spec.']

        plt_name = 'param_constraints_lcdm_cmb_cmb_prior.pdf'
        which_indices_to_keep = 2


    if fish_pond_number == 9:
        fish_matrices = [
            f_sky * (viscmb_f + visp_s),
            f_sky * (viscmb_f + visb_s),
            f_sky * (viscmb_f + visp_s + visb_s),
            f_sky * (viscmb_f + visp_f + visb_f)
        ]

        labels = [r'CMB $T+E$ prior + $C_\ell^{\psi_{\text{gal}}}$', r'CMB $T+E$ prior + $B_{\ell_1\ell_2\ell_3}^{\psi_{\text{gal}}}$', r'CMB $T+E$ prior + $C_\ell^{\psi_{\text{gal}}} + B_{\ell_1\ell_2\ell_3}^{\psi_{\text{gal}}}$', 'CMB $T+E$ prior + All lensing spec.']

        plt_name = 'param_constraints_lcdm_gal_cmb_prior.pdf'
        which_indices_to_keep = 2
    
    return fish_matrices, labels, plt_name, which_indices_to_keep

def create_plots(plot_type_numbers):
    for fish_pond_number in plot_type_numbers:
        fish_matrices, labels, plt_name, which_indices_to_keep = select_plot_type(fish_pond_number)
        g = plot_confidence_ellipses(fish_matrices, labels, which_indices_to_keep)
        plt.savefig('/home3/p319950/ResearchProject/fisher_calc_weak_lensing/paper/figures/' + plt_name, dpi = 300, bbox_inches="tight")
        print('created: ', plt_name)

if __name__ == "__main__":
    plots_to_make = [2, 3, 4, 5, 6, 7, 8, 9]
    create_plots(plots_to_make)

# from tabulate import tabulate

# def save_table(param_names, param_vals, which_pars, constraints, labels):
#     table = [[0 for i in range(len(constraints) + 1)] for j in range(len(which_pars))]
#     for i in range(len(which_pars)):
#         # param names and fiducial values
#         table[i][0] = param_names[i]
#         #table[i][1] = param_vals[i]
#         #constraints
#         for j in range(len(constraints)):
#             table[i][j+1] = np.abs(np.sqrt(constraints[j][i][i]))
#     labels = ['Par'] + labels
#     print(tabulate(table, headers=labels, tablefmt='latex_raw', floatfmt='.3e'))
        
    

# if __name__ == "__main__":
#     # Parameter info
#     param_names_latex = [
#         r"$H_0$",         
#         r"$\Omega_b h^2$",
#         r"$\Omega_c h^2$",
#         r"$n_s$",         
#         r"$m_\nu$", 
#         r"$\tau$",      
#         r"$A_s$",         
#         r"$w_0$",         
#         r"$\sigma_8$",
#         r"$S_8$",
#         r"$\Omega_m$",
#         r"$\sigma_8\Omega_m^{0.25}$"
#     ]

#     param_values = [
#         67.4,        # H_0
#         0.0224,      # Omega_b h^2
#         0.120,       # Omega_c h^2
#         0.965,       # n_s
#         0.06,        # m_nu
#         0.063,       # tau
#         2.1e-9,      # A_s
#         -1.0,        # w_0
#         0.811,       # sigma8
#         0.811 * y**(1/2), # S8
#         (0.120 + 0.0224) / 0.674**2, # Omega_m
#         0.811 * ((0.120 + 0.0224) / 0.674**2)**0.25
#     ]