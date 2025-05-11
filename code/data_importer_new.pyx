print('wahoo!')
import json
cimport interpolation as interp
import math
import numpy as np
from libc.stdio cimport printf
from libc.math cimport fmax, fabs

from libc.stdlib cimport getenv

from interp_settings import exp_par

# declare interpolation parameters, NEEDS TO BE SAME AS INTERPOLATION 
cdef double k_min = exp_par['k_min']
cdef double k_max = exp_par['k_max'] 
cdef int k_num = exp_par['k_num'] # * 2 # for finer data option
cdef int k_num_fine = k_num * 25
cdef double chi_min = exp_par['chi_min']
cdef double chi_max = exp_par['chi_max']
cdef int chi_num = exp_par['chi_num'] # * 2 # for finer data option
cdef double z_min = exp_par['z_min']
cdef double z_max = exp_par['z_max']
cdef int z_num = exp_par['z_num'] # * 2 # for finer data option
cdef int z_num_fine = z_num * 25

def get_k_max():
    return k_max

# specifies the noise that will be used for cmb lensing
cdef int cmb_noise_type = 2
cdef int galaxy_noise_type = 2
# CMB type values correspond to:
# 0: S0 noise curves (old, remove)
# 1: sigma = 1, Delta_P = 6 (stage 3 toshiya)
# 2: sigma = 3, Delta_P = 1 (stage 4 toshiya)
# 3: sigma = 5, Delta_T = 30, Delta_P = 52 (planck (double check))

# (units are in microKelvin arcmin and arcmin)

cdef int lmax_cmbps = 3000

# SO noise
cdef str folder_file_path = '/scratch/p319950/data/'
cdef str filepath_convergence_noise_file_path = '/scratch/p319950/data_rough/' + 'conv_noise.dat'
conv_noise_data_array = np.loadtxt(filepath_convergence_noise_file_path)
cdef double[:, :] conv_noise_data = conv_noise_data_array

# conv noise from quadratic estimator
cdef double[:] ls_cmbn = np.loadtxt('cmb_noise_files/ls_1_3000_64.txt')
ls_cmbn_np_array = np.loadtxt('cmb_noise_files/ls_1_3000_64.txt') # need to also have a np array version to easily convert noise values from convergence to lens potential below
cdef double lmin_cmbn = ls_cmbn[0]
cdef int lnum_cmbn = len(ls_cmbn)
cdef double lmax_cmbn = ls_cmbn[lnum_cmbn - 1]

cdef double[:] cmbn_301 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma3_DeltaT0_DeltaP1.txt')) / (ls_cmbn_np_array * (ls_cmbn_np_array + 1))
cdef double[:] cmbn_106 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma1_DeltaT0_DeltaP6.txt')) / (ls_cmbn_np_array * (ls_cmbn_np_array + 1))
cdef double[:] planck_noise = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma_planck.txt')) / (ls_cmbn_np_array * (ls_cmbn_np_array + 1))

# perhaps we can export data as c arrays instead of as memory views in the future so that we can specify return types like below and can avoid having to 
# declare types of all data before every data import
#cdef (dict[str, double], double, double[:, :], double[:, :], double[:], double[:], double[:], double[:], double[:], double[:], double[:, :], double[:]) data_importer(str folder_name):
def data_import_func(folder_name):
    cdef str filepath = folder_file_path + folder_name
    print(filepath)

    cdef double[:] cosm_par = np.load(filepath + '/cosm_par.npy')
    cdef double C
    cdef double[:, :] a_data = np.load(filepath + '/a.npy')
    cdef double[:, :] b_data = np.load(filepath + '/b.npy')
    cdef double[:, :] c_data = np.load(filepath + '/c.npy')
    cdef double[:] lps_cc_data = np.load(filepath + '/lensing_power_spectrum_cc.npy')
    cdef double[:] lps_cs_data = np.load(filepath + '/lensing_power_spectrum_cs.npy')
    cdef double[:] lps_ss_data = np.load(filepath + '/lensing_power_spectrum_ss.npy')
    cdef double[:] scale_factor_data = np.loadtxt(filepath + '/scale_factor')
    cdef double[:] window_c_data = np.load(filepath + '/window_func_c.npy')
    cdef double[:] window_s_data = np.load(filepath + '/window_func_s.npy')
    cdef double[:, :] mps_data = np.load(filepath + '/matter_power_spectrum.npy')
    cdef double[:] z_at_chi_data = np.load(filepath + '/z_at_chi.npy')
    cdef double[:, :] cmbps = np.load(filepath + '/cmb_ps_with_ls.npy')

    # all 2d funcs are made such that order of arguments is (k, z)
    C = 0.
    with open(filepath + '/rho_bar', 'r') as file:
        C = float(file.read())

    return cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data, cmbps



# turns data into functions 

cpdef double a(double k, double z, double[:, :] a_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, a_data)
cpdef double b(double k, double z, double[:, :] b_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, b_data)
cpdef double c(double k, double z, double[:, :] c_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, c_data)

# without gil cython doesn't like strings bc they are python objects, you can use char* C objects (C like strings)

cpdef double lps(double l, char* type1, char* type2, double[:] lps_cc_data, double[:] lps_cs_data, double[:] lps_ss_data, double[:, :] cmbps) noexcept nogil:
    # CAUTION: this function doesn't throw error messages if type1 or type2 are not of type convergence or shear
    # for optimization purposes
    if type1[0] == b'c' and type2[0] == b'c':
        return interp.logspace_linear_interp(l, k_min, k_max, k_num_fine, lps_cc_data)
    if (type1[0] == b'c' and type2[0] == b's') or (type1[0] == b's' and type2[0] == b'c'):
        return interp.logspace_linear_interp(l, k_min, k_max, k_num_fine, lps_cs_data)
    if type1[0] == b's' and type2[0] == b's':
        return interp.logspace_linear_interp(l, k_min, k_max, k_num_fine, lps_ss_data)
    if type1[0] == b't' and type2[0] == b't':
        return interp.linear_interp(l, cmbps[0, 0], lmax_cmbps, lmax_cmbps + 1, cmbps[:, 1])
    if type1[0] == b'e' and type2[0] == b'e':
        return interp.linear_interp(l, cmbps[0, 0], lmax_cmbps, lmax_cmbps + 1, cmbps[:, 2])
    if type1[0] == b'b' and type2[0] == b'b':
        return interp.linear_interp(l, cmbps[0, 0], lmax_cmbps, lmax_cmbps + 1, cmbps[:, 3])
    if (type1[0] == b't' and type2[0] == b'e') or (type1[0] == b'e' and type2[0] == b't'):
        return interp.linear_interp(l, cmbps[0, 0], lmax_cmbps, lmax_cmbps + 1, cmbps[:, 4])
    else:
        return 0

cdef double scale_factor(double chi, double[:] scale_factor_data) noexcept nogil:    
    return interp.linear_interp(chi, chi_min, chi_max, chi_num, scale_factor_data)

cdef double window_func(double chi, char* type, double[:] window_c_data, double[:] window_s_data) noexcept nogil:
    if type[0] == b'c':
        return fmax(0., interp.logspace_linear_interp(chi, chi_min, chi_max, chi_num, window_c_data))
    elif type[0] == b's':
        return fmax(0., interp.logspace_linear_interp(chi, chi_min, chi_max, chi_num, window_s_data))

# gives slightly different results than scipy RegularGridInterpolator for some obscure reason, however values on 
# grid points seem to agree with the imported data
cpdef double matter_power_spectrum(double k, double z, double[:, :] mps_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, mps_data)

cdef double z_at_chi(double chi, double[:] z_at_chi_data) noexcept nogil:
    return interp.linear_interp(chi, chi_min, chi_max, chi_num, z_at_chi_data)

# turns functions into bi- and powerspectra functions

cdef double law_cosines(double x, double y, double z) noexcept nogil:
    # gives cosine of angle between vector x and y, where we know the magnitudes of x, y, z and that x + y + z = 0 vector
    # 27/03/25 added -1 factor
    return -1 * (x**2 + y**2 - z**2) / (2 * x * y)

cpdef double F_2(double k1, double k2, double k3, double z_input, double[:, :] a_data, double[:, :] b_data, double[:, :] c_data) noexcept nogil:
    # cosine factors in term 2 and term 3 are determined by the delta^(3)(k_1, k_2, k_3) factor combined with law of cosines
    # this is why we also need k3, it is only used in the calculations to determine the angle between k1 and k2
    cdef double term_1 = (5. / 7.) * a(k1, z_input, a_data) * a(k2, z_input, a_data)
    cdef double term_2 = 0.5 * law_cosines(k1, k2, k3) * (k1 / k2 + k2 / k1) * b(k1, z_input, b_data) * b(k2, z_input, b_data)
    cdef double term_3 = (2. / 7.) * law_cosines(k1, k2, k3)**2 * c(k1, z_input, c_data) * c(k2, z_input, c_data)
    return term_1 + term_2 + term_3

cpdef double mbs(double k1, double k2, double k3, double z, double[:, :] mps_data, double[:, :] a_data, double[:, :] b_data, double[:, :] c_data) noexcept nogil:
    # precalculate to avoid having to calculate twice
    cdef double mps_k1 = matter_power_spectrum(k1, z, mps_data)
    cdef double mps_k2 = matter_power_spectrum(k2, z, mps_data)
    cdef double mps_k3 = matter_power_spectrum(k3, z, mps_data)

    cdef double term1 = 2. * F_2(k1, k2, k3, z, a_data, b_data, c_data) * mps_k1 * mps_k2
    cdef double term2 = 2. * F_2(k2, k3, k1, z, a_data, b_data, c_data) * mps_k2 * mps_k3
    cdef double term3 = 2. * F_2(k3, k1, k2, z, a_data, b_data, c_data) * mps_k3 * mps_k1
    return term1 + term2 + term3

cdef inline double lbs_integrand(
    double chi, double k1, double k2, double k3, 
    char* type1, char* type2, char* type3,
    double[:, :] a_data, double[:, :] b_data, double[:, :] c_data, 
    double[:] scale_factor_data, double[:] window_c_data, double[:] window_s_data, 
    double[:, :] mps_data, double[:] z_at_chi_data) noexcept nogil:

    cdef double integrand_val
    cdef double factor1, factor2, factor3, factor4, factor5, factor6

    # Compute individual factors
    factor1 = chi**2
    factor2 = scale_factor(chi, scale_factor_data)**(-3)
    factor3 = window_func(chi, type1, window_c_data, window_s_data)
    factor4 = window_func(chi, type2, window_c_data, window_s_data)
    factor5 = window_func(chi, type3, window_c_data, window_s_data)
    factor6 = mbs(k1 / chi, k2 / chi, k3 / chi, z_at_chi(chi, z_at_chi_data), 
                  mps_data, a_data, b_data, c_data)

    # Compute integrand value
    integrand_val = factor1 * factor2 * factor3 * factor4 * factor5 * factor6

    # Check if integrand is <0 (should never be the case) and print factors if true
    # if integrand_val < 0.0:
    # printf("Factor 1 (chi^2): %f\n", factor1)
    # printf("Factor 2 (scale_factor^(-3)): %f\n", factor2)
    # printf("Factor 3 (window_func type1): %f\n", factor3)
    # printf("Factor 4 (window_func type2): %f\n", factor4)
    # printf("Factor 5 (window_func type3): %f\n", factor5)
    # printf("Factor 6 (mbs): %f\n", factor6)
    # printf("Integrand value: %f\n", integrand_val)

    return integrand_val

# cdef inline double lbs_integrand(double chi, double k1, double k2, double k3, char* type1, char* type2, char* type3, double[:, :] a_data, double[:, :] b_data, double[:, :] c_data, double[:] scale_factor_data, double[:] window_c_data, double[:] window_s_data, double[:, :] mps_data, double[:] z_at_chi_data) noexcept nogil:
#     cdef double integrand_val
#     integrand_val = chi**2 * scale_factor(chi, scale_factor_data)**(-3) * window_func(chi, type1, window_c_data, window_s_data) * window_func(chi, type2, window_c_data, window_s_data) * window_func(chi, type3, window_c_data, window_s_data) * mbs(k1 / chi, k2 / chi, k3 / chi, z_at_chi(chi, z_at_chi_data), mps_data, a_data, b_data, c_data)
#     if integrand_val < 0.:

#     return integrand_val

cdef double lbs_integral(double k1, double k2, double k3, char* type1, char* type2, char* type3, int num_samples, double C, double[:, :] a_data, double[:, :] b_data, double[:, :] c_data, double[:] scale_factor_data, double[:] window_c_data, double[:] window_s_data, double[:, :] mps_data, double[:] z_at_chi_data) noexcept nogil:
    cdef double int_width = (chi_max - chi_min) / (num_samples - 1)
    cdef double result = 0
    cdef int i
    for i in range(1, num_samples - 1):
        result += lbs_integrand(chi_min + int_width * i, k1, k2, k3, type1, type2, type3, a_data, b_data, c_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data)
        # for boundary values use values that lie *just* inside the boundaries to prevent some nasty errors
    result += 0.5 * (
        lbs_integrand(chi_min * 1.01, k1, k2, k3, type1, type2, type3, a_data, b_data, c_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data) 
        + lbs_integrand(chi_max - 1.0, k1, k2, k3, type1, type2, type3, a_data, b_data, c_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data)
        )
    result *= int_width

    return result

cdef extern from "math.h":
    double sqrt(double x) nogil
    double log(double x) nogil
    double exp(double x) nogil

cdef extern from "complex.h":
    double creal(double complex x) nogil

# cdef double wigner_3j_approx_nocheck(int l1, int l2, int l3) noexcept nogil:
#     cdef double complex L_half, factor, term1, term2, term3, term1_pow, term2_pow, term3_pow, denominator
#     cdef int L = l1 + l2 + l3

#     # does not check if L is even or triangle inequalities, does so in lbs_f and lbs_der directly instead to save on computing bispec if result should be zero anyway

#     L_half = L / 2.0
    
#     # Common factors
#     factor = (-1)**L_half * (2 * 3.141592653589793)**(-0.5) * exp(3.0 / 2) * (L + 2)**(-0.25)
    
#     # Power term for each fraction
#     term1 = (L_half - l1 + 0.5) / (L_half - l1 + 1)
#     term2 = (L_half - l2 + 0.5) / (L_half - l2 + 1)
#     term3 = (L_half - l3 + 0.5) / (L_half - l3 + 1)

#     # Raising to required powers
#     term1_pow = term1 ** (L_half - l1 + 1.0 / 4.0)
#     term2_pow = term2 ** (L_half - l2 + 1.0 / 4.0)
#     term3_pow = term3 ** (L_half - l3 + 1.0 / 4.0)

#     # The denominator terms
#     denominator = ((L_half - l1 + 1)**0.25 * (L_half - l2 + 1)**0.25 * (L_half - l3 + 1)**0.25)
    
#     # The final result
#     return creal(factor * term1_pow * term2_pow * term3_pow / denominator)

cdef double wigner_3j_approx_nocheck(int l1, int l2, int l3) noexcept nogil:
    cdef double L = ( l1 + l2 + l3 ) / 2

    # does not check if 2L is even or triangle inequalities, does so in lbs_f and lbs_der directly instead to save on computing bispec if result should be zero anyway

    cdef complex factor = (-1)**L * sqrt(2.718**3 / (2 * 3.1415)) * (L + 1)**(-0.25)

    cdef complex term1 = (L-l1+1)**(-0.25) * ( (L-l1+0.5) / (L-l1+1) )**(L-l1+0.25)
    cdef complex term2 = (L-l2+1)**(-0.25) * ( (L-l2+0.5) / (L-l2+1) )**(L-l2+0.25)
    cdef complex term3 = (L-l3+1)**(-0.25) * ( (L-l3+0.5) / (L-l3+1) )**(L-l3+0.25)

    return creal(factor * term1 * term2 * term3)


# fiducial (!) full sky lensing bispectrum
cpdef double lbs(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, double C, double[:, :] a_data, double[:, :] b_data, double[:, :] c_data, double[:] scale_factor_data, double[:] window_c_data, double[:] window_s_data, double[:, :] mps_data, double[:] z_at_chi_data) noexcept nogil:
    cdef double wigner_factor
    cdef double sqrt_factor
    cdef double const_factor
    cdef double fraction_factor
    cdef double integration_result

    if (k1 + k2 + k3) % 2 == 0 and k3 <= k1 + k2 and k1 - k2 <= k3 and k2 - k1 <= k3:

        integration_result = lbs_integral(k1, k2, k3, type1, type2, type3, num_samples, C, a_data, b_data, c_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data)

        wigner_factor = fabs(wigner_3j_approx_nocheck(k1, k2, k3))
        sqrt_factor = sqrt((2.0*k1 + 1.0)*(2.0*k2 + 1.0)*(2.0*k3 + 1.0)/(4 * 3.14159)) # (pi)
        fraction_factor = 1 / (1.0 * k1 ** 2 * k2 ** 2 * k3 ** 2) # 1.0 factor to convert to floats
        const_factor = C**3 * 8 # sketchy

        # printf("interp: Tuple (l1, l2, l3): (%d, %d, %d)\n", k1, k2, k3)
        # printf("interp: Integration result: %f\n", integration_result)
        # printf("interp: Wigner Factor: %f\n", wigner_factor)
        # printf("interp: Sqrt Factor: %f\n", sqrt_factor)
        # printf("interp: Fraction Factor: %e\n", fraction_factor)
        # printf("interp: Const Factor: %e\n", const_factor)

        return wigner_factor * sqrt_factor * fraction_factor * const_factor * integration_result
    else:
        return 0.0

# fiducial values

cdef double[:] cosm_par_f
cdef double C_f
cdef double[:, :] a_data_f = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_f = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_f = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_f = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_f = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_f = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_f = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_f = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_f = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_f = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_f = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_f = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data, cmbps_f

cosm_par_f, C_f, a_data_f, b_data_f, c_data_f, lps_cc_data_f, lps_cs_data_f, lps_ss_data_f, scale_factor_data_f, window_c_data_f, window_s_data_f, mps_data_f, z_at_chi_data_f, cmbps_f = data_import_func('data_fiducial')

cdef double lbs_f(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_f, a_data_f, b_data_f, c_data_f, scale_factor_data_f, window_c_data_f, window_s_data_f, mps_data_f, z_at_chi_data_f)

cdef double lps_f(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_f, lps_cs_data_f, lps_ss_data_f, cmbps_f)

print(lmax_cmbps)
print(cmbps_f[0, 0], cmbps_f[3000, 0], cmbps_f[50, 1])

cdef double cmbps_noise(double l, double sigma, double Delta_X) noexcept nogil:
    # units to input:
    # sigma: arcmin
    # Delta_T, Delta_P: microKelvin arcmin

    # constant time!
    cdef double Tcmb = 2.728e6 # in micro kelvin
    cdef double arcmintorad = 3.14 / 10800
    sigma = sigma * arcmintorad # in radians
    Delta_X = Delta_X * arcmintorad
    
    return (Delta_X / Tcmb)**2 * exp(l * (l + 1) * sigma**2 / (8 * log(2)))

cdef double lps_noise(int l, char* type1, char* type2) noexcept nogil:
    cdef float noise = 0.

    if type1[0] == b'c' and type2[0] == b'c':
        if cmb_noise_type == 0:
            noise = 4. * (l * 1.0)**(-2) * (l + 1.0)**(-2) * conv_noise_data[l-2, 7]
        elif cmb_noise_type == 1:
            # stage 3 wide toshiya
            noise = interp.logspace_linear_interp(l, lmin_cmbn, lmax_cmbn, lnum_cmbn, cmbn_106)
        elif cmb_noise_type == 2:
            # stage 4 toshiya
            noise = interp.logspace_linear_interp(l, lmin_cmbn, lmax_cmbn, lnum_cmbn, cmbn_301)
        elif cmb_noise_type == 3:
            noise = interp.logspace_linear_interp(l, lmin_cmbn, lmax_cmbn, lnum_cmbn, planck_noise)

    if type1[0] == b's' and type2[0] == b's':
        if galaxy_noise_type == 1:
            noise = 4. * ((l - 1.) * l * (l + 1.) * (l + 2.))**(-1) * 0.3**2 / 5
        if galaxy_noise_type == 2:
            noise = 4. * ((l - 1.) * l * (l + 1.) * (l + 2.))**(-1) * 0.3**2 / 30 # 8 * 10.0 ** (-10) # this should actually be in sterradain, but that gives wrong results so for now we have it like this
        if galaxy_noise_type == 3:
            noise = 4. * ((l - 1.) * l * (l + 1.) * (l + 2.))**(-1) * 0.4**2 / 100 # noise from https://arxiv.org/pdf/astro-ph/0310125

    cdef double sigma
    cdef double Delta_X

    if type1[0] == b't' and type2[0] == b't':
        if cmb_noise_type == 1:
            # stage 3 wide toshiya
            sigma = 1
            Delta_X = 6
        elif cmb_noise_type == 2:
            # stage 4 toshiya
            sigma = 3
            Delta_X = 1
        elif cmb_noise_type == 3:
            sigma = 5
            Delta_X = 30
        
        noise = cmbps_noise(l, sigma, Delta_X)

    if type1[0] == b'e' and type2[0] == b'e':
        if cmb_noise_type == 1:
            # stage 3 wide toshiya
            sigma = 1
            Delta_X = 6 * 1.41
        elif cmb_noise_type == 2:
            # stage 4 toshiya
            sigma = 3
            Delta_X = 1 * 1.41
        elif cmb_noise_type == 3:
            sigma = 5
            Delta_X = 52
        
        noise = cmbps_noise(l, sigma, Delta_X)

    if type1[0] == b'b' and type2[0] == b'b':
        if cmb_noise_type == 1:
            # stage 3 wide toshiya
            sigma = 1
            Delta_X = 6 * 1.41
        elif cmb_noise_type == 2:
            # stage 4 toshiya
            sigma = 3
            Delta_X = 1 * 1.41
        elif cmb_noise_type == 3:
            sigma = 5
            Delta_X = 52
        
        noise = cmbps_noise(l, sigma, Delta_X)

    return noise

cdef double lps_f_obs(int l, char* type1, char* type2) noexcept nogil:

    #########################
    # NOISE: ON
    #########################

    return lps(l, type1, type2, lps_cc_data_f, lps_cs_data_f, lps_ss_data_f, cmbps_f) + lps_noise(l, type1, type2)

##########################################



cdef double[:] cosm_par_H_p_2m
cdef double C_H_p_2m
cdef double[:, :] a_data_H_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_2m, C_H_p_2m, a_data_H_p_2m, b_data_H_p_2m, c_data_H_p_2m, lps_cc_data_H_p_2m, lps_cs_data_H_p_2m, lps_ss_data_H_p_2m, scale_factor_data_H_p_2m, window_c_data_H_p_2m, window_s_data_H_p_2m, mps_data_H_p_2m, z_at_chi_data_H_p_2m, cmbps_H_p_2m = data_import_func('data_H_p_2m')

cdef double lbs_H_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_2m, a_data_H_p_2m, b_data_H_p_2m, c_data_H_p_2m, scale_factor_data_H_p_2m, window_c_data_H_p_2m, window_s_data_H_p_2m, mps_data_H_p_2m, z_at_chi_data_H_p_2m)

cdef double lps_H_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_2m, lps_cs_data_H_p_2m, lps_ss_data_H_p_2m, cmbps_H_p_2m)
            

cdef double[:] cosm_par_H_p_1m
cdef double C_H_p_1m
cdef double[:, :] a_data_H_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_1m, C_H_p_1m, a_data_H_p_1m, b_data_H_p_1m, c_data_H_p_1m, lps_cc_data_H_p_1m, lps_cs_data_H_p_1m, lps_ss_data_H_p_1m, scale_factor_data_H_p_1m, window_c_data_H_p_1m, window_s_data_H_p_1m, mps_data_H_p_1m, z_at_chi_data_H_p_1m, cmbps_H_p_1m = data_import_func('data_H_p_1m')

cdef double lbs_H_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_1m, a_data_H_p_1m, b_data_H_p_1m, c_data_H_p_1m, scale_factor_data_H_p_1m, window_c_data_H_p_1m, window_s_data_H_p_1m, mps_data_H_p_1m, z_at_chi_data_H_p_1m)

cdef double lps_H_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_1m, lps_cs_data_H_p_1m, lps_ss_data_H_p_1m, cmbps_H_p_1m)
            

cdef double[:] cosm_par_H_p_0
cdef double C_H_p_0
cdef double[:, :] a_data_H_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_0, C_H_p_0, a_data_H_p_0, b_data_H_p_0, c_data_H_p_0, lps_cc_data_H_p_0, lps_cs_data_H_p_0, lps_ss_data_H_p_0, scale_factor_data_H_p_0, window_c_data_H_p_0, window_s_data_H_p_0, mps_data_H_p_0, z_at_chi_data_H_p_0, cmbps_H_p_0 = data_import_func('data_H_p_0')

cdef double lbs_H_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_0, a_data_H_p_0, b_data_H_p_0, c_data_H_p_0, scale_factor_data_H_p_0, window_c_data_H_p_0, window_s_data_H_p_0, mps_data_H_p_0, z_at_chi_data_H_p_0)

cdef double lps_H_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_0, lps_cs_data_H_p_0, lps_ss_data_H_p_0, cmbps_H_p_0)
            

cdef double[:] cosm_par_H_p_1p
cdef double C_H_p_1p
cdef double[:, :] a_data_H_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_1p, C_H_p_1p, a_data_H_p_1p, b_data_H_p_1p, c_data_H_p_1p, lps_cc_data_H_p_1p, lps_cs_data_H_p_1p, lps_ss_data_H_p_1p, scale_factor_data_H_p_1p, window_c_data_H_p_1p, window_s_data_H_p_1p, mps_data_H_p_1p, z_at_chi_data_H_p_1p, cmbps_H_p_1p = data_import_func('data_H_p_1p')

cdef double lbs_H_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_1p, a_data_H_p_1p, b_data_H_p_1p, c_data_H_p_1p, scale_factor_data_H_p_1p, window_c_data_H_p_1p, window_s_data_H_p_1p, mps_data_H_p_1p, z_at_chi_data_H_p_1p)

cdef double lps_H_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_1p, lps_cs_data_H_p_1p, lps_ss_data_H_p_1p, cmbps_H_p_1p)
            

cdef double[:] cosm_par_H_p_2p
cdef double C_H_p_2p
cdef double[:, :] a_data_H_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_2p, C_H_p_2p, a_data_H_p_2p, b_data_H_p_2p, c_data_H_p_2p, lps_cc_data_H_p_2p, lps_cs_data_H_p_2p, lps_ss_data_H_p_2p, scale_factor_data_H_p_2p, window_c_data_H_p_2p, window_s_data_H_p_2p, mps_data_H_p_2p, z_at_chi_data_H_p_2p, cmbps_H_p_2p = data_import_func('data_H_p_2p')

cdef double lbs_H_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_2p, a_data_H_p_2p, b_data_H_p_2p, c_data_H_p_2p, scale_factor_data_H_p_2p, window_c_data_H_p_2p, window_s_data_H_p_2p, mps_data_H_p_2p, z_at_chi_data_H_p_2p)

cdef double lps_H_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_2p, lps_cs_data_H_p_2p, lps_ss_data_H_p_2p, cmbps_H_p_2p)
            

cdef double[:] cosm_par_H_m_2m
cdef double C_H_m_2m
cdef double[:, :] a_data_H_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_2m, C_H_m_2m, a_data_H_m_2m, b_data_H_m_2m, c_data_H_m_2m, lps_cc_data_H_m_2m, lps_cs_data_H_m_2m, lps_ss_data_H_m_2m, scale_factor_data_H_m_2m, window_c_data_H_m_2m, window_s_data_H_m_2m, mps_data_H_m_2m, z_at_chi_data_H_m_2m, cmbps_H_m_2m = data_import_func('data_H_m_2m')

cdef double lbs_H_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_2m, a_data_H_m_2m, b_data_H_m_2m, c_data_H_m_2m, scale_factor_data_H_m_2m, window_c_data_H_m_2m, window_s_data_H_m_2m, mps_data_H_m_2m, z_at_chi_data_H_m_2m)

cdef double lps_H_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_2m, lps_cs_data_H_m_2m, lps_ss_data_H_m_2m, cmbps_H_m_2m)
            

cdef double[:] cosm_par_H_m_1m
cdef double C_H_m_1m
cdef double[:, :] a_data_H_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_1m, C_H_m_1m, a_data_H_m_1m, b_data_H_m_1m, c_data_H_m_1m, lps_cc_data_H_m_1m, lps_cs_data_H_m_1m, lps_ss_data_H_m_1m, scale_factor_data_H_m_1m, window_c_data_H_m_1m, window_s_data_H_m_1m, mps_data_H_m_1m, z_at_chi_data_H_m_1m, cmbps_H_m_1m = data_import_func('data_H_m_1m')

cdef double lbs_H_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_1m, a_data_H_m_1m, b_data_H_m_1m, c_data_H_m_1m, scale_factor_data_H_m_1m, window_c_data_H_m_1m, window_s_data_H_m_1m, mps_data_H_m_1m, z_at_chi_data_H_m_1m)

cdef double lps_H_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_1m, lps_cs_data_H_m_1m, lps_ss_data_H_m_1m, cmbps_H_m_1m)
            

cdef double[:] cosm_par_H_m_0
cdef double C_H_m_0
cdef double[:, :] a_data_H_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_0, C_H_m_0, a_data_H_m_0, b_data_H_m_0, c_data_H_m_0, lps_cc_data_H_m_0, lps_cs_data_H_m_0, lps_ss_data_H_m_0, scale_factor_data_H_m_0, window_c_data_H_m_0, window_s_data_H_m_0, mps_data_H_m_0, z_at_chi_data_H_m_0, cmbps_H_m_0 = data_import_func('data_H_m_0')

cdef double lbs_H_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_0, a_data_H_m_0, b_data_H_m_0, c_data_H_m_0, scale_factor_data_H_m_0, window_c_data_H_m_0, window_s_data_H_m_0, mps_data_H_m_0, z_at_chi_data_H_m_0)

cdef double lps_H_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_0, lps_cs_data_H_m_0, lps_ss_data_H_m_0, cmbps_H_m_0)
            

cdef double[:] cosm_par_H_m_1p
cdef double C_H_m_1p
cdef double[:, :] a_data_H_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_1p, C_H_m_1p, a_data_H_m_1p, b_data_H_m_1p, c_data_H_m_1p, lps_cc_data_H_m_1p, lps_cs_data_H_m_1p, lps_ss_data_H_m_1p, scale_factor_data_H_m_1p, window_c_data_H_m_1p, window_s_data_H_m_1p, mps_data_H_m_1p, z_at_chi_data_H_m_1p, cmbps_H_m_1p = data_import_func('data_H_m_1p')

cdef double lbs_H_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_1p, a_data_H_m_1p, b_data_H_m_1p, c_data_H_m_1p, scale_factor_data_H_m_1p, window_c_data_H_m_1p, window_s_data_H_m_1p, mps_data_H_m_1p, z_at_chi_data_H_m_1p)

cdef double lps_H_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_1p, lps_cs_data_H_m_1p, lps_ss_data_H_m_1p, cmbps_H_m_1p)
            

cdef double[:] cosm_par_H_m_2p
cdef double C_H_m_2p
cdef double[:, :] a_data_H_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_2p, C_H_m_2p, a_data_H_m_2p, b_data_H_m_2p, c_data_H_m_2p, lps_cc_data_H_m_2p, lps_cs_data_H_m_2p, lps_ss_data_H_m_2p, scale_factor_data_H_m_2p, window_c_data_H_m_2p, window_s_data_H_m_2p, mps_data_H_m_2p, z_at_chi_data_H_m_2p, cmbps_H_m_2p = data_import_func('data_H_m_2p')

cdef double lbs_H_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_2p, a_data_H_m_2p, b_data_H_m_2p, c_data_H_m_2p, scale_factor_data_H_m_2p, window_c_data_H_m_2p, window_s_data_H_m_2p, mps_data_H_m_2p, z_at_chi_data_H_m_2p)

cdef double lps_H_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_2p, lps_cs_data_H_m_2p, lps_ss_data_H_m_2p, cmbps_H_m_2p)
            

cdef double[:] cosm_par_ombh2_p_2m
cdef double C_ombh2_p_2m
cdef double[:, :] a_data_ombh2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_2m, C_ombh2_p_2m, a_data_ombh2_p_2m, b_data_ombh2_p_2m, c_data_ombh2_p_2m, lps_cc_data_ombh2_p_2m, lps_cs_data_ombh2_p_2m, lps_ss_data_ombh2_p_2m, scale_factor_data_ombh2_p_2m, window_c_data_ombh2_p_2m, window_s_data_ombh2_p_2m, mps_data_ombh2_p_2m, z_at_chi_data_ombh2_p_2m, cmbps_ombh2_p_2m = data_import_func('data_ombh2_p_2m')

cdef double lbs_ombh2_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_2m, a_data_ombh2_p_2m, b_data_ombh2_p_2m, c_data_ombh2_p_2m, scale_factor_data_ombh2_p_2m, window_c_data_ombh2_p_2m, window_s_data_ombh2_p_2m, mps_data_ombh2_p_2m, z_at_chi_data_ombh2_p_2m)

cdef double lps_ombh2_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_2m, lps_cs_data_ombh2_p_2m, lps_ss_data_ombh2_p_2m, cmbps_ombh2_p_2m)
            

cdef double[:] cosm_par_ombh2_p_1m
cdef double C_ombh2_p_1m
cdef double[:, :] a_data_ombh2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_1m, C_ombh2_p_1m, a_data_ombh2_p_1m, b_data_ombh2_p_1m, c_data_ombh2_p_1m, lps_cc_data_ombh2_p_1m, lps_cs_data_ombh2_p_1m, lps_ss_data_ombh2_p_1m, scale_factor_data_ombh2_p_1m, window_c_data_ombh2_p_1m, window_s_data_ombh2_p_1m, mps_data_ombh2_p_1m, z_at_chi_data_ombh2_p_1m, cmbps_ombh2_p_1m = data_import_func('data_ombh2_p_1m')

cdef double lbs_ombh2_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_1m, a_data_ombh2_p_1m, b_data_ombh2_p_1m, c_data_ombh2_p_1m, scale_factor_data_ombh2_p_1m, window_c_data_ombh2_p_1m, window_s_data_ombh2_p_1m, mps_data_ombh2_p_1m, z_at_chi_data_ombh2_p_1m)

cdef double lps_ombh2_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_1m, lps_cs_data_ombh2_p_1m, lps_ss_data_ombh2_p_1m, cmbps_ombh2_p_1m)
            

cdef double[:] cosm_par_ombh2_p_0
cdef double C_ombh2_p_0
cdef double[:, :] a_data_ombh2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_0, C_ombh2_p_0, a_data_ombh2_p_0, b_data_ombh2_p_0, c_data_ombh2_p_0, lps_cc_data_ombh2_p_0, lps_cs_data_ombh2_p_0, lps_ss_data_ombh2_p_0, scale_factor_data_ombh2_p_0, window_c_data_ombh2_p_0, window_s_data_ombh2_p_0, mps_data_ombh2_p_0, z_at_chi_data_ombh2_p_0, cmbps_ombh2_p_0 = data_import_func('data_ombh2_p_0')

cdef double lbs_ombh2_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_0, a_data_ombh2_p_0, b_data_ombh2_p_0, c_data_ombh2_p_0, scale_factor_data_ombh2_p_0, window_c_data_ombh2_p_0, window_s_data_ombh2_p_0, mps_data_ombh2_p_0, z_at_chi_data_ombh2_p_0)

cdef double lps_ombh2_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_0, lps_cs_data_ombh2_p_0, lps_ss_data_ombh2_p_0, cmbps_ombh2_p_0)
            

cdef double[:] cosm_par_ombh2_p_1p
cdef double C_ombh2_p_1p
cdef double[:, :] a_data_ombh2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_1p, C_ombh2_p_1p, a_data_ombh2_p_1p, b_data_ombh2_p_1p, c_data_ombh2_p_1p, lps_cc_data_ombh2_p_1p, lps_cs_data_ombh2_p_1p, lps_ss_data_ombh2_p_1p, scale_factor_data_ombh2_p_1p, window_c_data_ombh2_p_1p, window_s_data_ombh2_p_1p, mps_data_ombh2_p_1p, z_at_chi_data_ombh2_p_1p, cmbps_ombh2_p_1p = data_import_func('data_ombh2_p_1p')

cdef double lbs_ombh2_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_1p, a_data_ombh2_p_1p, b_data_ombh2_p_1p, c_data_ombh2_p_1p, scale_factor_data_ombh2_p_1p, window_c_data_ombh2_p_1p, window_s_data_ombh2_p_1p, mps_data_ombh2_p_1p, z_at_chi_data_ombh2_p_1p)

cdef double lps_ombh2_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_1p, lps_cs_data_ombh2_p_1p, lps_ss_data_ombh2_p_1p, cmbps_ombh2_p_1p)
            

cdef double[:] cosm_par_ombh2_p_2p
cdef double C_ombh2_p_2p
cdef double[:, :] a_data_ombh2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_2p, C_ombh2_p_2p, a_data_ombh2_p_2p, b_data_ombh2_p_2p, c_data_ombh2_p_2p, lps_cc_data_ombh2_p_2p, lps_cs_data_ombh2_p_2p, lps_ss_data_ombh2_p_2p, scale_factor_data_ombh2_p_2p, window_c_data_ombh2_p_2p, window_s_data_ombh2_p_2p, mps_data_ombh2_p_2p, z_at_chi_data_ombh2_p_2p, cmbps_ombh2_p_2p = data_import_func('data_ombh2_p_2p')

cdef double lbs_ombh2_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_2p, a_data_ombh2_p_2p, b_data_ombh2_p_2p, c_data_ombh2_p_2p, scale_factor_data_ombh2_p_2p, window_c_data_ombh2_p_2p, window_s_data_ombh2_p_2p, mps_data_ombh2_p_2p, z_at_chi_data_ombh2_p_2p)

cdef double lps_ombh2_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_2p, lps_cs_data_ombh2_p_2p, lps_ss_data_ombh2_p_2p, cmbps_ombh2_p_2p)
            

cdef double[:] cosm_par_ombh2_m_2m
cdef double C_ombh2_m_2m
cdef double[:, :] a_data_ombh2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_2m, C_ombh2_m_2m, a_data_ombh2_m_2m, b_data_ombh2_m_2m, c_data_ombh2_m_2m, lps_cc_data_ombh2_m_2m, lps_cs_data_ombh2_m_2m, lps_ss_data_ombh2_m_2m, scale_factor_data_ombh2_m_2m, window_c_data_ombh2_m_2m, window_s_data_ombh2_m_2m, mps_data_ombh2_m_2m, z_at_chi_data_ombh2_m_2m, cmbps_ombh2_m_2m = data_import_func('data_ombh2_m_2m')

cdef double lbs_ombh2_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_2m, a_data_ombh2_m_2m, b_data_ombh2_m_2m, c_data_ombh2_m_2m, scale_factor_data_ombh2_m_2m, window_c_data_ombh2_m_2m, window_s_data_ombh2_m_2m, mps_data_ombh2_m_2m, z_at_chi_data_ombh2_m_2m)

cdef double lps_ombh2_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_2m, lps_cs_data_ombh2_m_2m, lps_ss_data_ombh2_m_2m, cmbps_ombh2_m_2m)
            

cdef double[:] cosm_par_ombh2_m_1m
cdef double C_ombh2_m_1m
cdef double[:, :] a_data_ombh2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_1m, C_ombh2_m_1m, a_data_ombh2_m_1m, b_data_ombh2_m_1m, c_data_ombh2_m_1m, lps_cc_data_ombh2_m_1m, lps_cs_data_ombh2_m_1m, lps_ss_data_ombh2_m_1m, scale_factor_data_ombh2_m_1m, window_c_data_ombh2_m_1m, window_s_data_ombh2_m_1m, mps_data_ombh2_m_1m, z_at_chi_data_ombh2_m_1m, cmbps_ombh2_m_1m = data_import_func('data_ombh2_m_1m')

cdef double lbs_ombh2_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_1m, a_data_ombh2_m_1m, b_data_ombh2_m_1m, c_data_ombh2_m_1m, scale_factor_data_ombh2_m_1m, window_c_data_ombh2_m_1m, window_s_data_ombh2_m_1m, mps_data_ombh2_m_1m, z_at_chi_data_ombh2_m_1m)

cdef double lps_ombh2_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_1m, lps_cs_data_ombh2_m_1m, lps_ss_data_ombh2_m_1m, cmbps_ombh2_m_1m)
            

cdef double[:] cosm_par_ombh2_m_0
cdef double C_ombh2_m_0
cdef double[:, :] a_data_ombh2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_0, C_ombh2_m_0, a_data_ombh2_m_0, b_data_ombh2_m_0, c_data_ombh2_m_0, lps_cc_data_ombh2_m_0, lps_cs_data_ombh2_m_0, lps_ss_data_ombh2_m_0, scale_factor_data_ombh2_m_0, window_c_data_ombh2_m_0, window_s_data_ombh2_m_0, mps_data_ombh2_m_0, z_at_chi_data_ombh2_m_0, cmbps_ombh2_m_0 = data_import_func('data_ombh2_m_0')

cdef double lbs_ombh2_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_0, a_data_ombh2_m_0, b_data_ombh2_m_0, c_data_ombh2_m_0, scale_factor_data_ombh2_m_0, window_c_data_ombh2_m_0, window_s_data_ombh2_m_0, mps_data_ombh2_m_0, z_at_chi_data_ombh2_m_0)

cdef double lps_ombh2_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_0, lps_cs_data_ombh2_m_0, lps_ss_data_ombh2_m_0, cmbps_ombh2_m_0)
            

cdef double[:] cosm_par_ombh2_m_1p
cdef double C_ombh2_m_1p
cdef double[:, :] a_data_ombh2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_1p, C_ombh2_m_1p, a_data_ombh2_m_1p, b_data_ombh2_m_1p, c_data_ombh2_m_1p, lps_cc_data_ombh2_m_1p, lps_cs_data_ombh2_m_1p, lps_ss_data_ombh2_m_1p, scale_factor_data_ombh2_m_1p, window_c_data_ombh2_m_1p, window_s_data_ombh2_m_1p, mps_data_ombh2_m_1p, z_at_chi_data_ombh2_m_1p, cmbps_ombh2_m_1p = data_import_func('data_ombh2_m_1p')

cdef double lbs_ombh2_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_1p, a_data_ombh2_m_1p, b_data_ombh2_m_1p, c_data_ombh2_m_1p, scale_factor_data_ombh2_m_1p, window_c_data_ombh2_m_1p, window_s_data_ombh2_m_1p, mps_data_ombh2_m_1p, z_at_chi_data_ombh2_m_1p)

cdef double lps_ombh2_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_1p, lps_cs_data_ombh2_m_1p, lps_ss_data_ombh2_m_1p, cmbps_ombh2_m_1p)
            

cdef double[:] cosm_par_ombh2_m_2p
cdef double C_ombh2_m_2p
cdef double[:, :] a_data_ombh2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_2p, C_ombh2_m_2p, a_data_ombh2_m_2p, b_data_ombh2_m_2p, c_data_ombh2_m_2p, lps_cc_data_ombh2_m_2p, lps_cs_data_ombh2_m_2p, lps_ss_data_ombh2_m_2p, scale_factor_data_ombh2_m_2p, window_c_data_ombh2_m_2p, window_s_data_ombh2_m_2p, mps_data_ombh2_m_2p, z_at_chi_data_ombh2_m_2p, cmbps_ombh2_m_2p = data_import_func('data_ombh2_m_2p')

cdef double lbs_ombh2_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_2p, a_data_ombh2_m_2p, b_data_ombh2_m_2p, c_data_ombh2_m_2p, scale_factor_data_ombh2_m_2p, window_c_data_ombh2_m_2p, window_s_data_ombh2_m_2p, mps_data_ombh2_m_2p, z_at_chi_data_ombh2_m_2p)

cdef double lps_ombh2_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_2p, lps_cs_data_ombh2_m_2p, lps_ss_data_ombh2_m_2p, cmbps_ombh2_m_2p)
            

cdef double[:] cosm_par_omch2_p_2m
cdef double C_omch2_p_2m
cdef double[:, :] a_data_omch2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_2m, C_omch2_p_2m, a_data_omch2_p_2m, b_data_omch2_p_2m, c_data_omch2_p_2m, lps_cc_data_omch2_p_2m, lps_cs_data_omch2_p_2m, lps_ss_data_omch2_p_2m, scale_factor_data_omch2_p_2m, window_c_data_omch2_p_2m, window_s_data_omch2_p_2m, mps_data_omch2_p_2m, z_at_chi_data_omch2_p_2m, cmbps_omch2_p_2m = data_import_func('data_omch2_p_2m')

cdef double lbs_omch2_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_2m, a_data_omch2_p_2m, b_data_omch2_p_2m, c_data_omch2_p_2m, scale_factor_data_omch2_p_2m, window_c_data_omch2_p_2m, window_s_data_omch2_p_2m, mps_data_omch2_p_2m, z_at_chi_data_omch2_p_2m)

cdef double lps_omch2_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_2m, lps_cs_data_omch2_p_2m, lps_ss_data_omch2_p_2m, cmbps_omch2_p_2m)
            

cdef double[:] cosm_par_omch2_p_1m
cdef double C_omch2_p_1m
cdef double[:, :] a_data_omch2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_1m, C_omch2_p_1m, a_data_omch2_p_1m, b_data_omch2_p_1m, c_data_omch2_p_1m, lps_cc_data_omch2_p_1m, lps_cs_data_omch2_p_1m, lps_ss_data_omch2_p_1m, scale_factor_data_omch2_p_1m, window_c_data_omch2_p_1m, window_s_data_omch2_p_1m, mps_data_omch2_p_1m, z_at_chi_data_omch2_p_1m, cmbps_omch2_p_1m = data_import_func('data_omch2_p_1m')

cdef double lbs_omch2_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_1m, a_data_omch2_p_1m, b_data_omch2_p_1m, c_data_omch2_p_1m, scale_factor_data_omch2_p_1m, window_c_data_omch2_p_1m, window_s_data_omch2_p_1m, mps_data_omch2_p_1m, z_at_chi_data_omch2_p_1m)

cdef double lps_omch2_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_1m, lps_cs_data_omch2_p_1m, lps_ss_data_omch2_p_1m, cmbps_omch2_p_1m)
            

cdef double[:] cosm_par_omch2_p_0
cdef double C_omch2_p_0
cdef double[:, :] a_data_omch2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_0, C_omch2_p_0, a_data_omch2_p_0, b_data_omch2_p_0, c_data_omch2_p_0, lps_cc_data_omch2_p_0, lps_cs_data_omch2_p_0, lps_ss_data_omch2_p_0, scale_factor_data_omch2_p_0, window_c_data_omch2_p_0, window_s_data_omch2_p_0, mps_data_omch2_p_0, z_at_chi_data_omch2_p_0, cmbps_omch2_p_0 = data_import_func('data_omch2_p_0')

cdef double lbs_omch2_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_0, a_data_omch2_p_0, b_data_omch2_p_0, c_data_omch2_p_0, scale_factor_data_omch2_p_0, window_c_data_omch2_p_0, window_s_data_omch2_p_0, mps_data_omch2_p_0, z_at_chi_data_omch2_p_0)

cdef double lps_omch2_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_0, lps_cs_data_omch2_p_0, lps_ss_data_omch2_p_0, cmbps_omch2_p_0)
            

cdef double[:] cosm_par_omch2_p_1p
cdef double C_omch2_p_1p
cdef double[:, :] a_data_omch2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_1p, C_omch2_p_1p, a_data_omch2_p_1p, b_data_omch2_p_1p, c_data_omch2_p_1p, lps_cc_data_omch2_p_1p, lps_cs_data_omch2_p_1p, lps_ss_data_omch2_p_1p, scale_factor_data_omch2_p_1p, window_c_data_omch2_p_1p, window_s_data_omch2_p_1p, mps_data_omch2_p_1p, z_at_chi_data_omch2_p_1p, cmbps_omch2_p_1p = data_import_func('data_omch2_p_1p')

cdef double lbs_omch2_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_1p, a_data_omch2_p_1p, b_data_omch2_p_1p, c_data_omch2_p_1p, scale_factor_data_omch2_p_1p, window_c_data_omch2_p_1p, window_s_data_omch2_p_1p, mps_data_omch2_p_1p, z_at_chi_data_omch2_p_1p)

cdef double lps_omch2_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_1p, lps_cs_data_omch2_p_1p, lps_ss_data_omch2_p_1p, cmbps_omch2_p_1p)
            

cdef double[:] cosm_par_omch2_p_2p
cdef double C_omch2_p_2p
cdef double[:, :] a_data_omch2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_2p, C_omch2_p_2p, a_data_omch2_p_2p, b_data_omch2_p_2p, c_data_omch2_p_2p, lps_cc_data_omch2_p_2p, lps_cs_data_omch2_p_2p, lps_ss_data_omch2_p_2p, scale_factor_data_omch2_p_2p, window_c_data_omch2_p_2p, window_s_data_omch2_p_2p, mps_data_omch2_p_2p, z_at_chi_data_omch2_p_2p, cmbps_omch2_p_2p = data_import_func('data_omch2_p_2p')

cdef double lbs_omch2_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_2p, a_data_omch2_p_2p, b_data_omch2_p_2p, c_data_omch2_p_2p, scale_factor_data_omch2_p_2p, window_c_data_omch2_p_2p, window_s_data_omch2_p_2p, mps_data_omch2_p_2p, z_at_chi_data_omch2_p_2p)

cdef double lps_omch2_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_2p, lps_cs_data_omch2_p_2p, lps_ss_data_omch2_p_2p, cmbps_omch2_p_2p)
            

cdef double[:] cosm_par_omch2_m_2m
cdef double C_omch2_m_2m
cdef double[:, :] a_data_omch2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_2m, C_omch2_m_2m, a_data_omch2_m_2m, b_data_omch2_m_2m, c_data_omch2_m_2m, lps_cc_data_omch2_m_2m, lps_cs_data_omch2_m_2m, lps_ss_data_omch2_m_2m, scale_factor_data_omch2_m_2m, window_c_data_omch2_m_2m, window_s_data_omch2_m_2m, mps_data_omch2_m_2m, z_at_chi_data_omch2_m_2m, cmbps_omch2_m_2m = data_import_func('data_omch2_m_2m')

cdef double lbs_omch2_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_2m, a_data_omch2_m_2m, b_data_omch2_m_2m, c_data_omch2_m_2m, scale_factor_data_omch2_m_2m, window_c_data_omch2_m_2m, window_s_data_omch2_m_2m, mps_data_omch2_m_2m, z_at_chi_data_omch2_m_2m)

cdef double lps_omch2_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_2m, lps_cs_data_omch2_m_2m, lps_ss_data_omch2_m_2m, cmbps_omch2_m_2m)
            

cdef double[:] cosm_par_omch2_m_1m
cdef double C_omch2_m_1m
cdef double[:, :] a_data_omch2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_1m, C_omch2_m_1m, a_data_omch2_m_1m, b_data_omch2_m_1m, c_data_omch2_m_1m, lps_cc_data_omch2_m_1m, lps_cs_data_omch2_m_1m, lps_ss_data_omch2_m_1m, scale_factor_data_omch2_m_1m, window_c_data_omch2_m_1m, window_s_data_omch2_m_1m, mps_data_omch2_m_1m, z_at_chi_data_omch2_m_1m, cmbps_omch2_m_1m = data_import_func('data_omch2_m_1m')

cdef double lbs_omch2_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_1m, a_data_omch2_m_1m, b_data_omch2_m_1m, c_data_omch2_m_1m, scale_factor_data_omch2_m_1m, window_c_data_omch2_m_1m, window_s_data_omch2_m_1m, mps_data_omch2_m_1m, z_at_chi_data_omch2_m_1m)

cdef double lps_omch2_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_1m, lps_cs_data_omch2_m_1m, lps_ss_data_omch2_m_1m, cmbps_omch2_m_1m)
            

cdef double[:] cosm_par_omch2_m_0
cdef double C_omch2_m_0
cdef double[:, :] a_data_omch2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_0, C_omch2_m_0, a_data_omch2_m_0, b_data_omch2_m_0, c_data_omch2_m_0, lps_cc_data_omch2_m_0, lps_cs_data_omch2_m_0, lps_ss_data_omch2_m_0, scale_factor_data_omch2_m_0, window_c_data_omch2_m_0, window_s_data_omch2_m_0, mps_data_omch2_m_0, z_at_chi_data_omch2_m_0, cmbps_omch2_m_0 = data_import_func('data_omch2_m_0')

cdef double lbs_omch2_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_0, a_data_omch2_m_0, b_data_omch2_m_0, c_data_omch2_m_0, scale_factor_data_omch2_m_0, window_c_data_omch2_m_0, window_s_data_omch2_m_0, mps_data_omch2_m_0, z_at_chi_data_omch2_m_0)

cdef double lps_omch2_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_0, lps_cs_data_omch2_m_0, lps_ss_data_omch2_m_0, cmbps_omch2_m_0)
            

cdef double[:] cosm_par_omch2_m_1p
cdef double C_omch2_m_1p
cdef double[:, :] a_data_omch2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_1p, C_omch2_m_1p, a_data_omch2_m_1p, b_data_omch2_m_1p, c_data_omch2_m_1p, lps_cc_data_omch2_m_1p, lps_cs_data_omch2_m_1p, lps_ss_data_omch2_m_1p, scale_factor_data_omch2_m_1p, window_c_data_omch2_m_1p, window_s_data_omch2_m_1p, mps_data_omch2_m_1p, z_at_chi_data_omch2_m_1p, cmbps_omch2_m_1p = data_import_func('data_omch2_m_1p')

cdef double lbs_omch2_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_1p, a_data_omch2_m_1p, b_data_omch2_m_1p, c_data_omch2_m_1p, scale_factor_data_omch2_m_1p, window_c_data_omch2_m_1p, window_s_data_omch2_m_1p, mps_data_omch2_m_1p, z_at_chi_data_omch2_m_1p)

cdef double lps_omch2_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_1p, lps_cs_data_omch2_m_1p, lps_ss_data_omch2_m_1p, cmbps_omch2_m_1p)
            

cdef double[:] cosm_par_omch2_m_2p
cdef double C_omch2_m_2p
cdef double[:, :] a_data_omch2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_2p, C_omch2_m_2p, a_data_omch2_m_2p, b_data_omch2_m_2p, c_data_omch2_m_2p, lps_cc_data_omch2_m_2p, lps_cs_data_omch2_m_2p, lps_ss_data_omch2_m_2p, scale_factor_data_omch2_m_2p, window_c_data_omch2_m_2p, window_s_data_omch2_m_2p, mps_data_omch2_m_2p, z_at_chi_data_omch2_m_2p, cmbps_omch2_m_2p = data_import_func('data_omch2_m_2p')

cdef double lbs_omch2_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_2p, a_data_omch2_m_2p, b_data_omch2_m_2p, c_data_omch2_m_2p, scale_factor_data_omch2_m_2p, window_c_data_omch2_m_2p, window_s_data_omch2_m_2p, mps_data_omch2_m_2p, z_at_chi_data_omch2_m_2p)

cdef double lps_omch2_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_2p, lps_cs_data_omch2_m_2p, lps_ss_data_omch2_m_2p, cmbps_omch2_m_2p)
            

cdef double[:] cosm_par_ns_p_2m
cdef double C_ns_p_2m
cdef double[:, :] a_data_ns_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_2m, C_ns_p_2m, a_data_ns_p_2m, b_data_ns_p_2m, c_data_ns_p_2m, lps_cc_data_ns_p_2m, lps_cs_data_ns_p_2m, lps_ss_data_ns_p_2m, scale_factor_data_ns_p_2m, window_c_data_ns_p_2m, window_s_data_ns_p_2m, mps_data_ns_p_2m, z_at_chi_data_ns_p_2m, cmbps_ns_p_2m = data_import_func('data_ns_p_2m')

cdef double lbs_ns_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_2m, a_data_ns_p_2m, b_data_ns_p_2m, c_data_ns_p_2m, scale_factor_data_ns_p_2m, window_c_data_ns_p_2m, window_s_data_ns_p_2m, mps_data_ns_p_2m, z_at_chi_data_ns_p_2m)

cdef double lps_ns_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_2m, lps_cs_data_ns_p_2m, lps_ss_data_ns_p_2m, cmbps_ns_p_2m)
            

cdef double[:] cosm_par_ns_p_1m
cdef double C_ns_p_1m
cdef double[:, :] a_data_ns_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_1m, C_ns_p_1m, a_data_ns_p_1m, b_data_ns_p_1m, c_data_ns_p_1m, lps_cc_data_ns_p_1m, lps_cs_data_ns_p_1m, lps_ss_data_ns_p_1m, scale_factor_data_ns_p_1m, window_c_data_ns_p_1m, window_s_data_ns_p_1m, mps_data_ns_p_1m, z_at_chi_data_ns_p_1m, cmbps_ns_p_1m = data_import_func('data_ns_p_1m')

cdef double lbs_ns_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_1m, a_data_ns_p_1m, b_data_ns_p_1m, c_data_ns_p_1m, scale_factor_data_ns_p_1m, window_c_data_ns_p_1m, window_s_data_ns_p_1m, mps_data_ns_p_1m, z_at_chi_data_ns_p_1m)

cdef double lps_ns_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_1m, lps_cs_data_ns_p_1m, lps_ss_data_ns_p_1m, cmbps_ns_p_1m)
            

cdef double[:] cosm_par_ns_p_0
cdef double C_ns_p_0
cdef double[:, :] a_data_ns_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_0, C_ns_p_0, a_data_ns_p_0, b_data_ns_p_0, c_data_ns_p_0, lps_cc_data_ns_p_0, lps_cs_data_ns_p_0, lps_ss_data_ns_p_0, scale_factor_data_ns_p_0, window_c_data_ns_p_0, window_s_data_ns_p_0, mps_data_ns_p_0, z_at_chi_data_ns_p_0, cmbps_ns_p_0 = data_import_func('data_ns_p_0')

cdef double lbs_ns_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_0, a_data_ns_p_0, b_data_ns_p_0, c_data_ns_p_0, scale_factor_data_ns_p_0, window_c_data_ns_p_0, window_s_data_ns_p_0, mps_data_ns_p_0, z_at_chi_data_ns_p_0)

cdef double lps_ns_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_0, lps_cs_data_ns_p_0, lps_ss_data_ns_p_0, cmbps_ns_p_0)
            

cdef double[:] cosm_par_ns_p_1p
cdef double C_ns_p_1p
cdef double[:, :] a_data_ns_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_1p, C_ns_p_1p, a_data_ns_p_1p, b_data_ns_p_1p, c_data_ns_p_1p, lps_cc_data_ns_p_1p, lps_cs_data_ns_p_1p, lps_ss_data_ns_p_1p, scale_factor_data_ns_p_1p, window_c_data_ns_p_1p, window_s_data_ns_p_1p, mps_data_ns_p_1p, z_at_chi_data_ns_p_1p, cmbps_ns_p_1p = data_import_func('data_ns_p_1p')

cdef double lbs_ns_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_1p, a_data_ns_p_1p, b_data_ns_p_1p, c_data_ns_p_1p, scale_factor_data_ns_p_1p, window_c_data_ns_p_1p, window_s_data_ns_p_1p, mps_data_ns_p_1p, z_at_chi_data_ns_p_1p)

cdef double lps_ns_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_1p, lps_cs_data_ns_p_1p, lps_ss_data_ns_p_1p, cmbps_ns_p_1p)
            

cdef double[:] cosm_par_ns_p_2p
cdef double C_ns_p_2p
cdef double[:, :] a_data_ns_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_2p, C_ns_p_2p, a_data_ns_p_2p, b_data_ns_p_2p, c_data_ns_p_2p, lps_cc_data_ns_p_2p, lps_cs_data_ns_p_2p, lps_ss_data_ns_p_2p, scale_factor_data_ns_p_2p, window_c_data_ns_p_2p, window_s_data_ns_p_2p, mps_data_ns_p_2p, z_at_chi_data_ns_p_2p, cmbps_ns_p_2p = data_import_func('data_ns_p_2p')

cdef double lbs_ns_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_2p, a_data_ns_p_2p, b_data_ns_p_2p, c_data_ns_p_2p, scale_factor_data_ns_p_2p, window_c_data_ns_p_2p, window_s_data_ns_p_2p, mps_data_ns_p_2p, z_at_chi_data_ns_p_2p)

cdef double lps_ns_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_2p, lps_cs_data_ns_p_2p, lps_ss_data_ns_p_2p, cmbps_ns_p_2p)
            

cdef double[:] cosm_par_ns_m_2m
cdef double C_ns_m_2m
cdef double[:, :] a_data_ns_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_2m, C_ns_m_2m, a_data_ns_m_2m, b_data_ns_m_2m, c_data_ns_m_2m, lps_cc_data_ns_m_2m, lps_cs_data_ns_m_2m, lps_ss_data_ns_m_2m, scale_factor_data_ns_m_2m, window_c_data_ns_m_2m, window_s_data_ns_m_2m, mps_data_ns_m_2m, z_at_chi_data_ns_m_2m, cmbps_ns_m_2m = data_import_func('data_ns_m_2m')

cdef double lbs_ns_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_2m, a_data_ns_m_2m, b_data_ns_m_2m, c_data_ns_m_2m, scale_factor_data_ns_m_2m, window_c_data_ns_m_2m, window_s_data_ns_m_2m, mps_data_ns_m_2m, z_at_chi_data_ns_m_2m)

cdef double lps_ns_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_2m, lps_cs_data_ns_m_2m, lps_ss_data_ns_m_2m, cmbps_ns_m_2m)
            

cdef double[:] cosm_par_ns_m_1m
cdef double C_ns_m_1m
cdef double[:, :] a_data_ns_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_1m, C_ns_m_1m, a_data_ns_m_1m, b_data_ns_m_1m, c_data_ns_m_1m, lps_cc_data_ns_m_1m, lps_cs_data_ns_m_1m, lps_ss_data_ns_m_1m, scale_factor_data_ns_m_1m, window_c_data_ns_m_1m, window_s_data_ns_m_1m, mps_data_ns_m_1m, z_at_chi_data_ns_m_1m, cmbps_ns_m_1m = data_import_func('data_ns_m_1m')

cdef double lbs_ns_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_1m, a_data_ns_m_1m, b_data_ns_m_1m, c_data_ns_m_1m, scale_factor_data_ns_m_1m, window_c_data_ns_m_1m, window_s_data_ns_m_1m, mps_data_ns_m_1m, z_at_chi_data_ns_m_1m)

cdef double lps_ns_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_1m, lps_cs_data_ns_m_1m, lps_ss_data_ns_m_1m, cmbps_ns_m_1m)
            

cdef double[:] cosm_par_ns_m_0
cdef double C_ns_m_0
cdef double[:, :] a_data_ns_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_0, C_ns_m_0, a_data_ns_m_0, b_data_ns_m_0, c_data_ns_m_0, lps_cc_data_ns_m_0, lps_cs_data_ns_m_0, lps_ss_data_ns_m_0, scale_factor_data_ns_m_0, window_c_data_ns_m_0, window_s_data_ns_m_0, mps_data_ns_m_0, z_at_chi_data_ns_m_0, cmbps_ns_m_0 = data_import_func('data_ns_m_0')

cdef double lbs_ns_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_0, a_data_ns_m_0, b_data_ns_m_0, c_data_ns_m_0, scale_factor_data_ns_m_0, window_c_data_ns_m_0, window_s_data_ns_m_0, mps_data_ns_m_0, z_at_chi_data_ns_m_0)

cdef double lps_ns_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_0, lps_cs_data_ns_m_0, lps_ss_data_ns_m_0, cmbps_ns_m_0)
            

cdef double[:] cosm_par_ns_m_1p
cdef double C_ns_m_1p
cdef double[:, :] a_data_ns_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_1p, C_ns_m_1p, a_data_ns_m_1p, b_data_ns_m_1p, c_data_ns_m_1p, lps_cc_data_ns_m_1p, lps_cs_data_ns_m_1p, lps_ss_data_ns_m_1p, scale_factor_data_ns_m_1p, window_c_data_ns_m_1p, window_s_data_ns_m_1p, mps_data_ns_m_1p, z_at_chi_data_ns_m_1p, cmbps_ns_m_1p = data_import_func('data_ns_m_1p')

cdef double lbs_ns_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_1p, a_data_ns_m_1p, b_data_ns_m_1p, c_data_ns_m_1p, scale_factor_data_ns_m_1p, window_c_data_ns_m_1p, window_s_data_ns_m_1p, mps_data_ns_m_1p, z_at_chi_data_ns_m_1p)

cdef double lps_ns_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_1p, lps_cs_data_ns_m_1p, lps_ss_data_ns_m_1p, cmbps_ns_m_1p)
            

cdef double[:] cosm_par_ns_m_2p
cdef double C_ns_m_2p
cdef double[:, :] a_data_ns_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_2p, C_ns_m_2p, a_data_ns_m_2p, b_data_ns_m_2p, c_data_ns_m_2p, lps_cc_data_ns_m_2p, lps_cs_data_ns_m_2p, lps_ss_data_ns_m_2p, scale_factor_data_ns_m_2p, window_c_data_ns_m_2p, window_s_data_ns_m_2p, mps_data_ns_m_2p, z_at_chi_data_ns_m_2p, cmbps_ns_m_2p = data_import_func('data_ns_m_2p')

cdef double lbs_ns_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_2p, a_data_ns_m_2p, b_data_ns_m_2p, c_data_ns_m_2p, scale_factor_data_ns_m_2p, window_c_data_ns_m_2p, window_s_data_ns_m_2p, mps_data_ns_m_2p, z_at_chi_data_ns_m_2p)

cdef double lps_ns_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_2p, lps_cs_data_ns_m_2p, lps_ss_data_ns_m_2p, cmbps_ns_m_2p)
            

cdef double[:] cosm_par_As_p_2m
cdef double C_As_p_2m
cdef double[:, :] a_data_As_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_2m, C_As_p_2m, a_data_As_p_2m, b_data_As_p_2m, c_data_As_p_2m, lps_cc_data_As_p_2m, lps_cs_data_As_p_2m, lps_ss_data_As_p_2m, scale_factor_data_As_p_2m, window_c_data_As_p_2m, window_s_data_As_p_2m, mps_data_As_p_2m, z_at_chi_data_As_p_2m, cmbps_As_p_2m = data_import_func('data_As_p_2m')

cdef double lbs_As_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_2m, a_data_As_p_2m, b_data_As_p_2m, c_data_As_p_2m, scale_factor_data_As_p_2m, window_c_data_As_p_2m, window_s_data_As_p_2m, mps_data_As_p_2m, z_at_chi_data_As_p_2m)

cdef double lps_As_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_2m, lps_cs_data_As_p_2m, lps_ss_data_As_p_2m, cmbps_As_p_2m)
            

cdef double[:] cosm_par_As_p_1m
cdef double C_As_p_1m
cdef double[:, :] a_data_As_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_1m, C_As_p_1m, a_data_As_p_1m, b_data_As_p_1m, c_data_As_p_1m, lps_cc_data_As_p_1m, lps_cs_data_As_p_1m, lps_ss_data_As_p_1m, scale_factor_data_As_p_1m, window_c_data_As_p_1m, window_s_data_As_p_1m, mps_data_As_p_1m, z_at_chi_data_As_p_1m, cmbps_As_p_1m = data_import_func('data_As_p_1m')

cdef double lbs_As_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_1m, a_data_As_p_1m, b_data_As_p_1m, c_data_As_p_1m, scale_factor_data_As_p_1m, window_c_data_As_p_1m, window_s_data_As_p_1m, mps_data_As_p_1m, z_at_chi_data_As_p_1m)

cdef double lps_As_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_1m, lps_cs_data_As_p_1m, lps_ss_data_As_p_1m, cmbps_As_p_1m)
            

cdef double[:] cosm_par_As_p_0
cdef double C_As_p_0
cdef double[:, :] a_data_As_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_0, C_As_p_0, a_data_As_p_0, b_data_As_p_0, c_data_As_p_0, lps_cc_data_As_p_0, lps_cs_data_As_p_0, lps_ss_data_As_p_0, scale_factor_data_As_p_0, window_c_data_As_p_0, window_s_data_As_p_0, mps_data_As_p_0, z_at_chi_data_As_p_0, cmbps_As_p_0 = data_import_func('data_As_p_0')

cdef double lbs_As_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_0, a_data_As_p_0, b_data_As_p_0, c_data_As_p_0, scale_factor_data_As_p_0, window_c_data_As_p_0, window_s_data_As_p_0, mps_data_As_p_0, z_at_chi_data_As_p_0)

cdef double lps_As_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_0, lps_cs_data_As_p_0, lps_ss_data_As_p_0, cmbps_As_p_0)
            

cdef double[:] cosm_par_As_p_1p
cdef double C_As_p_1p
cdef double[:, :] a_data_As_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_1p, C_As_p_1p, a_data_As_p_1p, b_data_As_p_1p, c_data_As_p_1p, lps_cc_data_As_p_1p, lps_cs_data_As_p_1p, lps_ss_data_As_p_1p, scale_factor_data_As_p_1p, window_c_data_As_p_1p, window_s_data_As_p_1p, mps_data_As_p_1p, z_at_chi_data_As_p_1p, cmbps_As_p_1p = data_import_func('data_As_p_1p')

cdef double lbs_As_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_1p, a_data_As_p_1p, b_data_As_p_1p, c_data_As_p_1p, scale_factor_data_As_p_1p, window_c_data_As_p_1p, window_s_data_As_p_1p, mps_data_As_p_1p, z_at_chi_data_As_p_1p)

cdef double lps_As_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_1p, lps_cs_data_As_p_1p, lps_ss_data_As_p_1p, cmbps_As_p_1p)
            

cdef double[:] cosm_par_As_p_2p
cdef double C_As_p_2p
cdef double[:, :] a_data_As_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_2p, C_As_p_2p, a_data_As_p_2p, b_data_As_p_2p, c_data_As_p_2p, lps_cc_data_As_p_2p, lps_cs_data_As_p_2p, lps_ss_data_As_p_2p, scale_factor_data_As_p_2p, window_c_data_As_p_2p, window_s_data_As_p_2p, mps_data_As_p_2p, z_at_chi_data_As_p_2p, cmbps_As_p_2p = data_import_func('data_As_p_2p')

cdef double lbs_As_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_2p, a_data_As_p_2p, b_data_As_p_2p, c_data_As_p_2p, scale_factor_data_As_p_2p, window_c_data_As_p_2p, window_s_data_As_p_2p, mps_data_As_p_2p, z_at_chi_data_As_p_2p)

cdef double lps_As_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_2p, lps_cs_data_As_p_2p, lps_ss_data_As_p_2p, cmbps_As_p_2p)
            

cdef double[:] cosm_par_As_m_2m
cdef double C_As_m_2m
cdef double[:, :] a_data_As_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_2m, C_As_m_2m, a_data_As_m_2m, b_data_As_m_2m, c_data_As_m_2m, lps_cc_data_As_m_2m, lps_cs_data_As_m_2m, lps_ss_data_As_m_2m, scale_factor_data_As_m_2m, window_c_data_As_m_2m, window_s_data_As_m_2m, mps_data_As_m_2m, z_at_chi_data_As_m_2m, cmbps_As_m_2m = data_import_func('data_As_m_2m')

cdef double lbs_As_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_2m, a_data_As_m_2m, b_data_As_m_2m, c_data_As_m_2m, scale_factor_data_As_m_2m, window_c_data_As_m_2m, window_s_data_As_m_2m, mps_data_As_m_2m, z_at_chi_data_As_m_2m)

cdef double lps_As_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_2m, lps_cs_data_As_m_2m, lps_ss_data_As_m_2m, cmbps_As_m_2m)
            

cdef double[:] cosm_par_As_m_1m
cdef double C_As_m_1m
cdef double[:, :] a_data_As_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_1m, C_As_m_1m, a_data_As_m_1m, b_data_As_m_1m, c_data_As_m_1m, lps_cc_data_As_m_1m, lps_cs_data_As_m_1m, lps_ss_data_As_m_1m, scale_factor_data_As_m_1m, window_c_data_As_m_1m, window_s_data_As_m_1m, mps_data_As_m_1m, z_at_chi_data_As_m_1m, cmbps_As_m_1m = data_import_func('data_As_m_1m')

cdef double lbs_As_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_1m, a_data_As_m_1m, b_data_As_m_1m, c_data_As_m_1m, scale_factor_data_As_m_1m, window_c_data_As_m_1m, window_s_data_As_m_1m, mps_data_As_m_1m, z_at_chi_data_As_m_1m)

cdef double lps_As_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_1m, lps_cs_data_As_m_1m, lps_ss_data_As_m_1m, cmbps_As_m_1m)
            

cdef double[:] cosm_par_As_m_0
cdef double C_As_m_0
cdef double[:, :] a_data_As_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_0, C_As_m_0, a_data_As_m_0, b_data_As_m_0, c_data_As_m_0, lps_cc_data_As_m_0, lps_cs_data_As_m_0, lps_ss_data_As_m_0, scale_factor_data_As_m_0, window_c_data_As_m_0, window_s_data_As_m_0, mps_data_As_m_0, z_at_chi_data_As_m_0, cmbps_As_m_0 = data_import_func('data_As_m_0')

cdef double lbs_As_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_0, a_data_As_m_0, b_data_As_m_0, c_data_As_m_0, scale_factor_data_As_m_0, window_c_data_As_m_0, window_s_data_As_m_0, mps_data_As_m_0, z_at_chi_data_As_m_0)

cdef double lps_As_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_0, lps_cs_data_As_m_0, lps_ss_data_As_m_0, cmbps_As_m_0)
            

cdef double[:] cosm_par_As_m_1p
cdef double C_As_m_1p
cdef double[:, :] a_data_As_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_1p, C_As_m_1p, a_data_As_m_1p, b_data_As_m_1p, c_data_As_m_1p, lps_cc_data_As_m_1p, lps_cs_data_As_m_1p, lps_ss_data_As_m_1p, scale_factor_data_As_m_1p, window_c_data_As_m_1p, window_s_data_As_m_1p, mps_data_As_m_1p, z_at_chi_data_As_m_1p, cmbps_As_m_1p = data_import_func('data_As_m_1p')

cdef double lbs_As_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_1p, a_data_As_m_1p, b_data_As_m_1p, c_data_As_m_1p, scale_factor_data_As_m_1p, window_c_data_As_m_1p, window_s_data_As_m_1p, mps_data_As_m_1p, z_at_chi_data_As_m_1p)

cdef double lps_As_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_1p, lps_cs_data_As_m_1p, lps_ss_data_As_m_1p, cmbps_As_m_1p)
            

cdef double[:] cosm_par_As_m_2p
cdef double C_As_m_2p
cdef double[:, :] a_data_As_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_2p, C_As_m_2p, a_data_As_m_2p, b_data_As_m_2p, c_data_As_m_2p, lps_cc_data_As_m_2p, lps_cs_data_As_m_2p, lps_ss_data_As_m_2p, scale_factor_data_As_m_2p, window_c_data_As_m_2p, window_s_data_As_m_2p, mps_data_As_m_2p, z_at_chi_data_As_m_2p, cmbps_As_m_2p = data_import_func('data_As_m_2p')

cdef double lbs_As_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_2p, a_data_As_m_2p, b_data_As_m_2p, c_data_As_m_2p, scale_factor_data_As_m_2p, window_c_data_As_m_2p, window_s_data_As_m_2p, mps_data_As_m_2p, z_at_chi_data_As_m_2p)

cdef double lps_As_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_2p, lps_cs_data_As_m_2p, lps_ss_data_As_m_2p, cmbps_As_m_2p)
            

cdef double[:] cosm_par_mnu_p_2m
cdef double C_mnu_p_2m
cdef double[:, :] a_data_mnu_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_2m, C_mnu_p_2m, a_data_mnu_p_2m, b_data_mnu_p_2m, c_data_mnu_p_2m, lps_cc_data_mnu_p_2m, lps_cs_data_mnu_p_2m, lps_ss_data_mnu_p_2m, scale_factor_data_mnu_p_2m, window_c_data_mnu_p_2m, window_s_data_mnu_p_2m, mps_data_mnu_p_2m, z_at_chi_data_mnu_p_2m, cmbps_mnu_p_2m = data_import_func('data_mnu_p_2m')

cdef double lbs_mnu_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_2m, a_data_mnu_p_2m, b_data_mnu_p_2m, c_data_mnu_p_2m, scale_factor_data_mnu_p_2m, window_c_data_mnu_p_2m, window_s_data_mnu_p_2m, mps_data_mnu_p_2m, z_at_chi_data_mnu_p_2m)

cdef double lps_mnu_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_2m, lps_cs_data_mnu_p_2m, lps_ss_data_mnu_p_2m, cmbps_mnu_p_2m)
            

cdef double[:] cosm_par_mnu_p_1m
cdef double C_mnu_p_1m
cdef double[:, :] a_data_mnu_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_1m, C_mnu_p_1m, a_data_mnu_p_1m, b_data_mnu_p_1m, c_data_mnu_p_1m, lps_cc_data_mnu_p_1m, lps_cs_data_mnu_p_1m, lps_ss_data_mnu_p_1m, scale_factor_data_mnu_p_1m, window_c_data_mnu_p_1m, window_s_data_mnu_p_1m, mps_data_mnu_p_1m, z_at_chi_data_mnu_p_1m, cmbps_mnu_p_1m = data_import_func('data_mnu_p_1m')

cdef double lbs_mnu_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_1m, a_data_mnu_p_1m, b_data_mnu_p_1m, c_data_mnu_p_1m, scale_factor_data_mnu_p_1m, window_c_data_mnu_p_1m, window_s_data_mnu_p_1m, mps_data_mnu_p_1m, z_at_chi_data_mnu_p_1m)

cdef double lps_mnu_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_1m, lps_cs_data_mnu_p_1m, lps_ss_data_mnu_p_1m, cmbps_mnu_p_1m)
            

cdef double[:] cosm_par_mnu_p_0
cdef double C_mnu_p_0
cdef double[:, :] a_data_mnu_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_0, C_mnu_p_0, a_data_mnu_p_0, b_data_mnu_p_0, c_data_mnu_p_0, lps_cc_data_mnu_p_0, lps_cs_data_mnu_p_0, lps_ss_data_mnu_p_0, scale_factor_data_mnu_p_0, window_c_data_mnu_p_0, window_s_data_mnu_p_0, mps_data_mnu_p_0, z_at_chi_data_mnu_p_0, cmbps_mnu_p_0 = data_import_func('data_mnu_p_0')

cdef double lbs_mnu_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_0, a_data_mnu_p_0, b_data_mnu_p_0, c_data_mnu_p_0, scale_factor_data_mnu_p_0, window_c_data_mnu_p_0, window_s_data_mnu_p_0, mps_data_mnu_p_0, z_at_chi_data_mnu_p_0)

cdef double lps_mnu_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_0, lps_cs_data_mnu_p_0, lps_ss_data_mnu_p_0, cmbps_mnu_p_0)
            

cdef double[:] cosm_par_mnu_p_1p
cdef double C_mnu_p_1p
cdef double[:, :] a_data_mnu_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_1p, C_mnu_p_1p, a_data_mnu_p_1p, b_data_mnu_p_1p, c_data_mnu_p_1p, lps_cc_data_mnu_p_1p, lps_cs_data_mnu_p_1p, lps_ss_data_mnu_p_1p, scale_factor_data_mnu_p_1p, window_c_data_mnu_p_1p, window_s_data_mnu_p_1p, mps_data_mnu_p_1p, z_at_chi_data_mnu_p_1p, cmbps_mnu_p_1p = data_import_func('data_mnu_p_1p')

cdef double lbs_mnu_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_1p, a_data_mnu_p_1p, b_data_mnu_p_1p, c_data_mnu_p_1p, scale_factor_data_mnu_p_1p, window_c_data_mnu_p_1p, window_s_data_mnu_p_1p, mps_data_mnu_p_1p, z_at_chi_data_mnu_p_1p)

cdef double lps_mnu_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_1p, lps_cs_data_mnu_p_1p, lps_ss_data_mnu_p_1p, cmbps_mnu_p_1p)
            

cdef double[:] cosm_par_mnu_p_2p
cdef double C_mnu_p_2p
cdef double[:, :] a_data_mnu_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_2p, C_mnu_p_2p, a_data_mnu_p_2p, b_data_mnu_p_2p, c_data_mnu_p_2p, lps_cc_data_mnu_p_2p, lps_cs_data_mnu_p_2p, lps_ss_data_mnu_p_2p, scale_factor_data_mnu_p_2p, window_c_data_mnu_p_2p, window_s_data_mnu_p_2p, mps_data_mnu_p_2p, z_at_chi_data_mnu_p_2p, cmbps_mnu_p_2p = data_import_func('data_mnu_p_2p')

cdef double lbs_mnu_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_2p, a_data_mnu_p_2p, b_data_mnu_p_2p, c_data_mnu_p_2p, scale_factor_data_mnu_p_2p, window_c_data_mnu_p_2p, window_s_data_mnu_p_2p, mps_data_mnu_p_2p, z_at_chi_data_mnu_p_2p)

cdef double lps_mnu_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_2p, lps_cs_data_mnu_p_2p, lps_ss_data_mnu_p_2p, cmbps_mnu_p_2p)
            

cdef double[:] cosm_par_mnu_m_2m
cdef double C_mnu_m_2m
cdef double[:, :] a_data_mnu_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_2m, C_mnu_m_2m, a_data_mnu_m_2m, b_data_mnu_m_2m, c_data_mnu_m_2m, lps_cc_data_mnu_m_2m, lps_cs_data_mnu_m_2m, lps_ss_data_mnu_m_2m, scale_factor_data_mnu_m_2m, window_c_data_mnu_m_2m, window_s_data_mnu_m_2m, mps_data_mnu_m_2m, z_at_chi_data_mnu_m_2m, cmbps_mnu_m_2m = data_import_func('data_mnu_m_2m')

cdef double lbs_mnu_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_2m, a_data_mnu_m_2m, b_data_mnu_m_2m, c_data_mnu_m_2m, scale_factor_data_mnu_m_2m, window_c_data_mnu_m_2m, window_s_data_mnu_m_2m, mps_data_mnu_m_2m, z_at_chi_data_mnu_m_2m)

cdef double lps_mnu_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_2m, lps_cs_data_mnu_m_2m, lps_ss_data_mnu_m_2m, cmbps_mnu_m_2m)
            

cdef double[:] cosm_par_mnu_m_1m
cdef double C_mnu_m_1m
cdef double[:, :] a_data_mnu_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_1m, C_mnu_m_1m, a_data_mnu_m_1m, b_data_mnu_m_1m, c_data_mnu_m_1m, lps_cc_data_mnu_m_1m, lps_cs_data_mnu_m_1m, lps_ss_data_mnu_m_1m, scale_factor_data_mnu_m_1m, window_c_data_mnu_m_1m, window_s_data_mnu_m_1m, mps_data_mnu_m_1m, z_at_chi_data_mnu_m_1m, cmbps_mnu_m_1m = data_import_func('data_mnu_m_1m')

cdef double lbs_mnu_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_1m, a_data_mnu_m_1m, b_data_mnu_m_1m, c_data_mnu_m_1m, scale_factor_data_mnu_m_1m, window_c_data_mnu_m_1m, window_s_data_mnu_m_1m, mps_data_mnu_m_1m, z_at_chi_data_mnu_m_1m)

cdef double lps_mnu_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_1m, lps_cs_data_mnu_m_1m, lps_ss_data_mnu_m_1m, cmbps_mnu_m_1m)
            

cdef double[:] cosm_par_mnu_m_0
cdef double C_mnu_m_0
cdef double[:, :] a_data_mnu_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_0, C_mnu_m_0, a_data_mnu_m_0, b_data_mnu_m_0, c_data_mnu_m_0, lps_cc_data_mnu_m_0, lps_cs_data_mnu_m_0, lps_ss_data_mnu_m_0, scale_factor_data_mnu_m_0, window_c_data_mnu_m_0, window_s_data_mnu_m_0, mps_data_mnu_m_0, z_at_chi_data_mnu_m_0, cmbps_mnu_m_0 = data_import_func('data_mnu_m_0')

cdef double lbs_mnu_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_0, a_data_mnu_m_0, b_data_mnu_m_0, c_data_mnu_m_0, scale_factor_data_mnu_m_0, window_c_data_mnu_m_0, window_s_data_mnu_m_0, mps_data_mnu_m_0, z_at_chi_data_mnu_m_0)

cdef double lps_mnu_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_0, lps_cs_data_mnu_m_0, lps_ss_data_mnu_m_0, cmbps_mnu_m_0)
            

cdef double[:] cosm_par_mnu_m_1p
cdef double C_mnu_m_1p
cdef double[:, :] a_data_mnu_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_1p, C_mnu_m_1p, a_data_mnu_m_1p, b_data_mnu_m_1p, c_data_mnu_m_1p, lps_cc_data_mnu_m_1p, lps_cs_data_mnu_m_1p, lps_ss_data_mnu_m_1p, scale_factor_data_mnu_m_1p, window_c_data_mnu_m_1p, window_s_data_mnu_m_1p, mps_data_mnu_m_1p, z_at_chi_data_mnu_m_1p, cmbps_mnu_m_1p = data_import_func('data_mnu_m_1p')

cdef double lbs_mnu_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_1p, a_data_mnu_m_1p, b_data_mnu_m_1p, c_data_mnu_m_1p, scale_factor_data_mnu_m_1p, window_c_data_mnu_m_1p, window_s_data_mnu_m_1p, mps_data_mnu_m_1p, z_at_chi_data_mnu_m_1p)

cdef double lps_mnu_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_1p, lps_cs_data_mnu_m_1p, lps_ss_data_mnu_m_1p, cmbps_mnu_m_1p)
            

cdef double[:] cosm_par_mnu_m_2p
cdef double C_mnu_m_2p
cdef double[:, :] a_data_mnu_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_2p, C_mnu_m_2p, a_data_mnu_m_2p, b_data_mnu_m_2p, c_data_mnu_m_2p, lps_cc_data_mnu_m_2p, lps_cs_data_mnu_m_2p, lps_ss_data_mnu_m_2p, scale_factor_data_mnu_m_2p, window_c_data_mnu_m_2p, window_s_data_mnu_m_2p, mps_data_mnu_m_2p, z_at_chi_data_mnu_m_2p, cmbps_mnu_m_2p = data_import_func('data_mnu_m_2p')

cdef double lbs_mnu_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_2p, a_data_mnu_m_2p, b_data_mnu_m_2p, c_data_mnu_m_2p, scale_factor_data_mnu_m_2p, window_c_data_mnu_m_2p, window_s_data_mnu_m_2p, mps_data_mnu_m_2p, z_at_chi_data_mnu_m_2p)

cdef double lps_mnu_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_2p, lps_cs_data_mnu_m_2p, lps_ss_data_mnu_m_2p, cmbps_mnu_m_2p)
            

cdef double[:] cosm_par_w0_p_2m
cdef double C_w0_p_2m
cdef double[:, :] a_data_w0_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_2m, C_w0_p_2m, a_data_w0_p_2m, b_data_w0_p_2m, c_data_w0_p_2m, lps_cc_data_w0_p_2m, lps_cs_data_w0_p_2m, lps_ss_data_w0_p_2m, scale_factor_data_w0_p_2m, window_c_data_w0_p_2m, window_s_data_w0_p_2m, mps_data_w0_p_2m, z_at_chi_data_w0_p_2m, cmbps_w0_p_2m = data_import_func('data_w0_p_2m')

cdef double lbs_w0_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_2m, a_data_w0_p_2m, b_data_w0_p_2m, c_data_w0_p_2m, scale_factor_data_w0_p_2m, window_c_data_w0_p_2m, window_s_data_w0_p_2m, mps_data_w0_p_2m, z_at_chi_data_w0_p_2m)

cdef double lps_w0_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_2m, lps_cs_data_w0_p_2m, lps_ss_data_w0_p_2m, cmbps_w0_p_2m)
            

cdef double[:] cosm_par_w0_p_1m
cdef double C_w0_p_1m
cdef double[:, :] a_data_w0_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_1m, C_w0_p_1m, a_data_w0_p_1m, b_data_w0_p_1m, c_data_w0_p_1m, lps_cc_data_w0_p_1m, lps_cs_data_w0_p_1m, lps_ss_data_w0_p_1m, scale_factor_data_w0_p_1m, window_c_data_w0_p_1m, window_s_data_w0_p_1m, mps_data_w0_p_1m, z_at_chi_data_w0_p_1m, cmbps_w0_p_1m = data_import_func('data_w0_p_1m')

cdef double lbs_w0_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_1m, a_data_w0_p_1m, b_data_w0_p_1m, c_data_w0_p_1m, scale_factor_data_w0_p_1m, window_c_data_w0_p_1m, window_s_data_w0_p_1m, mps_data_w0_p_1m, z_at_chi_data_w0_p_1m)

cdef double lps_w0_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_1m, lps_cs_data_w0_p_1m, lps_ss_data_w0_p_1m, cmbps_w0_p_1m)
            

cdef double[:] cosm_par_w0_p_0
cdef double C_w0_p_0
cdef double[:, :] a_data_w0_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_0, C_w0_p_0, a_data_w0_p_0, b_data_w0_p_0, c_data_w0_p_0, lps_cc_data_w0_p_0, lps_cs_data_w0_p_0, lps_ss_data_w0_p_0, scale_factor_data_w0_p_0, window_c_data_w0_p_0, window_s_data_w0_p_0, mps_data_w0_p_0, z_at_chi_data_w0_p_0, cmbps_w0_p_0 = data_import_func('data_w0_p_0')

cdef double lbs_w0_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_0, a_data_w0_p_0, b_data_w0_p_0, c_data_w0_p_0, scale_factor_data_w0_p_0, window_c_data_w0_p_0, window_s_data_w0_p_0, mps_data_w0_p_0, z_at_chi_data_w0_p_0)

cdef double lps_w0_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_0, lps_cs_data_w0_p_0, lps_ss_data_w0_p_0, cmbps_w0_p_0)
            

cdef double[:] cosm_par_w0_p_1p
cdef double C_w0_p_1p
cdef double[:, :] a_data_w0_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_1p, C_w0_p_1p, a_data_w0_p_1p, b_data_w0_p_1p, c_data_w0_p_1p, lps_cc_data_w0_p_1p, lps_cs_data_w0_p_1p, lps_ss_data_w0_p_1p, scale_factor_data_w0_p_1p, window_c_data_w0_p_1p, window_s_data_w0_p_1p, mps_data_w0_p_1p, z_at_chi_data_w0_p_1p, cmbps_w0_p_1p = data_import_func('data_w0_p_1p')

cdef double lbs_w0_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_1p, a_data_w0_p_1p, b_data_w0_p_1p, c_data_w0_p_1p, scale_factor_data_w0_p_1p, window_c_data_w0_p_1p, window_s_data_w0_p_1p, mps_data_w0_p_1p, z_at_chi_data_w0_p_1p)

cdef double lps_w0_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_1p, lps_cs_data_w0_p_1p, lps_ss_data_w0_p_1p, cmbps_w0_p_1p)
            

cdef double[:] cosm_par_w0_p_2p
cdef double C_w0_p_2p
cdef double[:, :] a_data_w0_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_2p, C_w0_p_2p, a_data_w0_p_2p, b_data_w0_p_2p, c_data_w0_p_2p, lps_cc_data_w0_p_2p, lps_cs_data_w0_p_2p, lps_ss_data_w0_p_2p, scale_factor_data_w0_p_2p, window_c_data_w0_p_2p, window_s_data_w0_p_2p, mps_data_w0_p_2p, z_at_chi_data_w0_p_2p, cmbps_w0_p_2p = data_import_func('data_w0_p_2p')

cdef double lbs_w0_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_2p, a_data_w0_p_2p, b_data_w0_p_2p, c_data_w0_p_2p, scale_factor_data_w0_p_2p, window_c_data_w0_p_2p, window_s_data_w0_p_2p, mps_data_w0_p_2p, z_at_chi_data_w0_p_2p)

cdef double lps_w0_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_2p, lps_cs_data_w0_p_2p, lps_ss_data_w0_p_2p, cmbps_w0_p_2p)
            

cdef double[:] cosm_par_w0_m_2m
cdef double C_w0_m_2m
cdef double[:, :] a_data_w0_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_2m, C_w0_m_2m, a_data_w0_m_2m, b_data_w0_m_2m, c_data_w0_m_2m, lps_cc_data_w0_m_2m, lps_cs_data_w0_m_2m, lps_ss_data_w0_m_2m, scale_factor_data_w0_m_2m, window_c_data_w0_m_2m, window_s_data_w0_m_2m, mps_data_w0_m_2m, z_at_chi_data_w0_m_2m, cmbps_w0_m_2m = data_import_func('data_w0_m_2m')

cdef double lbs_w0_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_2m, a_data_w0_m_2m, b_data_w0_m_2m, c_data_w0_m_2m, scale_factor_data_w0_m_2m, window_c_data_w0_m_2m, window_s_data_w0_m_2m, mps_data_w0_m_2m, z_at_chi_data_w0_m_2m)

cdef double lps_w0_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_2m, lps_cs_data_w0_m_2m, lps_ss_data_w0_m_2m, cmbps_w0_m_2m)
            

cdef double[:] cosm_par_w0_m_1m
cdef double C_w0_m_1m
cdef double[:, :] a_data_w0_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_1m, C_w0_m_1m, a_data_w0_m_1m, b_data_w0_m_1m, c_data_w0_m_1m, lps_cc_data_w0_m_1m, lps_cs_data_w0_m_1m, lps_ss_data_w0_m_1m, scale_factor_data_w0_m_1m, window_c_data_w0_m_1m, window_s_data_w0_m_1m, mps_data_w0_m_1m, z_at_chi_data_w0_m_1m, cmbps_w0_m_1m = data_import_func('data_w0_m_1m')

cdef double lbs_w0_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_1m, a_data_w0_m_1m, b_data_w0_m_1m, c_data_w0_m_1m, scale_factor_data_w0_m_1m, window_c_data_w0_m_1m, window_s_data_w0_m_1m, mps_data_w0_m_1m, z_at_chi_data_w0_m_1m)

cdef double lps_w0_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_1m, lps_cs_data_w0_m_1m, lps_ss_data_w0_m_1m, cmbps_w0_m_1m)
            

cdef double[:] cosm_par_w0_m_0
cdef double C_w0_m_0
cdef double[:, :] a_data_w0_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_0, C_w0_m_0, a_data_w0_m_0, b_data_w0_m_0, c_data_w0_m_0, lps_cc_data_w0_m_0, lps_cs_data_w0_m_0, lps_ss_data_w0_m_0, scale_factor_data_w0_m_0, window_c_data_w0_m_0, window_s_data_w0_m_0, mps_data_w0_m_0, z_at_chi_data_w0_m_0, cmbps_w0_m_0 = data_import_func('data_w0_m_0')

cdef double lbs_w0_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_0, a_data_w0_m_0, b_data_w0_m_0, c_data_w0_m_0, scale_factor_data_w0_m_0, window_c_data_w0_m_0, window_s_data_w0_m_0, mps_data_w0_m_0, z_at_chi_data_w0_m_0)

cdef double lps_w0_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_0, lps_cs_data_w0_m_0, lps_ss_data_w0_m_0, cmbps_w0_m_0)
            

cdef double[:] cosm_par_w0_m_1p
cdef double C_w0_m_1p
cdef double[:, :] a_data_w0_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_1p, C_w0_m_1p, a_data_w0_m_1p, b_data_w0_m_1p, c_data_w0_m_1p, lps_cc_data_w0_m_1p, lps_cs_data_w0_m_1p, lps_ss_data_w0_m_1p, scale_factor_data_w0_m_1p, window_c_data_w0_m_1p, window_s_data_w0_m_1p, mps_data_w0_m_1p, z_at_chi_data_w0_m_1p, cmbps_w0_m_1p = data_import_func('data_w0_m_1p')

cdef double lbs_w0_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_1p, a_data_w0_m_1p, b_data_w0_m_1p, c_data_w0_m_1p, scale_factor_data_w0_m_1p, window_c_data_w0_m_1p, window_s_data_w0_m_1p, mps_data_w0_m_1p, z_at_chi_data_w0_m_1p)

cdef double lps_w0_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_1p, lps_cs_data_w0_m_1p, lps_ss_data_w0_m_1p, cmbps_w0_m_1p)
            

cdef double[:] cosm_par_w0_m_2p
cdef double C_w0_m_2p
cdef double[:, :] a_data_w0_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_2p, C_w0_m_2p, a_data_w0_m_2p, b_data_w0_m_2p, c_data_w0_m_2p, lps_cc_data_w0_m_2p, lps_cs_data_w0_m_2p, lps_ss_data_w0_m_2p, scale_factor_data_w0_m_2p, window_c_data_w0_m_2p, window_s_data_w0_m_2p, mps_data_w0_m_2p, z_at_chi_data_w0_m_2p, cmbps_w0_m_2p = data_import_func('data_w0_m_2p')

cdef double lbs_w0_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_2p, a_data_w0_m_2p, b_data_w0_m_2p, c_data_w0_m_2p, scale_factor_data_w0_m_2p, window_c_data_w0_m_2p, window_s_data_w0_m_2p, mps_data_w0_m_2p, z_at_chi_data_w0_m_2p)

cdef double lps_w0_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_2p, lps_cs_data_w0_m_2p, lps_ss_data_w0_m_2p, cmbps_w0_m_2p)
            
            
##########################################


cdef double der_2o(double f2p, double f1p, double f1m, double f2m, double dx) noexcept nogil:
    return (-1 * f2p + 8 * f1p - 8 * f1m + f2m) / (12 * dx)

cpdef double der(double fp, double fm, double dx) noexcept nogil:
    return (fp - fm) / (2 * dx)

# cdef double[:] fiducial_cosm_par = np.array([67.4, 0.0224, 0.120, 0.965, 2.1e-9, 0.06])


cdef double lps_der(int k, char* type1, char* type2, char* par, int delta_delta) noexcept nogil:

    if delta_delta == -2:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lps_f(k, type1, type2)

        if par[0] == b'H':
            return der(lps_H_p_2m(k, type1, type2), lps_H_m_2m(k, type1, type2), cosm_par_H_p_2m[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lps_ombh2_p_2m(k, type1, type2), lps_ombh2_m_2m(k, type1, type2), cosm_par_ombh2_p_2m[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lps_omch2_p_2m(k, type1, type2), lps_omch2_m_2m(k, type1, type2), cosm_par_omch2_p_2m[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lps_ns_p_2m(k, type1, type2), lps_ns_m_2m(k, type1, type2), cosm_par_ns_p_2m[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lps_As_p_2m(k, type1, type2), lps_As_m_2m(k, type1, type2), cosm_par_As_p_2m[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lps_mnu_p_2m(k, type1, type2), lps_mnu_m_2m(k, type1, type2), cosm_par_mnu_p_2m[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lps_w0_p_2m(k, type1, type2), lps_w0_m_2m(k, type1, type2), cosm_par_w0_p_2m[6] - cosm_par_f[6])
    if delta_delta == -1:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lps_f(k, type1, type2)

        if par[0] == b'H':
            return der(lps_H_p_1m(k, type1, type2), lps_H_m_1m(k, type1, type2), cosm_par_H_p_1m[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lps_ombh2_p_1m(k, type1, type2), lps_ombh2_m_1m(k, type1, type2), cosm_par_ombh2_p_1m[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lps_omch2_p_1m(k, type1, type2), lps_omch2_m_1m(k, type1, type2), cosm_par_omch2_p_1m[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lps_ns_p_1m(k, type1, type2), lps_ns_m_1m(k, type1, type2), cosm_par_ns_p_1m[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lps_As_p_1m(k, type1, type2), lps_As_m_1m(k, type1, type2), cosm_par_As_p_1m[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lps_mnu_p_1m(k, type1, type2), lps_mnu_m_1m(k, type1, type2), cosm_par_mnu_p_1m[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lps_w0_p_1m(k, type1, type2), lps_w0_m_1m(k, type1, type2), cosm_par_w0_p_1m[6] - cosm_par_f[6])
    if delta_delta == 0:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lps_f(k, type1, type2)

        if par[0] == b'H':
            return der(lps_H_p_0(k, type1, type2), lps_H_m_0(k, type1, type2), cosm_par_H_p_0[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lps_ombh2_p_0(k, type1, type2), lps_ombh2_m_0(k, type1, type2), cosm_par_ombh2_p_0[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lps_omch2_p_0(k, type1, type2), lps_omch2_m_0(k, type1, type2), cosm_par_omch2_p_0[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lps_ns_p_0(k, type1, type2), lps_ns_m_0(k, type1, type2), cosm_par_ns_p_0[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lps_As_p_0(k, type1, type2), lps_As_m_0(k, type1, type2), cosm_par_As_p_0[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lps_mnu_p_0(k, type1, type2), lps_mnu_m_0(k, type1, type2), cosm_par_mnu_p_0[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lps_w0_p_0(k, type1, type2), lps_w0_m_0(k, type1, type2), cosm_par_w0_p_0[6] - cosm_par_f[6])
    if delta_delta == 1:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lps_f(k, type1, type2)

        if par[0] == b'H':
            return der(lps_H_p_1p(k, type1, type2), lps_H_m_1p(k, type1, type2), cosm_par_H_p_1p[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lps_ombh2_p_1p(k, type1, type2), lps_ombh2_m_1p(k, type1, type2), cosm_par_ombh2_p_1p[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lps_omch2_p_1p(k, type1, type2), lps_omch2_m_1p(k, type1, type2), cosm_par_omch2_p_1p[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lps_ns_p_1p(k, type1, type2), lps_ns_m_1p(k, type1, type2), cosm_par_ns_p_1p[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lps_As_p_1p(k, type1, type2), lps_As_m_1p(k, type1, type2), cosm_par_As_p_1p[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lps_mnu_p_1p(k, type1, type2), lps_mnu_m_1p(k, type1, type2), cosm_par_mnu_p_1p[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lps_w0_p_1p(k, type1, type2), lps_w0_m_1p(k, type1, type2), cosm_par_w0_p_1p[6] - cosm_par_f[6])
    if delta_delta == 2:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lps_f(k, type1, type2)

        if par[0] == b'H':
            return der(lps_H_p_2p(k, type1, type2), lps_H_m_2p(k, type1, type2), cosm_par_H_p_2p[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lps_ombh2_p_2p(k, type1, type2), lps_ombh2_m_2p(k, type1, type2), cosm_par_ombh2_p_2p[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lps_omch2_p_2p(k, type1, type2), lps_omch2_m_2p(k, type1, type2), cosm_par_omch2_p_2p[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lps_ns_p_2p(k, type1, type2), lps_ns_m_2p(k, type1, type2), cosm_par_ns_p_2p[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lps_As_p_2p(k, type1, type2), lps_As_m_2p(k, type1, type2), cosm_par_As_p_2p[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lps_mnu_p_2p(k, type1, type2), lps_mnu_m_2p(k, type1, type2), cosm_par_mnu_p_2p[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lps_w0_p_2p(k, type1, type2), lps_w0_m_2p(k, type1, type2), cosm_par_w0_p_2p[6] - cosm_par_f[6])


cdef double lbs_der(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, char* par, int delta_delta) noexcept nogil:
    
    if delta_delta == -2:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lbs_f(k1, k2, k3, type1, type2, type3, num_samples)

        if par[0] == b'H':
            return der(lbs_H_p_2m(k1, k2, k3, type1, type2, type3, num_samples), lbs_H_m_2m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_H_p_2m[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lbs_ombh2_p_2m(k1, k2, k3, type1, type2, type3, num_samples), lbs_ombh2_m_2m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ombh2_p_2m[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lbs_omch2_p_2m(k1, k2, k3, type1, type2, type3, num_samples), lbs_omch2_m_2m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_omch2_p_2m[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lbs_ns_p_2m(k1, k2, k3, type1, type2, type3, num_samples), lbs_ns_m_2m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ns_p_2m[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lbs_As_p_2m(k1, k2, k3, type1, type2, type3, num_samples), lbs_As_m_2m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_As_p_2m[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lbs_mnu_p_2m(k1, k2, k3, type1, type2, type3, num_samples), lbs_mnu_m_2m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_mnu_p_2m[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lbs_w0_p_2m(k1, k2, k3, type1, type2, type3, num_samples), lbs_w0_m_2m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_w0_p_2m[6] - cosm_par_f[6])
    if delta_delta == -1:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lbs_f(k1, k2, k3, type1, type2, type3, num_samples)

        if par[0] == b'H':
            return der(lbs_H_p_1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_H_m_1m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_H_p_1m[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lbs_ombh2_p_1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_ombh2_m_1m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ombh2_p_1m[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lbs_omch2_p_1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_omch2_m_1m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_omch2_p_1m[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lbs_ns_p_1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_ns_m_1m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ns_p_1m[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lbs_As_p_1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_As_m_1m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_As_p_1m[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lbs_mnu_p_1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_mnu_m_1m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_mnu_p_1m[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lbs_w0_p_1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_w0_m_1m(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_w0_p_1m[6] - cosm_par_f[6])
    if delta_delta == 0:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lbs_f(k1, k2, k3, type1, type2, type3, num_samples)

        if par[0] == b'H':
            return der(lbs_H_p_0(k1, k2, k3, type1, type2, type3, num_samples), lbs_H_m_0(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_H_p_0[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lbs_ombh2_p_0(k1, k2, k3, type1, type2, type3, num_samples), lbs_ombh2_m_0(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ombh2_p_0[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lbs_omch2_p_0(k1, k2, k3, type1, type2, type3, num_samples), lbs_omch2_m_0(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_omch2_p_0[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lbs_ns_p_0(k1, k2, k3, type1, type2, type3, num_samples), lbs_ns_m_0(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ns_p_0[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lbs_As_p_0(k1, k2, k3, type1, type2, type3, num_samples), lbs_As_m_0(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_As_p_0[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lbs_mnu_p_0(k1, k2, k3, type1, type2, type3, num_samples), lbs_mnu_m_0(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_mnu_p_0[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lbs_w0_p_0(k1, k2, k3, type1, type2, type3, num_samples), lbs_w0_m_0(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_w0_p_0[6] - cosm_par_f[6])
    if delta_delta == 1:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lbs_f(k1, k2, k3, type1, type2, type3, num_samples)

        if par[0] == b'H':
            return der(lbs_H_p_1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_H_m_1p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_H_p_1p[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lbs_ombh2_p_1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_ombh2_m_1p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ombh2_p_1p[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lbs_omch2_p_1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_omch2_m_1p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_omch2_p_1p[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lbs_ns_p_1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_ns_m_1p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ns_p_1p[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lbs_As_p_1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_As_m_1p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_As_p_1p[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lbs_mnu_p_1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_mnu_m_1p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_mnu_p_1p[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lbs_w0_p_1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_w0_m_1p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_w0_p_1p[6] - cosm_par_f[6])
    if delta_delta == 2:
        # for b'snr' case, just returns the original function
        if par[0] == b's':
            return lbs_f(k1, k2, k3, type1, type2, type3, num_samples)

        if par[0] == b'H':
            return der(lbs_H_p_2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_H_m_2p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_H_p_2p[0] - cosm_par_f[0])
                
        if par[0] == b'o' and par[2] == b'b':
            return der(lbs_ombh2_p_2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_ombh2_m_2p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ombh2_p_2p[1] - cosm_par_f[1])

        if par[0] == b'o' and par[2] == b'c':
            return der(lbs_omch2_p_2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_omch2_m_2p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_omch2_p_2p[2] - cosm_par_f[2])

        if par[0] == b'n':
            return der(lbs_ns_p_2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_ns_m_2p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_ns_p_2p[3] - cosm_par_f[3])

        if par[0] == b'A':
            return der(lbs_As_p_2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_As_m_2p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_As_p_2p[4] - cosm_par_f[4])

        if par[0] == b'm':
            return der(lbs_mnu_p_2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_mnu_m_2p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_mnu_p_2p[5] - cosm_par_f[5])

        if par[0] == b'w':
            return der(lbs_w0_p_2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_w0_m_2p(k1, k2, k3, type1, type2, type3, num_samples), cosm_par_w0_p_2p[6] - cosm_par_f[6])

#######################

cpdef double a_test(double l, double z):
    return a(l, z, a_data_f)

cpdef double b_test(double l, double z):
    return b(l, z, b_data_f)

cpdef double c_test(double l, double z):
    return c(l, z, c_data_f)

def mps_test(l, z):
    return matter_power_spectrum(l, z, mps_data_f)

def mbs_test(k1, k2, k3, z):
    return mbs(k1, k2, k3, z, mps_data_f, a_data_f, b_data_f, c_data_f)

def window_func_test(chi, type):
    return window_func(chi, type, window_c_data=window_c_data_f, window_s_data=window_s_data_f)

def scale_factor_test(chi):
    return scale_factor(chi, scale_factor_data=scale_factor_data_f)

def z_at_chi_test(chi):
    return z_at_chi(chi, z_at_chi_data=z_at_chi_data_f)

def lbs_integrand_test(chi, l1, l2, l3, type1, type2, type3):
    return lbs_integrand(chi, l1, l2, l3, type1, type2, type3, 
    a_data_f, b_data_f, c_data_f, scale_factor_data_f, window_c_data_f, window_s_data_f, mps_data_f, z_at_chi_data_f)

def lps_f_obs_test(l, type1, type2):
    return lps_f_obs(l, type1, type2)

def lps_der_test(k, type1, type2, par, deltadelta):
    return lps_der(k, type1, type2, par, deltadelta)

def lbs_der_test(k1, k2, k3, type1, type2, type3, num_samples, par, deltadelta):
    return lbs_der(k1, k2, k3, type1, type2, type3, num_samples, par, deltadelta)


# some code to get shapes of all arrays for debugging purposes

# # Assuming the arrays are already defined
# arrays = {
#     "a_data": a_data_f,
#     "b_data": b_data_f,
#     "c_data": c_data_f,
#     "lps_cc_data": lps_cc_data_f,
#     "lps_cs_data": lps_cs_data_f,
#     "lps_ss_data": lps_ss_data_f,
#     "scale_factor_data": scale_factor_data_f,
#     "window_c_data": window_c_data_f,
#     "window_s_data": window_s_data_f,
#     "mps_data": mps_data_f,
#     "z_at_chi_data": z_at_chi_data_f,
# }

# # Print the shape of each array
# for name, array in arrays.items():
#     print(f"{name}: {array.shape}")
#     print(f'contig: {array.is_c_contig()}')
#     arr = np.array(array)
#     nan_indices = np.where(np.isnan(arr))
#     nan_coords = list(zip(*nan_indices))
#     print(f'max val: {nan_coords}')
