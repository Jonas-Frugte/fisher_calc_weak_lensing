print('wahoo!')
import json
cimport interpolation as interp
import math
import numpy as np
from libc.stdio cimport printf
from libc.math cimport fmax, fabs

from libc.stdlib cimport getenv

# declare interpolation parameters, NEEDS TO BE SAME AS INTERPOLATION 
cdef double k_min = 1e-4 # is now 1e-4
cdef double k_max = 2000.0 # 1000 for rough option, 2000 for fine option
cdef int k_num = 300 # * 2 # for finer data option
cdef int k_num_fine = k_num * 25
cdef double chi_min = 1e-8 # is now 1e-8
cdef double chi_max = 14000.0
cdef int chi_num = 100 * 50 # * 2 # for finer data option
cdef double z_min = 1e-12
cdef double z_max = 1100.0
cdef int z_num = 100 # * 2 # for finer data option
cdef int z_num_fine = z_num * 25

cdef float window_func_c_scale_factor = 1.

def nan_to_val(x, val = 0.0):
    if math.isnan(x):
        return val
    else:
        return x

def get_k_max():
    return k_max

# specifies the noise that will be used for cmb lensing
cdef int conv_noise_type = 1
# type values correspond to:
# 0: S0 noise curves (old, remove)
# 1: sigma = 1, Delta_P = 6 (stage 3 toshiya)
# 2: sigma = 3, Delta_P = 1 (stage 4 toshiya)

# (units are in microKelvin arcmin and arcmin)

# SO noise
cdef str folder_file_path = '/scratch/p319950/data_rough/'
cdef str filepath_convergence_noise_file_path = folder_file_path + 'conv_noise.dat'
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

# perhaps we can export data as c arrays instead of as memory views in the future so that we can specify return types like below and can avoid having to 
# declare types of all data before every data import
#cdef (dict[str, double], double, double[:, :], double[:, :], double[:], double[:], double[:], double[:], double[:], double[:], double[:, :], double[:]) data_importer(str folder_name):
def data_import_func(folder_name):
    # declare data types

    cdef dict[str, double] cosm_par
    cdef double C
    cdef double[:, :] a_data = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
    cdef double[:, :] b_data = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
    cdef double[:, :] c_data = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
    cdef double[:] lps_cc_data = np.zeros(k_num, dtype=np.float64)
    cdef double[:] lps_cs_data = np.zeros(k_num, dtype=np.float64)
    cdef double[:] lps_ss_data = np.zeros(k_num, dtype=np.float64)
    cdef double[:] scale_factor_data = np.zeros(chi_num, dtype=np.float64)
    cdef double[:] window_c_data = np.zeros(chi_num, dtype=np.float64)
    cdef double[:] window_s_data = np.zeros(chi_num, dtype=np.float64)
    cdef double[:, :] mps_data = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
    cdef double[:] z_at_chi_data = np.zeros(chi_num, dtype=np.float64)

    # import all relevant data

    # indexes used for populating C arrays

    cdef int i
    cdef int j

    # tmp_dir = getenv("TMPDIR")
    # cdef str filepath = tmp_dir.decode("utf-8") + '/ResearchProject/data/' + folder_name
    
    # cdef str filepath = '/home3/p319950/ResearchProject/data/' + folder_name
    cdef str filepath = folder_file_path + folder_name
    print(filepath)

    # order of inputs: H0, ombh2, omch2, ns, redshifts, As, mnu,
    with open(filepath + '/cosm_par.json', 'r') as file:
        cosm_par = json.load(file)

    # all 2d funcs are made such that order of arguments is (k, z)
    C = 0.
    with open(filepath + '/rho_bar', 'r') as file:
        C = float(file.read())

    # IMPORTANT: a, b, c require first index to refer to ks because that is logarithmically spaced and second index to refer to zs because that is evenly spaced
    # export of a, b, c looks like data = [[cosm_data.a(k, z) for k in ks] for z in zs] so we need to take the transpose, i.e. a_data[j][i] instead of [i][j]

    # For 2D arrays like a_data, b_data, c_data
    with open(filepath + '/a', 'r') as file:
        i = 0
        for line in file:
            #print(i)
            j = 0
            for num in line.split():
                a_data[j][i] = float(nan_to_val(float(num))) # second conversion to float is a good luck charm
                j += 1
            i += 1

    with open(filepath + '/b', 'r') as file:
        i = 0
        for line in file:
            j = 0
            for num in line.split():
                b_data[j][i] = np.float64(num)
                j += 1
            i += 1
    
    with open(filepath + '/c', 'r') as file:
        i = 0
        for line in file:
            j = 0
            for num in line.split():
                c_data[j][i] = float(num)
                j += 1
            i += 1

    # For 1D arrays like lps_cc_data, lps_cs_data, etc.
    with open(filepath + '/lensing_power_spectrum_cc', 'r') as file:
        i = 0
        for line in file:
            lps_cc_data[i] = float(line)
            i += 1

    with open(filepath + '/lensing_power_spectrum_cs', 'r') as file:
        i = 0
        for line in file:
            lps_cs_data[i] = float(line)
            i += 1

    with open(filepath + '/lensing_power_spectrum_ss', 'r') as file:
        i = 0
        for line in file:
            lps_ss_data[i] = float(line)
            i += 1

    # For window_c_data, window_s_data, etc.
    with open(filepath + '/window_func_c', 'r') as file:
        i = 0
        for line in file:
            window_c_data[i] = float(line)
            i += 1

    with open(filepath + '/window_func_s', 'r') as file:
        i = 0
        for line in file:
            window_s_data[i] = float(nan_to_val(float(line)))
            i += 1

    # For z_at_chi_data
    with open(filepath + '/z_at_chi', 'r') as file:
        i = 0
        for line in file:
            z_at_chi_data[i] = float(line)
            i += 1

    # For the matter power spectrum, using the same approach as discussed earlier
    with open(filepath + '/matter_power_spectrum', 'r') as file:
        i = 0
        for line in file:
            j = 0
            for num in line.split():
                mps_data[i][j] = float(num)
                j += 1
            i += 1

    with open(filepath + '/scale_factor', 'r') as file:
        i = 0
        for line in file:
            scale_factor_data[i] = float(line)
            i += 1

    return cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data



# turns data into functions 

cpdef double a(double k, double z, double[:, :] a_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, a_data)
cpdef double b(double k, double z, double[:, :] b_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, b_data)
cpdef double c(double k, double z, double[:, :] c_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, c_data)

# without gil cython doesn't like strings bc they are python objects, you can use char* C objects (C like strings)

cpdef double lps(double l, char* type1, char* type2, double[:] lps_cc_data, double[:] lps_cs_data, double[:] lps_ss_data) noexcept nogil:
    # CAUTION: this function doesn't throw error messages if type1 or type2 are not of type convergence or shear
    # for optimization purposes
    if type1 == type2:
        if type1[0] == b'c':
            return interp.logspace_linear_interp(l, k_min, k_max, k_num, lps_cc_data)
        elif type1[0] == b's':
            return interp.logspace_linear_interp(l, k_min, k_max, k_num, lps_ss_data)
    else:
        return interp.logspace_linear_interp(l, k_min, k_max, k_num, lps_cs_data)

cdef double scale_factor(double chi, double[:] scale_factor_data) noexcept nogil:    
    return interp.linear_interp(chi, chi_min, chi_max, chi_num, scale_factor_data)

cdef double window_func(double chi, char* type, double[:] window_c_data, double[:] window_s_data) noexcept nogil:
    if type[0] == b'c':
        return window_func_c_scale_factor * fmax(0., interp.logspace_linear_interp(chi, chi_min, chi_max, chi_num, window_c_data))
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

cdef dict[str, double] cosm_par_f
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
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_f, C_f, a_data_f, b_data_f, c_data_f, lps_cc_data_f, lps_cs_data_f, lps_ss_data_f, scale_factor_data_f, window_c_data_f, window_s_data_f, mps_data_f, z_at_chi_data_f = data_import_func('data_fiducial')

cdef double lbs_f(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_f, a_data_f, b_data_f, c_data_f, scale_factor_data_f, window_c_data_f, window_s_data_f, mps_data_f, z_at_chi_data_f)

cdef double lps_f(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_f, lps_cs_data_f, lps_ss_data_f)

cdef double lps_noise(int l, char* type1, char* type2) noexcept nogil:
    cdef float noise = 0.
    if type1[0] == b'c' and type2[0] == b'c':
        # if we rescale psi_kappa, then we must also rescale the associated noise in the same way
        if conv_noise_type == 0:
            noise = window_func_c_scale_factor**2 * 4. * (l * 1.0)**(-2) * (l + 1.0)**(-2) * conv_noise_data[l-2, 7]
        elif conv_noise_type == 1:
            # stage 3 wide toshiya
            noise = interp.logspace_linear_interp(l, lmin_cmbn, lmax_cmbn, lnum_cmbn, cmbn_106)
        elif conv_noise_type == 2:
            # stage 4 toshiya
            noise = interp.logspace_linear_interp(l, lmin_cmbn, lmax_cmbn, lnum_cmbn, cmbn_301)

    if type1[0] == b's' and type2[0] == b's':
        if conv_noise_type == 1:
            noise = 4. * ((l - 1.) * l * (l + 1.) * (l + 2.))**(-1) * 0.3**2 / 5
        if conv_noise_type == 2:
            noise = 4. * ((l - 1.) * l * (l + 1.) * (l + 2.))**(-1) * 0.3**2 / 30 # 8 * 10.0 ** (-10) # this should actually be in sterradain, but that gives wrong results so for now we have it like this
    return noise

cdef double lps_f_obs(int l, char* type1, char* type2) noexcept nogil:

    #########################
    # NOISE: ON
    #########################

    return lps(l, type1, type2, lps_cc_data_f, lps_cs_data_f, lps_ss_data_f) + lps_noise(l, type1, type2)

##########################################

cdef dict[str, double] cosm_par_H2p
cdef double C_H2p
cdef double[:, :] a_data_H2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H2p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_H2p, C_H2p, a_data_H2p, b_data_H2p, c_data_H2p, lps_cc_data_H2p, lps_cs_data_H2p, lps_ss_data_H2p, scale_factor_data_H2p, window_c_data_H2p, window_s_data_H2p, mps_data_H2p, z_at_chi_data_H2p = data_import_func('data_H2p')

cpdef double lbs_H2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H2p, a_data_H2p, b_data_H2p, c_data_H2p, scale_factor_data_H2p, window_c_data_H2p, window_s_data_H2p, mps_data_H2p, z_at_chi_data_H2p)

cpdef double lps_H2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H2p, lps_cs_data_H2p, lps_ss_data_H2p)


cdef dict[str, double] cosm_par_H1p
cdef double C_H1p
cdef double[:, :] a_data_H1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H1p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_H1p, C_H1p, a_data_H1p, b_data_H1p, c_data_H1p, lps_cc_data_H1p, lps_cs_data_H1p, lps_ss_data_H1p, scale_factor_data_H1p, window_c_data_H1p, window_s_data_H1p, mps_data_H1p, z_at_chi_data_H1p = data_import_func('data_H1p')

cpdef double lbs_H1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H1p, a_data_H1p, b_data_H1p, c_data_H1p, scale_factor_data_H1p, window_c_data_H1p, window_s_data_H1p, mps_data_H1p, z_at_chi_data_H1p)

cpdef double lps_H1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H1p, lps_cs_data_H1p, lps_ss_data_H1p)


cdef dict[str, double] cosm_par_H1m
cdef double C_H1m
cdef double[:, :] a_data_H1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H1m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_H1m, C_H1m, a_data_H1m, b_data_H1m, c_data_H1m, lps_cc_data_H1m, lps_cs_data_H1m, lps_ss_data_H1m, scale_factor_data_H1m, window_c_data_H1m, window_s_data_H1m, mps_data_H1m, z_at_chi_data_H1m = data_import_func('data_H1m')

cpdef double lbs_H1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H1m, a_data_H1m, b_data_H1m, c_data_H1m, scale_factor_data_H1m, window_c_data_H1m, window_s_data_H1m, mps_data_H1m, z_at_chi_data_H1m)

cpdef double lps_H1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H1m, lps_cs_data_H1m, lps_ss_data_H1m)


cdef dict[str, double] cosm_par_H2m
cdef double C_H2m
cdef double[:, :] a_data_H2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H2m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_H2m, C_H2m, a_data_H2m, b_data_H2m, c_data_H2m, lps_cc_data_H2m, lps_cs_data_H2m, lps_ss_data_H2m, scale_factor_data_H2m, window_c_data_H2m, window_s_data_H2m, mps_data_H2m, z_at_chi_data_H2m = data_import_func('data_H2m')

cpdef double lbs_H2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H2m, a_data_H2m, b_data_H2m, c_data_H2m, scale_factor_data_H2m, window_c_data_H2m, window_s_data_H2m, mps_data_H2m, z_at_chi_data_H2m)

cpdef double lps_H2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H2m, lps_cs_data_H2m, lps_ss_data_H2m)


cdef dict[str, double] cosm_par_ombh22p
cdef double C_ombh22p
cdef double[:, :] a_data_ombh22p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh22p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh22p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh22p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh22p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh22p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh22p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh22p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh22p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh22p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh22p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_ombh22p, C_ombh22p, a_data_ombh22p, b_data_ombh22p, c_data_ombh22p, lps_cc_data_ombh22p, lps_cs_data_ombh22p, lps_ss_data_ombh22p, scale_factor_data_ombh22p, window_c_data_ombh22p, window_s_data_ombh22p, mps_data_ombh22p, z_at_chi_data_ombh22p = data_import_func('data_ombh22p')

cpdef double lbs_ombh22p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh22p, a_data_ombh22p, b_data_ombh22p, c_data_ombh22p, scale_factor_data_ombh22p, window_c_data_ombh22p, window_s_data_ombh22p, mps_data_ombh22p, z_at_chi_data_ombh22p)

cpdef double lps_ombh22p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh22p, lps_cs_data_ombh22p, lps_ss_data_ombh22p)


cdef dict[str, double] cosm_par_ombh21p
cdef double C_ombh21p
cdef double[:, :] a_data_ombh21p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh21p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh21p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh21p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh21p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh21p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh21p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh21p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh21p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh21p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh21p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_ombh21p, C_ombh21p, a_data_ombh21p, b_data_ombh21p, c_data_ombh21p, lps_cc_data_ombh21p, lps_cs_data_ombh21p, lps_ss_data_ombh21p, scale_factor_data_ombh21p, window_c_data_ombh21p, window_s_data_ombh21p, mps_data_ombh21p, z_at_chi_data_ombh21p = data_import_func('data_ombh21p')

cpdef double lbs_ombh21p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh21p, a_data_ombh21p, b_data_ombh21p, c_data_ombh21p, scale_factor_data_ombh21p, window_c_data_ombh21p, window_s_data_ombh21p, mps_data_ombh21p, z_at_chi_data_ombh21p)

cpdef double lps_ombh21p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh21p, lps_cs_data_ombh21p, lps_ss_data_ombh21p)


cdef dict[str, double] cosm_par_ombh21m
cdef double C_ombh21m
cdef double[:, :] a_data_ombh21m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh21m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh21m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh21m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh21m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh21m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh21m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh21m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh21m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh21m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh21m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_ombh21m, C_ombh21m, a_data_ombh21m, b_data_ombh21m, c_data_ombh21m, lps_cc_data_ombh21m, lps_cs_data_ombh21m, lps_ss_data_ombh21m, scale_factor_data_ombh21m, window_c_data_ombh21m, window_s_data_ombh21m, mps_data_ombh21m, z_at_chi_data_ombh21m = data_import_func('data_ombh21m')

cpdef double lbs_ombh21m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh21m, a_data_ombh21m, b_data_ombh21m, c_data_ombh21m, scale_factor_data_ombh21m, window_c_data_ombh21m, window_s_data_ombh21m, mps_data_ombh21m, z_at_chi_data_ombh21m)

cpdef double lps_ombh21m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh21m, lps_cs_data_ombh21m, lps_ss_data_ombh21m)


cdef dict[str, double] cosm_par_ombh22m
cdef double C_ombh22m
cdef double[:, :] a_data_ombh22m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh22m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh22m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh22m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh22m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh22m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh22m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh22m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh22m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh22m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh22m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_ombh22m, C_ombh22m, a_data_ombh22m, b_data_ombh22m, c_data_ombh22m, lps_cc_data_ombh22m, lps_cs_data_ombh22m, lps_ss_data_ombh22m, scale_factor_data_ombh22m, window_c_data_ombh22m, window_s_data_ombh22m, mps_data_ombh22m, z_at_chi_data_ombh22m = data_import_func('data_ombh22m')

cpdef double lbs_ombh22m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh22m, a_data_ombh22m, b_data_ombh22m, c_data_ombh22m, scale_factor_data_ombh22m, window_c_data_ombh22m, window_s_data_ombh22m, mps_data_ombh22m, z_at_chi_data_ombh22m)

cpdef double lps_ombh22m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh22m, lps_cs_data_ombh22m, lps_ss_data_ombh22m)


cdef dict[str, double] cosm_par_omch22p
cdef double C_omch22p
cdef double[:, :] a_data_omch22p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch22p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch22p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch22p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch22p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch22p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch22p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch22p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch22p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch22p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch22p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_omch22p, C_omch22p, a_data_omch22p, b_data_omch22p, c_data_omch22p, lps_cc_data_omch22p, lps_cs_data_omch22p, lps_ss_data_omch22p, scale_factor_data_omch22p, window_c_data_omch22p, window_s_data_omch22p, mps_data_omch22p, z_at_chi_data_omch22p = data_import_func('data_omch22p')

cpdef double lbs_omch22p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch22p, a_data_omch22p, b_data_omch22p, c_data_omch22p, scale_factor_data_omch22p, window_c_data_omch22p, window_s_data_omch22p, mps_data_omch22p, z_at_chi_data_omch22p)

cpdef double lps_omch22p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch22p, lps_cs_data_omch22p, lps_ss_data_omch22p)


cdef dict[str, double] cosm_par_omch21p
cdef double C_omch21p
cdef double[:, :] a_data_omch21p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch21p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch21p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch21p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch21p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch21p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch21p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch21p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch21p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch21p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch21p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_omch21p, C_omch21p, a_data_omch21p, b_data_omch21p, c_data_omch21p, lps_cc_data_omch21p, lps_cs_data_omch21p, lps_ss_data_omch21p, scale_factor_data_omch21p, window_c_data_omch21p, window_s_data_omch21p, mps_data_omch21p, z_at_chi_data_omch21p = data_import_func('data_omch21p')

cpdef double lbs_omch21p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch21p, a_data_omch21p, b_data_omch21p, c_data_omch21p, scale_factor_data_omch21p, window_c_data_omch21p, window_s_data_omch21p, mps_data_omch21p, z_at_chi_data_omch21p)

cpdef double lps_omch21p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch21p, lps_cs_data_omch21p, lps_ss_data_omch21p)


cdef dict[str, double] cosm_par_omch21m
cdef double C_omch21m
cdef double[:, :] a_data_omch21m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch21m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch21m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch21m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch21m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch21m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch21m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch21m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch21m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch21m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch21m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_omch21m, C_omch21m, a_data_omch21m, b_data_omch21m, c_data_omch21m, lps_cc_data_omch21m, lps_cs_data_omch21m, lps_ss_data_omch21m, scale_factor_data_omch21m, window_c_data_omch21m, window_s_data_omch21m, mps_data_omch21m, z_at_chi_data_omch21m = data_import_func('data_omch21m')

cpdef double lbs_omch21m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch21m, a_data_omch21m, b_data_omch21m, c_data_omch21m, scale_factor_data_omch21m, window_c_data_omch21m, window_s_data_omch21m, mps_data_omch21m, z_at_chi_data_omch21m)

cpdef double lps_omch21m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch21m, lps_cs_data_omch21m, lps_ss_data_omch21m)


cdef dict[str, double] cosm_par_omch22m
cdef double C_omch22m
cdef double[:, :] a_data_omch22m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch22m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch22m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch22m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch22m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch22m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch22m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch22m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch22m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch22m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch22m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_omch22m, C_omch22m, a_data_omch22m, b_data_omch22m, c_data_omch22m, lps_cc_data_omch22m, lps_cs_data_omch22m, lps_ss_data_omch22m, scale_factor_data_omch22m, window_c_data_omch22m, window_s_data_omch22m, mps_data_omch22m, z_at_chi_data_omch22m = data_import_func('data_omch22m')

cpdef double lbs_omch22m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch22m, a_data_omch22m, b_data_omch22m, c_data_omch22m, scale_factor_data_omch22m, window_c_data_omch22m, window_s_data_omch22m, mps_data_omch22m, z_at_chi_data_omch22m)

cpdef double lps_omch22m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch22m, lps_cs_data_omch22m, lps_ss_data_omch22m)


cdef dict[str, double] cosm_par_ns2p
cdef double C_ns2p
cdef double[:, :] a_data_ns2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns2p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_ns2p, C_ns2p, a_data_ns2p, b_data_ns2p, c_data_ns2p, lps_cc_data_ns2p, lps_cs_data_ns2p, lps_ss_data_ns2p, scale_factor_data_ns2p, window_c_data_ns2p, window_s_data_ns2p, mps_data_ns2p, z_at_chi_data_ns2p = data_import_func('data_ns2p')

cpdef double lbs_ns2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns2p, a_data_ns2p, b_data_ns2p, c_data_ns2p, scale_factor_data_ns2p, window_c_data_ns2p, window_s_data_ns2p, mps_data_ns2p, z_at_chi_data_ns2p)

cpdef double lps_ns2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns2p, lps_cs_data_ns2p, lps_ss_data_ns2p)


cdef dict[str, double] cosm_par_ns1p
cdef double C_ns1p
cdef double[:, :] a_data_ns1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns1p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_ns1p, C_ns1p, a_data_ns1p, b_data_ns1p, c_data_ns1p, lps_cc_data_ns1p, lps_cs_data_ns1p, lps_ss_data_ns1p, scale_factor_data_ns1p, window_c_data_ns1p, window_s_data_ns1p, mps_data_ns1p, z_at_chi_data_ns1p = data_import_func('data_ns1p')

cpdef double lbs_ns1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns1p, a_data_ns1p, b_data_ns1p, c_data_ns1p, scale_factor_data_ns1p, window_c_data_ns1p, window_s_data_ns1p, mps_data_ns1p, z_at_chi_data_ns1p)

cpdef double lps_ns1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns1p, lps_cs_data_ns1p, lps_ss_data_ns1p)


cdef dict[str, double] cosm_par_ns1m
cdef double C_ns1m
cdef double[:, :] a_data_ns1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns1m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_ns1m, C_ns1m, a_data_ns1m, b_data_ns1m, c_data_ns1m, lps_cc_data_ns1m, lps_cs_data_ns1m, lps_ss_data_ns1m, scale_factor_data_ns1m, window_c_data_ns1m, window_s_data_ns1m, mps_data_ns1m, z_at_chi_data_ns1m = data_import_func('data_ns1m')

cpdef double lbs_ns1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns1m, a_data_ns1m, b_data_ns1m, c_data_ns1m, scale_factor_data_ns1m, window_c_data_ns1m, window_s_data_ns1m, mps_data_ns1m, z_at_chi_data_ns1m)

cpdef double lps_ns1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns1m, lps_cs_data_ns1m, lps_ss_data_ns1m)


cdef dict[str, double] cosm_par_ns2m
cdef double C_ns2m
cdef double[:, :] a_data_ns2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns2m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_ns2m, C_ns2m, a_data_ns2m, b_data_ns2m, c_data_ns2m, lps_cc_data_ns2m, lps_cs_data_ns2m, lps_ss_data_ns2m, scale_factor_data_ns2m, window_c_data_ns2m, window_s_data_ns2m, mps_data_ns2m, z_at_chi_data_ns2m = data_import_func('data_ns2m')

cpdef double lbs_ns2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns2m, a_data_ns2m, b_data_ns2m, c_data_ns2m, scale_factor_data_ns2m, window_c_data_ns2m, window_s_data_ns2m, mps_data_ns2m, z_at_chi_data_ns2m)

cpdef double lps_ns2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns2m, lps_cs_data_ns2m, lps_ss_data_ns2m)


cdef dict[str, double] cosm_par_As2p
cdef double C_As2p
cdef double[:, :] a_data_As2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As2p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_As2p, C_As2p, a_data_As2p, b_data_As2p, c_data_As2p, lps_cc_data_As2p, lps_cs_data_As2p, lps_ss_data_As2p, scale_factor_data_As2p, window_c_data_As2p, window_s_data_As2p, mps_data_As2p, z_at_chi_data_As2p = data_import_func('data_As2p')

cpdef double lbs_As2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As2p, a_data_As2p, b_data_As2p, c_data_As2p, scale_factor_data_As2p, window_c_data_As2p, window_s_data_As2p, mps_data_As2p, z_at_chi_data_As2p)

cpdef double lps_As2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As2p, lps_cs_data_As2p, lps_ss_data_As2p)


cdef dict[str, double] cosm_par_As1p
cdef double C_As1p
cdef double[:, :] a_data_As1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As1p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_As1p, C_As1p, a_data_As1p, b_data_As1p, c_data_As1p, lps_cc_data_As1p, lps_cs_data_As1p, lps_ss_data_As1p, scale_factor_data_As1p, window_c_data_As1p, window_s_data_As1p, mps_data_As1p, z_at_chi_data_As1p = data_import_func('data_As1p')

cpdef double lbs_As1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As1p, a_data_As1p, b_data_As1p, c_data_As1p, scale_factor_data_As1p, window_c_data_As1p, window_s_data_As1p, mps_data_As1p, z_at_chi_data_As1p)

cpdef double lps_As1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As1p, lps_cs_data_As1p, lps_ss_data_As1p)


cdef dict[str, double] cosm_par_As1m
cdef double C_As1m
cdef double[:, :] a_data_As1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As1m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_As1m, C_As1m, a_data_As1m, b_data_As1m, c_data_As1m, lps_cc_data_As1m, lps_cs_data_As1m, lps_ss_data_As1m, scale_factor_data_As1m, window_c_data_As1m, window_s_data_As1m, mps_data_As1m, z_at_chi_data_As1m = data_import_func('data_As1m')

cpdef double lbs_As1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As1m, a_data_As1m, b_data_As1m, c_data_As1m, scale_factor_data_As1m, window_c_data_As1m, window_s_data_As1m, mps_data_As1m, z_at_chi_data_As1m)

cpdef double lps_As1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As1m, lps_cs_data_As1m, lps_ss_data_As1m)


cdef dict[str, double] cosm_par_As2m
cdef double C_As2m
cdef double[:, :] a_data_As2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As2m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_As2m, C_As2m, a_data_As2m, b_data_As2m, c_data_As2m, lps_cc_data_As2m, lps_cs_data_As2m, lps_ss_data_As2m, scale_factor_data_As2m, window_c_data_As2m, window_s_data_As2m, mps_data_As2m, z_at_chi_data_As2m = data_import_func('data_As2m')

cpdef double lbs_As2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As2m, a_data_As2m, b_data_As2m, c_data_As2m, scale_factor_data_As2m, window_c_data_As2m, window_s_data_As2m, mps_data_As2m, z_at_chi_data_As2m)

cpdef double lps_As2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As2m, lps_cs_data_As2m, lps_ss_data_As2m)


cdef dict[str, double] cosm_par_mnu2p
cdef double C_mnu2p
cdef double[:, :] a_data_mnu2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu2p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_mnu2p, C_mnu2p, a_data_mnu2p, b_data_mnu2p, c_data_mnu2p, lps_cc_data_mnu2p, lps_cs_data_mnu2p, lps_ss_data_mnu2p, scale_factor_data_mnu2p, window_c_data_mnu2p, window_s_data_mnu2p, mps_data_mnu2p, z_at_chi_data_mnu2p = data_import_func('data_mnu2p')

cpdef double lbs_mnu2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu2p, a_data_mnu2p, b_data_mnu2p, c_data_mnu2p, scale_factor_data_mnu2p, window_c_data_mnu2p, window_s_data_mnu2p, mps_data_mnu2p, z_at_chi_data_mnu2p)

cpdef double lps_mnu2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu2p, lps_cs_data_mnu2p, lps_ss_data_mnu2p)


cdef dict[str, double] cosm_par_mnu1p
cdef double C_mnu1p
cdef double[:, :] a_data_mnu1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu1p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_mnu1p, C_mnu1p, a_data_mnu1p, b_data_mnu1p, c_data_mnu1p, lps_cc_data_mnu1p, lps_cs_data_mnu1p, lps_ss_data_mnu1p, scale_factor_data_mnu1p, window_c_data_mnu1p, window_s_data_mnu1p, mps_data_mnu1p, z_at_chi_data_mnu1p = data_import_func('data_mnu1p')

cpdef double lbs_mnu1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu1p, a_data_mnu1p, b_data_mnu1p, c_data_mnu1p, scale_factor_data_mnu1p, window_c_data_mnu1p, window_s_data_mnu1p, mps_data_mnu1p, z_at_chi_data_mnu1p)

cpdef double lps_mnu1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu1p, lps_cs_data_mnu1p, lps_ss_data_mnu1p)


cdef dict[str, double] cosm_par_mnu1m
cdef double C_mnu1m
cdef double[:, :] a_data_mnu1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu1m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_mnu1m, C_mnu1m, a_data_mnu1m, b_data_mnu1m, c_data_mnu1m, lps_cc_data_mnu1m, lps_cs_data_mnu1m, lps_ss_data_mnu1m, scale_factor_data_mnu1m, window_c_data_mnu1m, window_s_data_mnu1m, mps_data_mnu1m, z_at_chi_data_mnu1m = data_import_func('data_mnu1m')

cpdef double lbs_mnu1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu1m, a_data_mnu1m, b_data_mnu1m, c_data_mnu1m, scale_factor_data_mnu1m, window_c_data_mnu1m, window_s_data_mnu1m, mps_data_mnu1m, z_at_chi_data_mnu1m)

cpdef double lps_mnu1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu1m, lps_cs_data_mnu1m, lps_ss_data_mnu1m)


cdef dict[str, double] cosm_par_mnu2m
cdef double C_mnu2m
cdef double[:, :] a_data_mnu2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu2m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_mnu2m, C_mnu2m, a_data_mnu2m, b_data_mnu2m, c_data_mnu2m, lps_cc_data_mnu2m, lps_cs_data_mnu2m, lps_ss_data_mnu2m, scale_factor_data_mnu2m, window_c_data_mnu2m, window_s_data_mnu2m, mps_data_mnu2m, z_at_chi_data_mnu2m = data_import_func('data_mnu2m')

cpdef double lbs_mnu2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu2m, a_data_mnu2m, b_data_mnu2m, c_data_mnu2m, scale_factor_data_mnu2m, window_c_data_mnu2m, window_s_data_mnu2m, mps_data_mnu2m, z_at_chi_data_mnu2m)

cpdef double lps_mnu2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu2m, lps_cs_data_mnu2m, lps_ss_data_mnu2m)


cdef dict[str, double] cosm_par_w02p
cdef double C_w02p
cdef double[:, :] a_data_w02p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w02p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w02p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w02p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w02p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w02p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w02p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w02p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w02p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w02p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w02p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_w02p, C_w02p, a_data_w02p, b_data_w02p, c_data_w02p, lps_cc_data_w02p, lps_cs_data_w02p, lps_ss_data_w02p, scale_factor_data_w02p, window_c_data_w02p, window_s_data_w02p, mps_data_w02p, z_at_chi_data_w02p = data_import_func('data_w02p')

cpdef double lbs_w02p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w02p, a_data_w02p, b_data_w02p, c_data_w02p, scale_factor_data_w02p, window_c_data_w02p, window_s_data_w02p, mps_data_w02p, z_at_chi_data_w02p)

cpdef double lps_w02p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w02p, lps_cs_data_w02p, lps_ss_data_w02p)


cdef dict[str, double] cosm_par_w01p
cdef double C_w01p
cdef double[:, :] a_data_w01p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w01p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w01p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w01p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w01p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w01p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w01p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w01p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w01p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w01p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w01p = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_w01p, C_w01p, a_data_w01p, b_data_w01p, c_data_w01p, lps_cc_data_w01p, lps_cs_data_w01p, lps_ss_data_w01p, scale_factor_data_w01p, window_c_data_w01p, window_s_data_w01p, mps_data_w01p, z_at_chi_data_w01p = data_import_func('data_w01p')

cpdef double lbs_w01p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w01p, a_data_w01p, b_data_w01p, c_data_w01p, scale_factor_data_w01p, window_c_data_w01p, window_s_data_w01p, mps_data_w01p, z_at_chi_data_w01p)

cpdef double lps_w01p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w01p, lps_cs_data_w01p, lps_ss_data_w01p)


cdef dict[str, double] cosm_par_w01m
cdef double C_w01m
cdef double[:, :] a_data_w01m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w01m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w01m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w01m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w01m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w01m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w01m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w01m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w01m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w01m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w01m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_w01m, C_w01m, a_data_w01m, b_data_w01m, c_data_w01m, lps_cc_data_w01m, lps_cs_data_w01m, lps_ss_data_w01m, scale_factor_data_w01m, window_c_data_w01m, window_s_data_w01m, mps_data_w01m, z_at_chi_data_w01m = data_import_func('data_w01m')

cpdef double lbs_w01m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w01m, a_data_w01m, b_data_w01m, c_data_w01m, scale_factor_data_w01m, window_c_data_w01m, window_s_data_w01m, mps_data_w01m, z_at_chi_data_w01m)

cpdef double lps_w01m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w01m, lps_cs_data_w01m, lps_ss_data_w01m)


cdef dict[str, double] cosm_par_w02m
cdef double C_w02m
cdef double[:, :] a_data_w02m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w02m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w02m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w02m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w02m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w02m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w02m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w02m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w02m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w02m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w02m = np.zeros(chi_num, dtype=np.float64)
#cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

cosm_par_w02m, C_w02m, a_data_w02m, b_data_w02m, c_data_w02m, lps_cc_data_w02m, lps_cs_data_w02m, lps_ss_data_w02m, scale_factor_data_w02m, window_c_data_w02m, window_s_data_w02m, mps_data_w02m, z_at_chi_data_w02m = data_import_func('data_w02m')

cpdef double lbs_w02m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w02m, a_data_w02m, b_data_w02m, c_data_w02m, scale_factor_data_w02m, window_c_data_w02m, window_s_data_w02m, mps_data_w02m, z_at_chi_data_w02m)

cpdef double lps_w02m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w02m, lps_cs_data_w02m, lps_ss_data_w02m)


# cdef dict[str, double] cosm_par_wa2p
# cdef double C_wa2p
# cdef double[:, :] a_data_wa2p = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:, :] b_data_wa2p = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:, :] c_data_wa2p = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:] lps_cc_data_wa2p = np.zeros(k_num, dtype=np.float64)
# cdef double[:] lps_cs_data_wa2p = np.zeros(k_num, dtype=np.float64)
# cdef double[:] lps_ss_data_wa2p = np.zeros(k_num, dtype=np.float64)
# cdef double[:] scale_factor_data_wa2p = np.zeros(chi_num, dtype=np.float64)
# cdef double[:] window_c_data_wa2p = np.zeros(chi_num, dtype=np.float64)
# cdef double[:] window_s_data_wa2p = np.zeros(chi_num, dtype=np.float64)
# cdef double[:, :] mps_data_wa2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
# cdef double[:] z_at_chi_data_wa2p = np.zeros(chi_num, dtype=np.float64)
# #cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

# cosm_par_wa2p, C_wa2p, a_data_wa2p, b_data_wa2p, c_data_wa2p, lps_cc_data_wa2p, lps_cs_data_wa2p, lps_ss_data_wa2p, scale_factor_data_wa2p, window_c_data_wa2p, window_s_data_wa2p, mps_data_wa2p, z_at_chi_data_wa2p = data_import_func('data_wa2p')

# cpdef double lbs_wa2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
#     return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_wa2p, a_data_wa2p, b_data_wa2p, c_data_wa2p, scale_factor_data_wa2p, window_c_data_wa2p, window_s_data_wa2p, mps_data_wa2p, z_at_chi_data_wa2p)

# cpdef double lps_wa2p(int l, char* type1, char* type2) noexcept nogil:
#     return lps(l, type1, type2, lps_cc_data_wa2p, lps_cs_data_wa2p, lps_ss_data_wa2p)


# cdef dict[str, double] cosm_par_wa1p
# cdef double C_wa1p
# cdef double[:, :] a_data_wa1p = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:, :] b_data_wa1p = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:, :] c_data_wa1p = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:] lps_cc_data_wa1p = np.zeros(k_num, dtype=np.float64)
# cdef double[:] lps_cs_data_wa1p = np.zeros(k_num, dtype=np.float64)
# cdef double[:] lps_ss_data_wa1p = np.zeros(k_num, dtype=np.float64)
# cdef double[:] scale_factor_data_wa1p = np.zeros(chi_num, dtype=np.float64)
# cdef double[:] window_c_data_wa1p = np.zeros(chi_num, dtype=np.float64)
# cdef double[:] window_s_data_wa1p = np.zeros(chi_num, dtype=np.float64)
# cdef double[:, :] mps_data_wa1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
# cdef double[:] z_at_chi_data_wa1p = np.zeros(chi_num, dtype=np.float64)
# #cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

# cosm_par_wa1p, C_wa1p, a_data_wa1p, b_data_wa1p, c_data_wa1p, lps_cc_data_wa1p, lps_cs_data_wa1p, lps_ss_data_wa1p, scale_factor_data_wa1p, window_c_data_wa1p, window_s_data_wa1p, mps_data_wa1p, z_at_chi_data_wa1p = data_import_func('data_wa1p')

# cpdef double lbs_wa1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
#     return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_wa1p, a_data_wa1p, b_data_wa1p, c_data_wa1p, scale_factor_data_wa1p, window_c_data_wa1p, window_s_data_wa1p, mps_data_wa1p, z_at_chi_data_wa1p)

# cpdef double lps_wa1p(int l, char* type1, char* type2) noexcept nogil:
#     return lps(l, type1, type2, lps_cc_data_wa1p, lps_cs_data_wa1p, lps_ss_data_wa1p)


# cdef dict[str, double] cosm_par_wa1m
# cdef double C_wa1m
# cdef double[:, :] a_data_wa1m = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:, :] b_data_wa1m = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:, :] c_data_wa1m = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:] lps_cc_data_wa1m = np.zeros(k_num, dtype=np.float64)
# cdef double[:] lps_cs_data_wa1m = np.zeros(k_num, dtype=np.float64)
# cdef double[:] lps_ss_data_wa1m = np.zeros(k_num, dtype=np.float64)
# cdef double[:] scale_factor_data_wa1m = np.zeros(chi_num, dtype=np.float64)
# cdef double[:] window_c_data_wa1m = np.zeros(chi_num, dtype=np.float64)
# cdef double[:] window_s_data_wa1m = np.zeros(chi_num, dtype=np.float64)
# cdef double[:, :] mps_data_wa1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
# cdef double[:] z_at_chi_data_wa1m = np.zeros(chi_num, dtype=np.float64)
# #cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

# cosm_par_wa1m, C_wa1m, a_data_wa1m, b_data_wa1m, c_data_wa1m, lps_cc_data_wa1m, lps_cs_data_wa1m, lps_ss_data_wa1m, scale_factor_data_wa1m, window_c_data_wa1m, window_s_data_wa1m, mps_data_wa1m, z_at_chi_data_wa1m = data_import_func('data_wa1m')

# cpdef double lbs_wa1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
#     return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_wa1m, a_data_wa1m, b_data_wa1m, c_data_wa1m, scale_factor_data_wa1m, window_c_data_wa1m, window_s_data_wa1m, mps_data_wa1m, z_at_chi_data_wa1m)

# cpdef double lps_wa1m(int l, char* type1, char* type2) noexcept nogil:
#     return lps(l, type1, type2, lps_cc_data_wa1m, lps_cs_data_wa1m, lps_ss_data_wa1m)


# cdef dict[str, double] cosm_par_wa2m
# cdef double C_wa2m
# cdef double[:, :] a_data_wa2m = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:, :] b_data_wa2m = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:, :] c_data_wa2m = np.zeros((z_num, k_num), dtype=np.float64)
# cdef double[:] lps_cc_data_wa2m = np.zeros(k_num, dtype=np.float64)
# cdef double[:] lps_cs_data_wa2m = np.zeros(k_num, dtype=np.float64)
# cdef double[:] lps_ss_data_wa2m = np.zeros(k_num, dtype=np.float64)
# cdef double[:] scale_factor_data_wa2m = np.zeros(chi_num, dtype=np.float64)
# cdef double[:] window_c_data_wa2m = np.zeros(chi_num, dtype=np.float64)
# cdef double[:] window_s_data_wa2m = np.zeros(chi_num, dtype=np.float64)
# cdef double[:, :] mps_data_wa2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
# cdef double[:] z_at_chi_data_wa2m = np.zeros(chi_num, dtype=np.float64)
# #cosm_par, C, a_data, b_data, c_data, lps_cc_data, lps_cs_data, lps_ss_data, scale_factor_data, window_c_data, window_s_data, mps_data, z_at_chi_data

# cosm_par_wa2m, C_wa2m, a_data_wa2m, b_data_wa2m, c_data_wa2m, lps_cc_data_wa2m, lps_cs_data_wa2m, lps_ss_data_wa2m, scale_factor_data_wa2m, window_c_data_wa2m, window_s_data_wa2m, mps_data_wa2m, z_at_chi_data_wa2m = data_import_func('data_wa2m')

# cpdef double lbs_wa2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
#     return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_wa2m, a_data_wa2m, b_data_wa2m, c_data_wa2m, scale_factor_data_wa2m, window_c_data_wa2m, window_s_data_wa2m, mps_data_wa2m, z_at_chi_data_wa2m)

# cpdef double lps_wa2m(int l, char* type1, char* type2) noexcept nogil:
#     return lps(l, type1, type2, lps_cc_data_wa2m, lps_cs_data_wa2m, lps_ss_data_wa2m)

##########################################


cdef double der_2o(double f2p, double f1p, double f1m, double f2m, double dx) noexcept nogil:
    return (-1 * f2p + 8 * f1p - 8 * f1m + f2m) / (12 * dx)

cdef double dx_frac = 0.025
cdef double[:] fiducial_cosm_par = np.array([67.4, 0.0224, 0.120, 0.965, 2.1e-9, 0.06])


cdef double lps_der(int k, char* type1, char* type2, char* par) noexcept nogil:
    # for b'snr' case, just returns the original function
    if par[0] == b's':
        return lps_f(k, type1, type2)

    if par[0] == b'H':
        return der_2o(lps_H2p(k, type1, type2), lps_H1p(k, type1, type2), lps_H1m(k, type1, type2), lps_H2m(k, type1, type2), fiducial_cosm_par[0] * dx_frac)
            
    if par[0] == b'o' and par[2] == b'b':
        return der_2o(lps_ombh22p(k, type1, type2), lps_ombh21p(k, type1, type2), lps_ombh21m(k, type1, type2), lps_ombh22m(k, type1, type2), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'o' and par[2] == b'c':
        return der_2o(lps_omch22p(k, type1, type2), lps_omch21p(k, type1, type2), lps_omch21m(k, type1, type2), lps_omch22m(k, type1, type2), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'A':
        return der_2o(lps_As2p(k, type1, type2), lps_As1p(k, type1, type2), lps_As1m(k, type1, type2), lps_As2m(k, type1, type2), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'n':
        return der_2o(lps_ns2p(k, type1, type2), lps_ns1p(k, type1, type2), lps_ns1m(k, type1, type2), lps_ns2m(k, type1, type2), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'm':
        return der_2o(lps_mnu2p(k, type1, type2), lps_mnu1p(k, type1, type2), lps_mnu1m(k, type1, type2), lps_mnu2m(k, type1, type2), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'w' and par[1] == b'0':
        return der_2o(lps_w02p(k, type1, type2), lps_w01p(k, type1, type2), lps_w01m(k, type1, type2), lps_w02m(k, type1, type2), fiducial_cosm_par[0] * dx_frac)

    # if par[0] == b'w' and par[1] == b'a':
    #     return der_2o(lps_wa2p(k, type1, type2), lps_wa1p(k, type1, type2), lps_wa1m(k, type1, type2), lps_wa2m(k, type1, type2), fiducial_cosm_par[0] * dx_frac)


cdef double lbs_der(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, char* par) noexcept nogil:
    # for b'snr' case, just returns the original function
    if par[0] == b's':
        return lbs_f(k1, k2, k3, type1, type2, type3, num_samples)

    if par[0] == b'H':
        return der_2o(lbs_H2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_H1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_H1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_H2m(k1, k2, k3, type1, type2, type3, num_samples), fiducial_cosm_par[0] * dx_frac)
            
    if par[0] == b'o' and par[2] == b'b':
        return der_2o(lbs_ombh22p(k1, k2, k3, type1, type2, type3, num_samples), lbs_ombh21p(k1, k2, k3, type1, type2, type3, num_samples), lbs_ombh21m(k1, k2, k3, type1, type2, type3, num_samples), lbs_ombh22m(k1, k2, k3, type1, type2, type3, num_samples), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'o' and par[2] == b'c':
        return der_2o(lbs_omch22p(k1, k2, k3, type1, type2, type3, num_samples), lbs_omch21p(k1, k2, k3, type1, type2, type3, num_samples), lbs_omch21m(k1, k2, k3, type1, type2, type3, num_samples), lbs_omch22m(k1, k2, k3, type1, type2, type3, num_samples), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'A':
        return der_2o(lbs_As2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_As1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_As1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_As2m(k1, k2, k3, type1, type2, type3, num_samples), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'n':
        return der_2o(lbs_ns2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_ns1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_ns1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_ns2m(k1, k2, k3, type1, type2, type3, num_samples), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'm':
        return der_2o(lbs_mnu2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_mnu1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_mnu1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_mnu2m(k1, k2, k3, type1, type2, type3, num_samples), fiducial_cosm_par[0] * dx_frac)

    if par[0] == b'w' and par[1] == b'0':
        return der_2o(lbs_w02p(k1, k2, k3, type1, type2, type3, num_samples), lbs_w01p(k1, k2, k3, type1, type2, type3, num_samples), lbs_w01m(k1, k2, k3, type1, type2, type3, num_samples), lbs_w02m(k1, k2, k3, type1, type2, type3, num_samples), fiducial_cosm_par[0] * dx_frac)

    # if par[0] == b'w' and par[1] == b'a':
    #     return der_2o(lbs_wa2p(k1, k2, k3, type1, type2, type3, num_samples), lbs_wa1p(k1, k2, k3, type1, type2, type3, num_samples), lbs_wa1m(k1, k2, k3, type1, type2, type3, num_samples), lbs_wa2m(k1, k2, k3, type1, type2, type3, num_samples), fiducial_cosm_par[0] * dx_frac)

'''

cdef double lbs_der(int l1, int l2, int l3, char* type1, char* type2, char* type3, int num_samples, char* par) noexcept nogil:
    cdef double wigner_factor
    cdef double sqrt_factor
    if full_sky:
        if (l1 + l2 + l3) % 2 == 0 and l3 <= l1 + l2 and l1 - l2 <= l3 and l2 - l1 <= l3:
            wigner_factor = wigner_3j_approx_nocheck(l1, l2, l3)
            sqrt_factor = sqrt((2*l1 + 1)*(2*l2 + 1)*(2*l3 + 1)/(4 * 3.14159)) # (pi)
            return wigner_factor * sqrt_factor * lbs_der_flat(l1, l2, l3, type1, type2, type3, num_samples, par)
        else:
            return 0.
    else:
        return lbs_der_flat(l1, l2, l3, type1, type2, type3, num_samples, par)

# from here stuff is old


cdef extern from "math.h":
    double sqrt(double x) nogil
    double log(double x) nogil
    double exp(double x) nogil

cdef extern from "complex.h":
    double creal(double complex x) nogil

# factorial function got too big in wigner3j function and caused over/underflow problems
cdef double factorial(int n) nogil:
    cdef long result = 1
    cdef int i
    for i in range(2, n+1):
        result *= i
    return result

# gives e log of factorial instead
cdef double log_factorial(int n) nogil:
    cdef double result = 0
    cdef int i
    for i in range(2, n+1):
        result += log(i)
    return result

# exact version of wigner_3j symbol
# cdef double wigner_3j(int l1, int l2, int l3) noexcept nogil:
#     cdef int L = l1 + l2 + l3
#     if L % 2 != 0:
#         return 0.0  # The result is zero for odd L.

#     if L > 60:
#         return 0.0

#     cdef int L_half = L // 2

#     # Numerator
#     cdef double numerator = (-1) ** L_half

#     # Cancelling factorial terms:
#     # The factorials in the denominator are divided by the numerator factorial to reduce overflow.
    
#     cdef double denom1 = 1.0
#     cdef int k = 1

#     cdef int i
    
#     # Calculate the denominator terms while cancelling as much as possible with the numerator
#     for i in range(1, L_half - l1 + 1):
#         denom1 *= i

#     for i in range(1, L_half - l2 + 1):
#         denom1 *= i
    
#     for i in range(1, L_half - l3 + 1):
#         denom1 *= i

#     # Square root term with cancelling
#     cdef double sqrt_term = 1.0
#     for i in range(L_half - l1 + 1, L - 2 * l1 + 1):
#         sqrt_term *= i

#     for i in range(L_half - l2 + 1, L - 2 * l2 + 1):
#         sqrt_term *= i

#     for i in range(L_half - l3 + 1, L - 2 * l3 + 1):
#         sqrt_term *= i

#     return numerator / denom1 * sqrt(exp(log(sqrt_term) - log_factorial(L + 1)))

# approximate expression of wigner_3j symbol based on appendix in paper "Cosmological parameters from lensing power spectrum and bispectrum tomography"
cdef double wigner_3j_approx_nocheck(int l1, int l2, int l3) noexcept nogil:
    cdef double complex L_half, factor, term1, term2, term3, term1_pow, term2_pow, term3_pow, denominator
    cdef int L = l1 + l2 + l3

    # does not check if L is even or triangle inequalities, does so in lbs_f and lbs_der directly instead to save on computing bispec if result should be zero anyway

    L_half = L / 2.0
    
    # Common factors
    factor = (-1)**L_half * (2 * 3.141592653589793)**(-0.5) * exp(3.0 / 2) * (L + 2)**(-0.25)
    
    # Power term for each fraction
    term1 = (L_half - l1 + 0.5) / (L_half - l1 + 1)
    term2 = (L_half - l2 + 0.5) / (L_half - l2 + 1)
    term3 = (L_half - l3 + 0.5) / (L_half - l3 + 1)

    # Raising to required powers
    term1_pow = term1 ** (L_half - l1 + 1.0 / 4.0)
    term2_pow = term2 ** (L_half - l2 + 1.0 / 4.0)
    term3_pow = term3 ** (L_half - l3 + 1.0 / 4.0)

    # The denominator terms
    denominator = ((L_half - l1 + 1)**0.25 * (L_half - l2 + 1)**0.25 * (L_half - l3 + 1)**0.25)
    
    # The final result
    return creal(factor * term1_pow * term2_pow * term3_pow / denominator)

# fiducial (!) full sky lensing bispectrum
cpdef double lbs(int l1, int l2, int l3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    cdef double wigner_factor
    cdef double sqrt_factor
    if (l1 + l2 + l3) % 2 == 0 and l3 <= l1 + l2 and l1 - l2 <= l3 and l2 - l1 <= l3:
        wigner_factor = wigner_3j_approx_nocheck(l1, l2, l3)
        sqrt_factor = sqrt((2.0*l1 + 1.0)*(2.0*l2 + 1.0)*(2.0*l3 + 1.0)/(4 * 3.14159)) # (pi)
        return wigner_factor * sqrt_factor * lbs_f(l1, l2, l3, type1, type2, type3, num_samples)
    else:
        return 0.0

'''

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

def lps_der_test(k, type1, type2, par):
    return lps_der(k, type1, type2, par)

def lbs_der_test(k1, k2, k3, type1, type2, type3, num_samples,  par):
    return lbs_der(k1, k2, k3, type1, type2, type3, num_samples,  par)


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
