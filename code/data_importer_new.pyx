import json
cimport interpolation as interp
import math
import numpy as np
from libc.stdio cimport printf
from libc.math cimport fmax, fabs, pow

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

def set_noise_types(cmb_type, galaxy_type):
    '''
    CMB type values correspond to:
    0: S0 noise curves (old, remove)
    1: sigma = 1, Delta_P = 6 (stage 3 toshiya)
    2: sigma = 3, Delta_P = 1 (stage 4 toshiya)
    3: sigma = 5, Delta_T = 30, Delta_P = 52 (planck (double check))
    Galaxy type values correspond to:
    1: 5 arcmin^-2 galaxy density, eccentricity 0.3 (stage 3 like)
    2: 30 arcmin^-2 galaxy density, eccentricity 0.3 (stage 4 like)
    3: 100 arcmin^-2 galaxy density, eccentricity 0.4 (optimistic, from https://arxiv.org/pdf/astro-ph/0310125)
    '''
    global cmb_noise_type, galaxy_noise_type
    cmb_noise_type = cmb_type
    galaxy_noise_type = galaxy_type
    print(f'CMB noise type: {cmb_noise_type}\nGalaxy noise type: {galaxy_noise_type}')
    pass

print(f'CMB noise type: {cmb_noise_type}\nGalaxy noise type: {galaxy_noise_type}')

# (units are in microKelvin arcmin and arcmin)

cdef int lmax_cmbps = 3000

# SO noise
cdef str folder_file_path = '/scratch/p319950/data/'
# cdef str filepath_convergence_noise_file_path = '/scratch/p319950/data_rough/' + 'conv_noise.dat'
# conv_noise_data_array = np.loadtxt(filepath_convergence_noise_file_path)
# cdef double[:, :] conv_noise_data = conv_noise_data_array
# conv noise from quadratic estimator
cdef double[:] ls_cmbn = np.loadtxt('cmb_noise_files/ls_1_3000_64.txt')
ls_cmbn_np_array = np.loadtxt('cmb_noise_files/ls_1_3000_64.txt') # need to also have a np array version to easily convert noise values from convergence to lens potential below
cdef double lmin_cmbn = ls_cmbn[0]
cdef int lnum_cmbn = len(ls_cmbn)
cdef double lmax_cmbn = ls_cmbn[lnum_cmbn - 1]

cdef double[:] cmbn_301 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma3_DeltaT0.71_DeltaP1.txt')) / (ls_cmbn_np_array * (ls_cmbn_np_array + 1))
cdef double[:] cmbn_106 = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma1_DeltaT4.2_DeltaP6.txt')) / (ls_cmbn_np_array * (ls_cmbn_np_array + 1))
cdef double[:] planck_noise = np.abs(np.loadtxt('cmb_noise_files/Ns_sigma_planck.txt')) / (ls_cmbn_np_array * (ls_cmbn_np_array + 1))

toshiya_derivatives = True
# perhaps we can export data as c arrays instead of as memory views in the future so that we can specify return types like below and can avoid having to 
# declare types of all data before every data import
#cdef (dict[str, double], double, double[:, :], double[:, :], double[:], double[:], double[:], double[:], double[:], double[:], double[:, :], double[:]) data_importer(str folder_name):
def data_import_func(folder_name):
    cdef str filepath
    if toshiya_derivatives:# and (folder_name[-1] == '0' or folder_name == 'data_fiducial'):
        filepath = '/scratch/p319950/data_toshiya_like/' + folder_name

    else:
        filepath = folder_file_path + folder_name
    print(filepath)

    cdef double[:] cosm_par = np.load(filepath + '/cosm_par.npy')
    cdef double C
    cdef double[:, :] a_data = np.load(filepath + '/a.npy')
    cdef double[:, :] b_data = np.load(filepath + '/b.npy')
    cdef double[:, :] c_data = np.load(filepath + '/c.npy')

    cdef double[:, :] lps_data = np.load(filepath + '/lensing_power_spectrum.npy')

    cdef double[:] scale_factor_data = np.load(filepath + '/scale_factor.npy')
    cdef double[:, :] window_data = np.load(filepath + '/window_func.npy')
    cdef double[:, :] mps_data = np.load(filepath + '/matter_power_spectrum.npy')
    cdef double[:] z_at_chi_data = np.load(filepath + '/z_at_chi.npy')
    cdef double[:, :] cmbps = np.load(filepath + '/cmb_ps_with_ls.npy')
    cdef double[:, :] galaxy_density_chi_bins = np.load(filepath + '/galaxy_density_chi_bins.npy')

    # all 2d funcs are made such that order of arguments is (k, z)
    C = 0.
    with open(filepath + '/rho_bar', 'r') as file:
        C = float(file.read())

    return cosm_par, C, a_data, b_data, c_data, lps_data, scale_factor_data, window_data, mps_data, z_at_chi_data, cmbps, galaxy_density_chi_bins



# turns data into functions 

cpdef double a(double k, double z, double[:, :] a_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, a_data)
cpdef double b(double k, double z, double[:, :] b_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, b_data)
cpdef double c(double k, double z, double[:, :] c_data) noexcept nogil:
    return interp.logspace_linear_interp2d(k, z, k_min, k_max, k_num_fine, z_min, z_max, z_num_fine, c_data)

cdef int lps_tracer_pair_index(int idx1, int idx2) noexcept nogil:
    # Integer-indexed version: (0,0) (0,1) etc., where 0=c, 1=s1, 2=s2, 3=s3, 4=s4
    # Mapping: pairs ordered as [('c', 'c'), ('c', 's1'), ('c', 's2'), ('c', 's3'), ('c', 's4'), ('s1', 's1'), ('s1', 's2'), ('s1', 's3'), ('s1', 's4'), ('s2', 's2'), ('s2', 's3'), ('s2', 's4'), ('s3', 's3'), ('s3', 's4'), ('s4', 's4')]
    
    cdef int i1 = idx1
    cdef int i2 = idx2
    
    # Canonical ordering: c comes first, then s1, s2, s3, s4 in increasing order
    if i1 > i2:
        i1, i2 = i2, i1
    
    # Map to index
    if i1 == 0 and i2 == 0:
        return 0
    elif i1 == 0 and i2 == 1:
        return 1
    elif i1 == 0 and i2 == 2:
        return 2
    elif i1 == 0 and i2 == 3:
        return 3
    elif i1 == 0 and i2 == 4:
        return 4
    elif i1 == 1 and i2 == 1:
        return 5
    elif i1 == 1 and i2 == 2:
        return 6
    elif i1 == 1 and i2 == 3:
        return 7
    elif i1 == 1 and i2 == 4:
        return 8
    elif i1 == 2 and i2 == 2:
        return 9
    elif i1 == 2 and i2 == 3:
        return 10
    elif i1 == 2 and i2 == 4:
        return 11
    elif i1 == 3 and i2 == 3:
        return 12
    elif i1 == 3 and i2 == 4:
        return 13
    elif i1 == 4 and i2 == 4:
        return 14
    else:
        printf("Error: invalid indices in lps_tracer_pair_index\n")
        return -1

# without gil cython doesn't like strings bc they are python objects, you can use char* C objects (C like strings)
cdef double lps(double l, int type1, int type2, double[:, :] lps_data, double[:, :] cmbps) noexcept nogil:
    # CAUTION: this function doesn't throw error messages if type1 or type2 are not of type convergence or shear
    # for optimization purposes
    cdef int index
    # non lensing stuff
    if type1 == 5 and type2 == 5:  # 't' (temperature)
        return interp.linear_interp(l, cmbps[0, 0], lmax_cmbps, lmax_cmbps + 1, cmbps[:, 1])
    elif type1 == 6 and type2 == 6:  # 'e' (E-mode)
        return interp.linear_interp(l, cmbps[0, 0], lmax_cmbps, lmax_cmbps + 1, cmbps[:, 2])
    elif type1 == 7 and type2 == 7:  # 'b' (B-mode)
        return interp.linear_interp(l, cmbps[0, 0], lmax_cmbps, lmax_cmbps + 1, cmbps[:, 3])
    elif (type1 == 5 and type2 == 6) or (type1 == 6 and type2 == 5):  # 't' and 'e'
        return interp.linear_interp(l, cmbps[0, 0], lmax_cmbps, lmax_cmbps + 1, cmbps[:, 4])    
    else:
        index = lps_tracer_pair_index(type1, type2)
        return interp.logspace_linear_interp(l, k_min, k_max, k_num_fine, lps_data[index, :])

cdef double scale_factor(double chi, double[:] scale_factor_data) noexcept nogil:    
    return interp.linear_interp(chi, chi_min, chi_max, chi_num, scale_factor_data)


cdef double window_func(double chi, int type, double[:, :] window_data) noexcept nogil:
    if type == 0 or type == 1 or type == 2 or type == 3 or type == 4:
        return fmax(0., interp.logspace_linear_interp(chi, chi_min, chi_max, chi_num, window_data[type, :]))
    else:
        printf("Error: invalid type in window_func\n")
        return -1.

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

############

# post born corrections, based on https://arxiv.org/pdf/1605.05662

# old version for total gal population, not split into bins
# cdef double galaxy_density(double chi, double[:] z_at_chi_data) noexcept nogil: # TODO: adjust this for bin number
#     cdef double z = z_at_chi(chi, z_at_chi_data)
#     # because p(z) is a probability distribution, if you want it in chi coordinates instead of z coordinates, you
#     # need to take p(z(chi)) * dz/dchi. I have checked that the following simple approximation for the derivative works quite well
#     cdef double dzdchi = (z_at_chi(chi + 0.1, z_at_chi_data) - z_at_chi(chi, z_at_chi_data)) / 0.1
#     return dzdchi * pow(z, 2) * 0.64 * exp(-pow(z/0.64, 1.5)) / 0.111848  # TODO: check that exp here works, last factor is for normalization

cdef double galaxy_density(double chi, int bin_number, double[:, :] galaxy_density_chi_bins) noexcept nogil:
    return interp.linear_interp(chi, chi_min, chi_max, chi_num, galaxy_density_chi_bins[bin_number - 1, :])

cdef double pb_window_func(
    double chi,
    double chi_s,
    int source_type,
    double[:, :] galaxy_density_chi_bins_data
) noexcept nogil:
    cdef int N = 20 # TODO: convergence check in window_func_convergence.png, N=20 seems to be fine to a percent or so
    cdef double dchi, chi_prime
    cdef double sum = 0.0
    cdef double weight, kernel
    cdef int i

    if source_type == 0: # source_type is 0 (convergence)
        if chi < chi_s:
            return (chi_s - chi) / (chi * chi_s)
        else:
            return 0.

    else: # TODO: make warning signal for if source type is not recognized
        if chi >= chi_s:
            return 0.

        dchi = (chi_s - chi) / N
        for i in range(N + 1):
            chi_prime = chi + i * dchi
            weight = 0.5 if (i == 0 or i == N) else 1.0
            kernel = (chi_prime - chi) / (chi_prime * chi)
            sum += weight * galaxy_density(chi_prime, source_type, galaxy_density_chi_bins_data) * kernel

        return dchi * sum

# version from 21/06/2025
# cdef double lbs_pb_lps_integrand(
#     double chi,
#     double chi_s, # chi_s should be chi of last scattering in the M_s formula because that is the one associated with source 2
#     double chi_s_prime, 
#     double l, 
#     double C_data, 
#     double [:,:] mps_data, 
#     double [:] scale_factor_data,
#     double [:] z_at_chi_data,
#     char* source_type_2,
#     char* source_type_3
#     ) noexcept nogil:
#     # C = 3 * self.omega_m * H0**2 / (2 * self.lightspeed_kms**2)
#     if chi < chi_s and chi < chi_s_prime:
#         return 0.25 * 4 * C_data**2 * chi**2 * scale_factor(chi, scale_factor_data)**(-2) * pb_window_func(chi, chi_s, source_type_2, z_at_chi_data) * pb_window_func(chi, chi_s_prime, source_type_3, z_at_chi_data) * matter_power_spectrum(l / chi, z_at_chi(chi, z_at_chi_data), mps_data)
#     else:
#         return 0.

# new version
cdef double lbs_pb_lps_integrand(
    double chi,
    double chi_s, # chi_s should be chi of last scattering in the M_s formula because that is the one associated with source 2
    double chi_s_prime, 
    double l, 
    double C_data, 
    double [:,:] mps_data, 
    double [:] scale_factor_data,
    double [:] z_at_chi_data,
    double [:, :] galaxy_density_chi_bins_data,
    int source_type_2,
    int source_type_3
    ) noexcept nogil:
    # C = 3 * self.omega_m * H0**2 / (2 * self.lightspeed_kms**2)
    if chi < chi_s and chi < chi_s_prime:
        return 0.25 * 4 * C_data**2 * chi**2 * scale_factor(chi, scale_factor_data)**(-2) * pb_window_func(chi, chi_s, source_type_2, galaxy_density_chi_bins_data) * (chi_s_prime - chi) / (chi_s_prime * chi) * matter_power_spectrum(l / chi, z_at_chi(chi, z_at_chi_data), mps_data)
    else:
        return 0.

cpdef double lbs_pb_lps(
    double chi_s,
    double chi_s_prime, 
    double l,
    double C_data, 
    double [:,:] mps_data, 
    double [:] scale_factor_data,
    double [:] z_at_chi_data,
    double [:, :] galaxy_density_chi_bins_data,
    int num_samples,
    int source_type_2, # starting from 2 is to follow notation in the paper
    int source_type_3
    ) noexcept nogil:
    cdef double int_width = (chi_max - chi_min) / (num_samples - 1)
    cdef double result = 0
    cdef int i
    for i in range(1, num_samples - 1):
        result += lbs_pb_lps_integrand(chi_min + int_width * i, chi_s, chi_s_prime, l, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, source_type_2, source_type_3)
        # for boundary values use values that lie *just* inside the boundaries to prevent some nasty errors
    result += 0.5 * (
        lbs_pb_lps_integrand(chi_min * 1.01, chi_s, chi_s_prime, l, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, source_type_2, source_type_3) 
        + lbs_pb_lps_integrand(chi_max - 1.0, chi_s, chi_s_prime, l, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, source_type_2, source_type_3)
        )
    result *= int_width

    return result

cdef double lbs_pb_ms_integrand(
    double chi, # this is the integration variable, called chi_1 in my notes
    double l, 
    double l_prime, 
    double C_data, 
    double [:,:] mps_data, 
    double [:] scale_factor_data,
    double [:] z_at_chi_data,
    double [:, :] galaxy_density_chi_bins_data,
    double [:, :] window_data,
    int num_samples,
    int source_type_1,
    int source_type_2,
    int source_type_3
    ) noexcept nogil:

    # here l is l_1, l_prime is l_2, l_3 is not used (explicitly, at least)

    cdef double chi_last_scatter = 13912 # according to fiducial model

    return 0.25 * 4 * C_data**2 * chi**2 * scale_factor(chi, scale_factor_data)**(-2) * pb_window_func(chi, chi_last_scatter, source_type_1, galaxy_density_chi_bins_data) * pb_window_func(chi, chi_last_scatter, source_type_3, galaxy_density_chi_bins_data) * matter_power_spectrum(l / chi, z_at_chi(chi, z_at_chi_data), mps_data) * lbs_pb_lps(chi_last_scatter, chi, l_prime, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, num_samples, source_type_2, source_type_3)

cpdef double lbs_pb_ms(
    double l,
    double l_prime, 
    double C_data, 
    double [:,:] mps_data, 
    double [:] scale_factor_data,
    double [:] z_at_chi_data,
    double [:, :] galaxy_density_chi_bins_data,
    double [:, :] window_data,
    int num_samples,
    int type1,
    int type2,
    int type3
    ) noexcept nogil:
    cdef double int_width = (chi_max - chi_min) / (num_samples - 1)
    cdef double result = 0
    cdef int i
    for i in range(1, num_samples - 1):
        result += lbs_pb_ms_integrand(chi_min + int_width * i, l, l_prime, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, window_data, num_samples, type1, type2, type3)
        # for boundary values use values that lie *just* inside the boundaries to prevent some nasty errors
    result += 0.5 * (
        lbs_pb_ms_integrand(chi_min * 1.01, l, l_prime, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, window_data, num_samples, type1, type2, type3)
        + lbs_pb_ms_integrand(chi_max - 1.0, l, l_prime, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, window_data, num_samples, type1, type2, type3)
        )
    result *= int_width

    return result

cdef double lbs_pb_term(
    double l1, 
    double l2, 
    double l3, 
    double C_data, 
    double [:,:] mps_data, 
    double [:] scale_factor_data,
    double [:] z_at_chi_data,
    double [:, :] galaxy_density_chi_bins_data,
    double [:, :] window_data,
    int num_samples,
    int type1,
    int type2,
    int type3
    ) noexcept nogil:
    # cdef double law_cosines(double x, double y, double z) noexcept nogil:
    # gives cosine of angle between vector x and y, where we know the magnitudes of x, y, z and that x + y + z = 0 vector
    # return -1 * (x**2 + y**2 - z**2) / (2 * x * y)

    return 2 * law_cosines(l1, l2, l3) * l1**(-1) * l2**(-1) * (law_cosines(l1, l3, l2) * l1 * l3 * lbs_pb_ms(l1, l2, C_data, mps_data,scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, window_data, num_samples, type1, type2, type3) + law_cosines(l2, l3, l1) * l2 * l3 * lbs_pb_ms(l2, l1, C_data, mps_data,scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, window_data, num_samples, type2, type1, type3))

cdef double lbs_pb_flat(
    double l1, 
    double l2, 
    double l3, 
    double C_data, 
    double [:,:] mps_data, 
    double [:] scale_factor_data,
    double [:] z_at_chi_data,
    double [:, :] galaxy_density_chi_bins_data,
    double [:, :] window_data,
    int num_samples,
    int type1,
    int type2,
    int type3
    ) noexcept nogil:
    return lbs_pb_term(l1, l2, l3, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, window_data, num_samples, type1, type2, type3) + lbs_pb_term(l2, l3, l1, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, window_data, num_samples, type2, type3, type1) + lbs_pb_term(l3, l1, l2, C_data, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, window_data, num_samples, type3, type1, type2)

############

cdef inline double lbs_integrand(
    double chi, 
    double k1, 
    double k2, 
    double k3, 
    int type1, 
    int type2, 
    int type3,
    double[:, :] a_data, 
    double[:, :] b_data, 
    double[:, :] c_data, 
    double[:] scale_factor_data, 
    double[:, :] window_data, 
    double[:, :] mps_data, 
    double[:] z_at_chi_data
    ) noexcept nogil:

    cdef double integrand_val
    cdef double factor1, factor2, factor3, factor4, factor5, factor6

    # Compute individual factors
    factor1 = chi**2
    factor2 = scale_factor(chi, scale_factor_data)**(-3)
    factor3 = window_func(chi, type1, window_data)
    factor4 = window_func(chi, type2, window_data)
    factor5 = window_func(chi, type3, window_data)
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

# cdef inline double lbs_integrand(double chi, double k1, double k2, double k3, char* type1, char* type2, char* type3, double[:, :] a_data, double[:, :] b_data, double[:, :] c_data, double[:] scale_factor_data, double[:, :] window_data, double[:, :] mps_data, double[:] z_at_chi_data) noexcept nogil:
#     cdef double integrand_val
#     integrand_val = chi**2 * scale_factor(chi, scale_factor_data)**(-3) * window_func(chi, type1, window_data) * window_func(chi, type2, window_data) * window_func(chi, type3, window_data) * mbs(k1 / chi, k2 / chi, k3 / chi, z_at_chi(chi, z_at_chi_data), mps_data, a_data, b_data, c_data)
#     if integrand_val < 0.:

#     return integrand_val

cdef double lbs_flat(
    double k1,
    double k2,
    double k3,
    int type1,
    int type2,
    int type3,
    int num_samples,
    bint pb_correction,
    double C,
    double[:, :] a_data,
    double[:, :] b_data,
    double[:, :] c_data,
    double[:] scale_factor_data,
    double[:, :] window_data,
    double[:, :] mps_data,
    double[:] z_at_chi_data,
    double[:, :] galaxy_density_chi_bins_data
    ) noexcept nogil:

    # flat bisp is meant for testing and as such does not check for triangle inequalities!

    cdef double int_width = (chi_max - chi_min) / (num_samples - 1)
    cdef double result = 0
    cdef int i
    for i in range(1, num_samples - 1):
        result += lbs_integrand(chi_min + int_width * i, k1, k2, k3, type1, type2, type3, a_data, b_data, c_data, scale_factor_data, window_data, mps_data, z_at_chi_data)
        # for boundary values use values that lie *just* inside the boundaries to prevent some nasty errors
    result += 0.5 * (
        lbs_integrand(chi_min * 1.01, k1, k2, k3, type1, type2, type3, a_data, b_data, c_data, scale_factor_data, window_data, mps_data, z_at_chi_data) 
        + lbs_integrand(chi_max - 1.0, k1, k2, k3, type1, type2, type3, a_data, b_data, c_data, scale_factor_data, window_data, mps_data, z_at_chi_data)
        )
    result *= int_width

    cdef double fraction_factor = 1 / (1.0 * k1 ** 2 * k2 ** 2 * k3 ** 2) # 1.0 factor to convert to floats
    cdef double const_factor = C**3 * 8

    cdef double pb_correction_val = 0.
    cdef double convergence_to_potential
    if pb_correction:
        convergence_to_potential = 1 / ((k1 * 1.)**2 * (k2 * 1.)**2 * (k3 * 1.)**2 / 8)
        pb_correction_val = convergence_to_potential * lbs_pb_flat(k1, k2, k3, C, mps_data, scale_factor_data, z_at_chi_data, galaxy_density_chi_bins_data, window_data, num_samples, type1, type2, type3)

    return fraction_factor * const_factor * result + pb_correction_val

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

# full sky lensing bispectrum
cpdef double lbs(
    int k1,
    int k2,
    int k3,
    int type1,
    int type2,
    int type3,
    int num_samples,
    bint pb_correction,
    double C,
    double[:, :] a_data,
    double[:, :] b_data,
    double[:, :] c_data,
    double[:] scale_factor_data,
    double[:, :] window_data,
    double[:, :] mps_data,
    double[:] z_at_chi_data,
    double[:, :] galaxy_density_chi_bins_data,
    ) noexcept nogil:
    cdef double wigner_factor
    cdef double sqrt_factor
    cdef double const_factor
    cdef double fraction_factor
    cdef double integration_result

    if (k1 + k2 + k3) % 2 == 0 and k3 <= k1 + k2 and k1 - k2 <= k3 and k2 - k1 <= k3:
        integration_result = lbs_flat(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C, a_data, b_data, c_data, scale_factor_data, window_data, mps_data, z_at_chi_data, galaxy_density_chi_bins_data)

        wigner_factor = fabs(wigner_3j_approx_nocheck(k1, k2, k3))
        sqrt_factor = sqrt((2.0*k1 + 1.0)*(2.0*k2 + 1.0)*(2.0*k3 + 1.0)/(4 * 3.14159))

        return wigner_factor * sqrt_factor * integration_result
    else:
        return 0.0

# fiducial values

cdef double[:] cosm_par_f
cdef double C_f
cdef double[:, :] a_data_f = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_f = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_f = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_f = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_f = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_f = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_f = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_f = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_f = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_f = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_f, C_f, a_data_f, b_data_f, c_data_f, lps_data_f, scale_factor_data_f, window_data_f, mps_data_f, z_at_chi_data_f, cmbps_f, galaxy_density_chi_bins_f = data_import_func('data_fiducial')

cpdef double lbs_flat_f(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs_flat(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_f, a_data_f, b_data_f, c_data_f, scale_factor_data_f, window_data_f, mps_data_f, z_at_chi_data_f, galaxy_density_chi_bins_f)

cdef double lbs_f(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_f, a_data_f, b_data_f, c_data_f, scale_factor_data_f, window_data_f, mps_data_f, z_at_chi_data_f, galaxy_density_chi_bins_f)

cpdef double lbs_pb_lps_f(double chi_source, double chi_source_prime, int type1, int type2, double l, int num_samples) noexcept nogil:
    return lbs_pb_lps(chi_source, chi_source_prime, l, C_f, mps_data_f, scale_factor_data_f, z_at_chi_data_f, galaxy_density_chi_bins_f, num_samples, type1, type2)

cpdef double lbs_pb_flat_f(double l1, double l2, double l3, int type1, int type2, int type3, int num_samples) noexcept nogil:
    return lbs_pb_flat(l1, l2, l3, C_f, mps_data_f, scale_factor_data_f, z_at_chi_data_f, galaxy_density_chi_bins_f, window_data_f, num_samples, type1, type2, type3)
    
cdef double lps_f(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_f, cmbps_f)

cdef bint pink_noise = False
cdef double l_knee = 50
cdef int alpha = 2
print(f'Using pink noise: {pink_noise}')
if pink_noise:
    print(f'l_knee, alpha: {l_knee}, {alpha}')


cdef double cmbps_noise(double l, double sigma, double Delta_X) noexcept nogil:
    # units to input:
    # sigma: arcmin
    # Delta_T, Delta_P: microKelvin arcmin

    # constant time!
    cdef double Tcmb = 2.728e6 # in micro kelvin
    cdef double arcmintorad = 3.14 / 10800
    cdef double sigma_rad = sigma * arcmintorad # in radians
    cdef double Delta_X_rad = Delta_X * arcmintorad
    cdef double white_noise = (Delta_X_rad / Tcmb)**2 * exp(l * (l + 1) * sigma_rad**2 / (8 * log(2)))
    
    if pink_noise:
        return white_noise * (1 + (l_knee / l) ** alpha)
    else:
        return white_noise

cdef double ARCMIN2TOSTERRADIAN = (3.1415 / (180 * 60)) ** 2

cpdef double lps_noise(int l, int type1, int type2) noexcept nogil:
    cdef float noise = 0.

    cdef double sigma
    cdef double Delta_X

    if type1 == 0 and type2 == 0:  # 0 = 'c' (convergence)
        # if cmb_noise_type == 0:
        #     noise = 4. * (l * 1.0)**(-2) * (l + 1.0)**(-2) * conv_noise_data[l-2, 7]
        if cmb_noise_type == 1:
            # stage 3 wide toshiya
            noise = interp.logspace_linear_interp(l, lmin_cmbn, lmax_cmbn, lnum_cmbn, cmbn_106)
        elif cmb_noise_type == 2:
            # stage 4 toshiya
            noise = interp.logspace_linear_interp(l, lmin_cmbn, lmax_cmbn, lnum_cmbn, cmbn_301)
        elif cmb_noise_type == 3:
            noise = interp.logspace_linear_interp(l, lmin_cmbn, lmax_cmbn, lnum_cmbn, planck_noise)

    cdef int num_galaxy_bins = 4

    if (type1 == 1 and type2 == 1) or (type1 == 2 and type2 == 2) or (type1 == 3 and type2 == 3) or (type1 == 4 and type2 == 4):
        if galaxy_noise_type == 1:
            noise = 4. * ((l - 1.) * l * (l + 1.) * (l + 2.))**(-1) * 0.3**2 / 5 * ARCMIN2TOSTERRADIAN * num_galaxy_bins
        if galaxy_noise_type == 2:
            noise = 4. * ((l - 1.) * l * (l + 1.) * (l + 2.))**(-1) * 0.3**2 / 30 * ARCMIN2TOSTERRADIAN * num_galaxy_bins
        if galaxy_noise_type == 3:
            noise = 4. * ((l - 1.) * l * (l + 1.) * (l + 2.))**(-1) * 0.4**2 / 100 * ARCMIN2TOSTERRADIAN * num_galaxy_bins # noise from https://arxiv.org/pdf/astro-ph/0310125

    if type1 == 5 and type2 == 5:  # 5 = 't' (temperature)
        if cmb_noise_type == 1:
            # stage 3 wide toshiya
            sigma = 1
            Delta_X = 6 * 0.71
        elif cmb_noise_type == 2:
            # stage 4 toshiya
            sigma = 3
            Delta_X = 1 # 0.71
        elif cmb_noise_type == 3:
            sigma = 5
            Delta_X = 30
        
        noise = cmbps_noise(l, sigma, Delta_X)

    if type1 == 6 and type2 == 6:  # 6 = 'e' (E-mode)
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
            Delta_X = 52
        
        noise = cmbps_noise(l, sigma, Delta_X)

    if type1 == 7 and type2 == 7:  # 7 = 'b' (B-mode)
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
            Delta_X = 52
        
        noise = cmbps_noise(l, sigma, Delta_X)

    return noise

cdef double lps_f_obs(int l, int type1, int type2) noexcept nogil:

    #########################
    # NOISE: ON
    #########################

    return lps(l, type1, type2, lps_data_f, cmbps_f) + lps_noise(l, type1, type2)

##########################################

cdef double[:] cosm_par_H_p
cdef double C_H_p
cdef double[:, :] a_data_H_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_H_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_H_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_H_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_H_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_H_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_H_p, C_H_p, a_data_H_p, b_data_H_p, c_data_H_p, lps_data_H_p, scale_factor_data_H_p, window_data_H_p, mps_data_H_p, z_at_chi_data_H_p, cmbps_H_p, galaxy_density_chi_bins_H_p = data_import_func('data_H_p')

cdef double lbs_H_p(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_H_p, a_data_H_p, b_data_H_p, c_data_H_p, scale_factor_data_H_p, window_data_H_p, mps_data_H_p, z_at_chi_data_H_p, galaxy_density_chi_bins_H_p)

cdef double lps_H_p(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_H_p, cmbps_H_p)
            

cdef double[:] cosm_par_H_m
cdef double C_H_m
cdef double[:, :] a_data_H_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_H_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_H_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_H_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_H_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_H_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_H_m, C_H_m, a_data_H_m, b_data_H_m, c_data_H_m, lps_data_H_m, scale_factor_data_H_m, window_data_H_m, mps_data_H_m, z_at_chi_data_H_m, cmbps_H_m, galaxy_density_chi_bins_H_m = data_import_func('data_H_m')

cdef double lbs_H_m(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_H_m, a_data_H_m, b_data_H_m, c_data_H_m, scale_factor_data_H_m, window_data_H_m, mps_data_H_m, z_at_chi_data_H_m, galaxy_density_chi_bins_H_m)

cdef double lps_H_m(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_H_m, cmbps_H_m)
            

cdef double[:] cosm_par_ombh2_p
cdef double C_ombh2_p
cdef double[:, :] a_data_ombh2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_ombh2_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_ombh2_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_ombh2_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_ombh2_p, C_ombh2_p, a_data_ombh2_p, b_data_ombh2_p, c_data_ombh2_p, lps_data_ombh2_p, scale_factor_data_ombh2_p, window_data_ombh2_p, mps_data_ombh2_p, z_at_chi_data_ombh2_p, cmbps_ombh2_p, galaxy_density_chi_bins_ombh2_p = data_import_func('data_ombh2_p')

cdef double lbs_ombh2_p(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ombh2_p, a_data_ombh2_p, b_data_ombh2_p, c_data_ombh2_p, scale_factor_data_ombh2_p, window_data_ombh2_p, mps_data_ombh2_p, z_at_chi_data_ombh2_p, galaxy_density_chi_bins_ombh2_p)

cdef double lps_ombh2_p(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_ombh2_p, cmbps_ombh2_p)
            

cdef double[:] cosm_par_ombh2_m
cdef double C_ombh2_m
cdef double[:, :] a_data_ombh2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_ombh2_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_ombh2_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_ombh2_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_ombh2_m, C_ombh2_m, a_data_ombh2_m, b_data_ombh2_m, c_data_ombh2_m, lps_data_ombh2_m, scale_factor_data_ombh2_m, window_data_ombh2_m, mps_data_ombh2_m, z_at_chi_data_ombh2_m, cmbps_ombh2_m, galaxy_density_chi_bins_ombh2_m = data_import_func('data_ombh2_m')

cdef double lbs_ombh2_m(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ombh2_m, a_data_ombh2_m, b_data_ombh2_m, c_data_ombh2_m, scale_factor_data_ombh2_m, window_data_ombh2_m, mps_data_ombh2_m, z_at_chi_data_ombh2_m, galaxy_density_chi_bins_ombh2_m)

cdef double lps_ombh2_m(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_ombh2_m, cmbps_ombh2_m)
            

cdef double[:] cosm_par_omch2_p
cdef double C_omch2_p
cdef double[:, :] a_data_omch2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_omch2_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_omch2_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_omch2_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_omch2_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_omch2_p, C_omch2_p, a_data_omch2_p, b_data_omch2_p, c_data_omch2_p, lps_data_omch2_p, scale_factor_data_omch2_p, window_data_omch2_p, mps_data_omch2_p, z_at_chi_data_omch2_p, cmbps_omch2_p, galaxy_density_chi_bins_omch2_p = data_import_func('data_omch2_p')

cdef double lbs_omch2_p(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_omch2_p, a_data_omch2_p, b_data_omch2_p, c_data_omch2_p, scale_factor_data_omch2_p, window_data_omch2_p, mps_data_omch2_p, z_at_chi_data_omch2_p, galaxy_density_chi_bins_omch2_p)

cdef double lps_omch2_p(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_omch2_p, cmbps_omch2_p)
            

cdef double[:] cosm_par_omch2_m
cdef double C_omch2_m
cdef double[:, :] a_data_omch2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_omch2_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_omch2_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_omch2_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_omch2_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_omch2_m, C_omch2_m, a_data_omch2_m, b_data_omch2_m, c_data_omch2_m, lps_data_omch2_m, scale_factor_data_omch2_m, window_data_omch2_m, mps_data_omch2_m, z_at_chi_data_omch2_m, cmbps_omch2_m, galaxy_density_chi_bins_omch2_m = data_import_func('data_omch2_m')

cdef double lbs_omch2_m(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_omch2_m, a_data_omch2_m, b_data_omch2_m, c_data_omch2_m, scale_factor_data_omch2_m, window_data_omch2_m, mps_data_omch2_m, z_at_chi_data_omch2_m, galaxy_density_chi_bins_omch2_m)

cdef double lps_omch2_m(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_omch2_m, cmbps_omch2_m)
            

cdef double[:] cosm_par_ns_p
cdef double C_ns_p
cdef double[:, :] a_data_ns_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_ns_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_ns_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_ns_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_ns_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_ns_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_ns_p, C_ns_p, a_data_ns_p, b_data_ns_p, c_data_ns_p, lps_data_ns_p, scale_factor_data_ns_p, window_data_ns_p, mps_data_ns_p, z_at_chi_data_ns_p, cmbps_ns_p, galaxy_density_chi_bins_ns_p = data_import_func('data_ns_p')

cdef double lbs_ns_p(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ns_p, a_data_ns_p, b_data_ns_p, c_data_ns_p, scale_factor_data_ns_p, window_data_ns_p, mps_data_ns_p, z_at_chi_data_ns_p, galaxy_density_chi_bins_ns_p)

cdef double lps_ns_p(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_ns_p, cmbps_ns_p)
            

cdef double[:] cosm_par_ns_m
cdef double C_ns_m
cdef double[:, :] a_data_ns_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_ns_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_ns_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_ns_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_ns_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_ns_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_ns_m, C_ns_m, a_data_ns_m, b_data_ns_m, c_data_ns_m, lps_data_ns_m, scale_factor_data_ns_m, window_data_ns_m, mps_data_ns_m, z_at_chi_data_ns_m, cmbps_ns_m, galaxy_density_chi_bins_ns_m = data_import_func('data_ns_m')

cdef double lbs_ns_m(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ns_m, a_data_ns_m, b_data_ns_m, c_data_ns_m, scale_factor_data_ns_m, window_data_ns_m, mps_data_ns_m, z_at_chi_data_ns_m, galaxy_density_chi_bins_ns_m)

cdef double lps_ns_m(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_ns_m, cmbps_ns_m)
            

cdef double[:] cosm_par_As_p
cdef double C_As_p
cdef double[:, :] a_data_As_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_As_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_As_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_As_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_As_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_As_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_As_p, C_As_p, a_data_As_p, b_data_As_p, c_data_As_p, lps_data_As_p, scale_factor_data_As_p, window_data_As_p, mps_data_As_p, z_at_chi_data_As_p, cmbps_As_p, galaxy_density_chi_bins_As_p = data_import_func('data_As_p')

cdef double lbs_As_p(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_As_p, a_data_As_p, b_data_As_p, c_data_As_p, scale_factor_data_As_p, window_data_As_p, mps_data_As_p, z_at_chi_data_As_p, galaxy_density_chi_bins_As_p)

cdef double lps_As_p(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_As_p, cmbps_As_p)
            

cdef double[:] cosm_par_As_m
cdef double C_As_m
cdef double[:, :] a_data_As_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_As_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_As_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_As_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_As_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_As_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_As_m, C_As_m, a_data_As_m, b_data_As_m, c_data_As_m, lps_data_As_m, scale_factor_data_As_m, window_data_As_m, mps_data_As_m, z_at_chi_data_As_m, cmbps_As_m, galaxy_density_chi_bins_As_m = data_import_func('data_As_m')

cdef double lbs_As_m(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_As_m, a_data_As_m, b_data_As_m, c_data_As_m, scale_factor_data_As_m, window_data_As_m, mps_data_As_m, z_at_chi_data_As_m, galaxy_density_chi_bins_As_m)

cdef double lps_As_m(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_As_m, cmbps_As_m)
            

cdef double[:] cosm_par_tau_p
cdef double C_tau_p
cdef double[:, :] a_data_tau_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_tau_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_tau_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_tau_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_tau_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_tau_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_tau_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_tau_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_tau_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_tau_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_tau_p, C_tau_p, a_data_tau_p, b_data_tau_p, c_data_tau_p, lps_data_tau_p, scale_factor_data_tau_p, window_data_tau_p, mps_data_tau_p, z_at_chi_data_tau_p, cmbps_tau_p, galaxy_density_chi_bins_tau_p = data_import_func('data_tau_p')

cdef double lbs_tau_p(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_tau_p, a_data_tau_p, b_data_tau_p, c_data_tau_p, scale_factor_data_tau_p, window_data_tau_p, mps_data_tau_p, z_at_chi_data_tau_p, galaxy_density_chi_bins_tau_p)

cdef double lps_tau_p(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_tau_p, cmbps_tau_p)
            

cdef double[:] cosm_par_tau_m
cdef double C_tau_m
cdef double[:, :] a_data_tau_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_tau_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_tau_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_tau_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_tau_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_tau_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_tau_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_tau_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_tau_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_tau_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_tau_m, C_tau_m, a_data_tau_m, b_data_tau_m, c_data_tau_m, lps_data_tau_m, scale_factor_data_tau_m, window_data_tau_m, mps_data_tau_m, z_at_chi_data_tau_m, cmbps_tau_m, galaxy_density_chi_bins_tau_m = data_import_func('data_tau_m')

cdef double lbs_tau_m(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_tau_m, a_data_tau_m, b_data_tau_m, c_data_tau_m, scale_factor_data_tau_m, window_data_tau_m, mps_data_tau_m, z_at_chi_data_tau_m, galaxy_density_chi_bins_tau_m)

cdef double lps_tau_m(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_tau_m, cmbps_tau_m)
            

cdef double[:] cosm_par_mnu_p
cdef double C_mnu_p
cdef double[:, :] a_data_mnu_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_mnu_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_mnu_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_mnu_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_mnu_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_mnu_p, C_mnu_p, a_data_mnu_p, b_data_mnu_p, c_data_mnu_p, lps_data_mnu_p, scale_factor_data_mnu_p, window_data_mnu_p, mps_data_mnu_p, z_at_chi_data_mnu_p, cmbps_mnu_p, galaxy_density_chi_bins_mnu_p = data_import_func('data_mnu_p')

cdef double lbs_mnu_p(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_mnu_p, a_data_mnu_p, b_data_mnu_p, c_data_mnu_p, scale_factor_data_mnu_p, window_data_mnu_p, mps_data_mnu_p, z_at_chi_data_mnu_p, galaxy_density_chi_bins_mnu_p)

cdef double lps_mnu_p(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_mnu_p, cmbps_mnu_p)
            

cdef double[:] cosm_par_mnu_m
cdef double C_mnu_m
cdef double[:, :] a_data_mnu_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_mnu_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_mnu_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_mnu_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_mnu_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_mnu_m, C_mnu_m, a_data_mnu_m, b_data_mnu_m, c_data_mnu_m, lps_data_mnu_m, scale_factor_data_mnu_m, window_data_mnu_m, mps_data_mnu_m, z_at_chi_data_mnu_m, cmbps_mnu_m, galaxy_density_chi_bins_mnu_m = data_import_func('data_mnu_m')

cdef double lbs_mnu_m(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_mnu_m, a_data_mnu_m, b_data_mnu_m, c_data_mnu_m, scale_factor_data_mnu_m, window_data_mnu_m, mps_data_mnu_m, z_at_chi_data_mnu_m, galaxy_density_chi_bins_mnu_m)

cdef double lps_mnu_m(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_mnu_m, cmbps_mnu_m)
            

cdef double[:] cosm_par_w0_p
cdef double C_w0_p
cdef double[:, :] a_data_w0_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_w0_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_w0_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_w0_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_w0_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_w0_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_w0_p, C_w0_p, a_data_w0_p, b_data_w0_p, c_data_w0_p, lps_data_w0_p, scale_factor_data_w0_p, window_data_w0_p, mps_data_w0_p, z_at_chi_data_w0_p, cmbps_w0_p, galaxy_density_chi_bins_w0_p = data_import_func('data_w0_p')

cdef double lbs_w0_p(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_w0_p, a_data_w0_p, b_data_w0_p, c_data_w0_p, scale_factor_data_w0_p, window_data_w0_p, mps_data_w0_p, z_at_chi_data_w0_p, galaxy_density_chi_bins_w0_p)

cdef double lps_w0_p(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_w0_p, cmbps_w0_p)
            

cdef double[:] cosm_par_w0_m
cdef double C_w0_m
cdef double[:, :] a_data_w0_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_w0_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_w0_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_w0_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_w0_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_w0_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_w0_m, C_w0_m, a_data_w0_m, b_data_w0_m, c_data_w0_m, lps_data_w0_m, scale_factor_data_w0_m, window_data_w0_m, mps_data_w0_m, z_at_chi_data_w0_m, cmbps_w0_m, galaxy_density_chi_bins_w0_m = data_import_func('data_w0_m')

cdef double lbs_w0_m(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_w0_m, a_data_w0_m, b_data_w0_m, c_data_w0_m, scale_factor_data_w0_m, window_data_w0_m, mps_data_w0_m, z_at_chi_data_w0_m, galaxy_density_chi_bins_w0_m)

cdef double lps_w0_m(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_w0_m, cmbps_w0_m)
            

cdef double[:] cosm_par_logT_AGN_p
cdef double C_logT_AGN_p
cdef double[:, :] a_data_logT_AGN_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_logT_AGN_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_logT_AGN_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_logT_AGN_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_logT_AGN_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_logT_AGN_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_logT_AGN_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_logT_AGN_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_logT_AGN_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_logT_AGN_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_logT_AGN_p, C_logT_AGN_p, a_data_logT_AGN_p, b_data_logT_AGN_p, c_data_logT_AGN_p, lps_data_logT_AGN_p, scale_factor_data_logT_AGN_p, window_data_logT_AGN_p, mps_data_logT_AGN_p, z_at_chi_data_logT_AGN_p, cmbps_logT_AGN_p, galaxy_density_chi_bins_logT_AGN_p = data_import_func('data_logT_AGN_p')

cdef double lbs_logT_AGN_p(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_logT_AGN_p, a_data_logT_AGN_p, b_data_logT_AGN_p, c_data_logT_AGN_p, scale_factor_data_logT_AGN_p, window_data_logT_AGN_p, mps_data_logT_AGN_p, z_at_chi_data_logT_AGN_p, galaxy_density_chi_bins_logT_AGN_p)

cdef double lps_logT_AGN_p(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_logT_AGN_p, cmbps_logT_AGN_p)
            

cdef double[:] cosm_par_logT_AGN_m
cdef double C_logT_AGN_m
cdef double[:, :] a_data_logT_AGN_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_logT_AGN_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_logT_AGN_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] lps_data_logT_AGN_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_logT_AGN_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] window_data_logT_AGN_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_logT_AGN_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_logT_AGN_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_logT_AGN_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_logT_AGN_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_logT_AGN_m, C_logT_AGN_m, a_data_logT_AGN_m, b_data_logT_AGN_m, c_data_logT_AGN_m, lps_data_logT_AGN_m, scale_factor_data_logT_AGN_m, window_data_logT_AGN_m, mps_data_logT_AGN_m, z_at_chi_data_logT_AGN_m, cmbps_logT_AGN_m, galaxy_density_chi_bins_logT_AGN_m = data_import_func('data_logT_AGN_m')

cdef double lbs_logT_AGN_m(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_logT_AGN_m, a_data_logT_AGN_m, b_data_logT_AGN_m, c_data_logT_AGN_m, scale_factor_data_logT_AGN_m, window_data_logT_AGN_m, mps_data_logT_AGN_m, z_at_chi_data_logT_AGN_m, galaxy_density_chi_bins_logT_AGN_m)

cdef double lps_logT_AGN_m(int l, int type1, int type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_logT_AGN_m, cmbps_logT_AGN_m)

##########################################


# cdef double der_2o(double f2p, double f1p, double f1m, double f2m, double dx) noexcept nogil:
#     return (-1 * f2p + 8 * f1p - 8 * f1m + f2m) / (12 * dx)

cpdef double der(double fp, double fm, double dx) noexcept nogil:
    return (fp - fm) / (2 * dx)

# cdef double[:] fiducial_cosm_par = np.array([67.4, 0.0224, 0.120, 0.965, 2.1e-9, 0.06])

print((cosm_par_H_p[0] - cosm_par_f[0]) / cosm_par_f[0])
print((cosm_par_ombh2_p[1] - cosm_par_f[1]) / cosm_par_f[1])
print((cosm_par_omch2_p[2] - cosm_par_f[2]) / cosm_par_f[2])
print((cosm_par_ns_p[3] - cosm_par_f[3]) / cosm_par_f[3])
print((cosm_par_As_p[4] - cosm_par_f[4]) / cosm_par_f[4])
print((cosm_par_tau_p[5] - cosm_par_f[5]) / cosm_par_f[5])
print((cosm_par_mnu_p[6] - cosm_par_f[6]) / cosm_par_f[6])
print((cosm_par_w0_p[7] - cosm_par_f[7]) / cosm_par_f[7])
print((cosm_par_logT_AGN_p[8] - cosm_par_f[8]) / cosm_par_f[8])


cdef double lps_der(int k, int type1, int type2, char* par) noexcept nogil:
    # for b'snr' case, just returns the original function
    if par[0] == b's':
        return lps_f(k, type1, type2)

    if par[0] == b'H':
        return der(lps_H_p(k, type1, type2), lps_H_m(k, type1, type2), cosm_par_H_p[0] - cosm_par_f[0]) * cosm_par_f[0]
            
    if par[0] == b'o' and par[2] == b'b':
        return der(lps_ombh2_p(k, type1, type2), lps_ombh2_m(k, type1, type2), cosm_par_ombh2_p[1] - cosm_par_f[1]) * cosm_par_f[1]

    if par[0] == b'o' and par[2] == b'c':
        return der(lps_omch2_p(k, type1, type2), lps_omch2_m(k, type1, type2), cosm_par_omch2_p[2] - cosm_par_f[2]) * cosm_par_f[2]

    if par[0] == b'n':
        return der(lps_ns_p(k, type1, type2), lps_ns_m(k, type1, type2), cosm_par_ns_p[3] - cosm_par_f[3]) * cosm_par_f[3]

    if par[0] == b'A':
        return der(lps_As_p(k, type1, type2), lps_As_m(k, type1, type2), cosm_par_As_p[4] - cosm_par_f[4]) * cosm_par_f[4]
    
    if par[0] == b't':
        return der(lps_tau_p(k, type1, type2), lps_tau_m(k, type1, type2), cosm_par_tau_p[5] - cosm_par_f[5]) * cosm_par_f[5]

    if par[0] == b'm':
        return der(lps_mnu_p(k, type1, type2), lps_mnu_m(k, type1, type2), cosm_par_mnu_p[6] - cosm_par_f[6]) * cosm_par_f[6]

    if par[0] == b'w':
        return der(lps_w0_p(k, type1, type2), lps_w0_m(k, type1, type2), cosm_par_w0_p[7] - cosm_par_f[7]) * cosm_par_f[7]

    if par[0] == b'l':
        return der(lps_logT_AGN_p(k, type1, type2), lps_logT_AGN_m(k, type1, type2), cosm_par_logT_AGN_p[8] - cosm_par_f[8]) * cosm_par_f[8]

cdef double lbs_der(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction, char* par) noexcept nogil:
    # for b'snr' case, just returns the original function
    if par[0] == b's':
        return lbs_f(k1, k2, k3, type1, type2, type3, num_samples, pb_correction)

    if par[0] == b'H':
        return der(lbs_H_p(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), lbs_H_m(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), cosm_par_H_p[0] - cosm_par_f[0]) * cosm_par_f[0]
            
    if par[0] == b'o' and par[2] == b'b':
        return der(lbs_ombh2_p(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), lbs_ombh2_m(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), cosm_par_ombh2_p[1] - cosm_par_f[1]) * cosm_par_f[1]

    if par[0] == b'o' and par[2] == b'c':
        return der(lbs_omch2_p(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), lbs_omch2_m(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), cosm_par_omch2_p[2] - cosm_par_f[2]) * cosm_par_f[2]

    if par[0] == b'n':
        return der(lbs_ns_p(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), lbs_ns_m(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), cosm_par_ns_p[3] - cosm_par_f[3]) * cosm_par_f[3]

    if par[0] == b'A':
        return der(lbs_As_p(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), lbs_As_m(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), cosm_par_As_p[4] - cosm_par_f[4]) * cosm_par_f[4]
    
    if par[0] == b't':
        return der(lbs_tau_p(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), lbs_tau_m(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), cosm_par_tau_p[5] - cosm_par_f[5]) * cosm_par_f[5]

    if par[0] == b'm':
        return der(lbs_mnu_p(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), lbs_mnu_m(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), cosm_par_mnu_p[6] - cosm_par_f[6]) * cosm_par_f[6]

    if par[0] == b'w':
        return der(lbs_w0_p(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), lbs_w0_m(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), cosm_par_w0_p[7] - cosm_par_f[7]) * cosm_par_f[7]

    if par[0] == b'l':
        return der(lbs_logT_AGN_p(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), lbs_logT_AGN_m(k1, k2, k3, type1, type2, type3, num_samples, pb_correction), cosm_par_logT_AGN_p[8] - cosm_par_f[8]) * cosm_par_f[8]
        
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
    return window_func(chi, type, window_data=window_data_f)

def window_func_pb_test(chi, chi_s, source_type):
    return pb_window_func(chi, chi_s, source_type, galaxy_density_chi_bins_f)

def scale_factor_test(chi):
    return scale_factor(chi, scale_factor_data=scale_factor_data_f)

def z_at_chi_test(chi):
    return z_at_chi(chi, z_at_chi_data=z_at_chi_data_f)

def lbs_integrand_test(chi, l1, l2, l3, type1, type2, type3):
    return lbs_integrand(chi, l1, l2, l3, type1, type2, type3, 
    a_data_f, b_data_f, c_data_f, scale_factor_data_f, window_data_f, mps_data_f, z_at_chi_data_f)

def lps_f_obs_test(l, type1, type2):
    return lps_f_obs(l, type1, type2)

def lps_der_test(k, type1, type2, par):
    return lps_der(k, type1, type2, par)

def lbs_der_test(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, par):
    return lbs_der(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, par)

def cmbps_noise_test(l, sigma, Delta_X):
    return cmbps_noise(l, sigma, Delta_X)

def lps_logT_AGN_p_test(l, type1, type2):
    return lps_logT_AGN_p(l, type1, type2)

def lps_logT_AGN_m_test(l, type1, type2):
    return lps_logT_AGN_m(l, type1, type2)


# some code to get shapes of all arrays for debugging purposes

# # Assuming the arrays are already defined
# arrays = {
#     "a_data": a_data_f,
#     "b_data": b_data_f,
#     "c_data": c_data_f,
#     "lps_data": lps_data_f,
#     "scale_factor_data": scale_factor_data_f,
#     "window_data": window_data_f,
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
