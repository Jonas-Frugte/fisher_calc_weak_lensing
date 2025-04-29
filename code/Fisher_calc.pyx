import numpy as np
import itertools
from cython.parallel import prange, threadid
from libc.string cimport strncmp

from libc.math cimport isnan, isinf
from libc.stdio cimport printf

from data_importer_new cimport lbs_f, lbs_der, lps_f, lps_f_obs, lps_der
from data_importer_new import get_k_max
from itertools import *

cdef int dd = 0

k_max = get_k_max()

def generate_triangles(lmin, lminbin, lmax, stepsize, ordermatters = False):
    total_num_triangles = 0
    triangles = []
    
    if ordermatters:

        for l1 in range(lmin, lmax + 1):
            for l2 in range(lmin, lmax + 1):
                for l3 in range(lmin, lmax + 1):
                    if max(l1, l2, l3) > lminbin:
                        total_num_triangles += 1

        for l1 in range(lmin, lmax + 1, stepsize):
            for l2 in range(lmin, lmax + 1, stepsize):
                for l3 in range(lmin, lmax + 1):
                    if max(l1, l2, l3) > lminbin:
                        triangles.append([l1, l2, l3])
    
    else:
        for l1 in range(lmin, lmax + 1):
            for l2 in range(l1, lmax + 1):
                for l3 in range(l2, lmax + 1):
                    if max(l1, l2, l3) > lminbin:
                        total_num_triangles += 1

        for l1 in range(lmin, lmax + 1, stepsize):
            for l2 in range(l1, lmax + 1, stepsize):
                for l3 in range(l2, lmax + 1):
                    if max(l1, l2, l3) > lminbin:
                        triangles.append([l1, l2, l3])
        
        # total_num_triangles = math.comb(lmax - lmin + 1 + 2, lmax - lmin + 1 - 1) - math.comb(lminbin - lmin + 1 + 2, lminbin - lmin + 1 - 1)

    prop_calc = total_num_triangles / len(triangles)

    return np.array(triangles, dtype = np.int32), prop_calc

cdef int Delta(int l1, int l2, int l3) noexcept nogil:
    if l1 == l2 == l3:
        return 6
    elif l1 == l2 or l1 == l3 or l2 == l3:
        return 2
    else:
        return 1

cdef double Fisher_mat_single_term(int l1, int l2, int l3, char* type, char* par1, char* par2, int num_samples) noexcept nogil:
        return lbs_der(l1, l2, l3, type, type, type, num_samples, par1, dd) * lps_f_obs(l1, type, type)**(-1) \
        * lps_f_obs(l2, type, type)**(-1) * lps_f_obs(l3, type, type)**(-1) * lbs_der(l1, l2, l3, type, type, type, num_samples, par2, dd) / Delta(l1, l2, l3)

cpdef double Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type):
    if lmax > k_max:
        raise Exception('lmax should be less than upper bound for interpolation (k_max)')

    cdef int[:, :] triangles
    cdef double prop_calc
    triangles, prop_calc = generate_triangles(lmin, lminbin, lmax, triangle_step_size, ordermatters = False)

    cdef int num_samples = int(len(triangles))

    # Create a local result array for each thread
    #cdef double[:] local_results = np.zeros(num_cores, dtype=np.float64)
    cdef double result = 0

    cdef int k1, k2, k3
    cdef int index

    for index in prange(0, num_samples, schedule='static', num_threads=num_cores, nogil=True):
        k1 = triangles[index, 0]
        k2 = triangles[index, 1]
        k3 = triangles[index, 2]
        result += Fisher_mat_single_term(k1, k2, k3, type, par1, par2, num_bispec_samples)

    return result * prop_calc

# CAUTION: function doesn't check for invalid type1, 2 input to optimize performance
cdef double lensing_power_spectrum_inv(int l, char* type1, char* type2) noexcept nogil:
        cdef double det = lps_f_obs(l, b'c', b'c') * lps_f_obs(l, b's', b's') - lps_f_obs(l, b'c', b's')**2
        if type1[0] == b'c' and type2[0] == b'c':
            return lps_f_obs(l, b's', b's') / det
        if type1[0] == b's' and type2[0] == b's':
            return lps_f_obs(l, b'c', b'c') / det
        else:
            return -1 * lps_f_obs(l, b'c', b's') / det

cdef double Fisher_mat_full_term(int l1, int l2, int l3, char* type11, char* type12, char* type13, char* type21, char* type22, char* type23, char* par1, char* par2, int num_samples) noexcept nogil:
        return lbs_der(l1, l2, l3, type11, type12, type13, num_samples, par1, dd) * lensing_power_spectrum_inv(l1, type11, type21) \
        * lensing_power_spectrum_inv(l2, type12, type22) * lensing_power_spectrum_inv(l3, type13, type23) * lbs_der(l1, l2, l3, type21, type22, type23, num_samples, par2, dd) / Delta(l1, l2, l3)

cpdef double Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores):
    if lmax > k_max:
        raise Exception('lmax should be less than upper bound for interpolation (k_max)')

    cdef int[:, :] triangles
    cdef double prop_calc
    triangles, prop_calc = generate_triangles(lmin, lminbin, lmax, triangle_step_size, ordermatters = False)

    cdef int num_samples = int(len(triangles))

    # Create a local result array for each thread
    # cdef double[:] local_results = np.zeros(num_cores, dtype=np.float64)
    cdef double result = 0

    cdef int k1, k2, k3
    cdef char* i
    cdef char* j
    cdef char* k
    cdef char* p
    cdef char* q
    cdef char* r
    cdef int index

    for i, j, k, p, q, r in itertools.product([b'c', b's'], repeat = 6):
        for index in prange(0, num_samples, schedule='static', num_threads=num_cores, nogil=True):
            k1 = triangles[index, 0]
            k2 = triangles[index, 1]
            k3 = triangles[index, 2]
            result += Fisher_mat_full_term(k1, k2, k3, i, j, k, p, q, r, par1, par2, num_bispec_samples)

    return result * prop_calc

########################################
# BELOW WE COMPARE B^2/C^3 and C^2/C^2 #
########################################

def B2P3_dist(lmin, lmax, stepsize, tracer = b's', num_samples = 100, printsamples= False):
    triangles, prop_calc = generate_triangles(lmin, 0, lmax, stepsize, ordermatters = False)

    samples = []
    for triangle in triangles:
        samples.append(Fisher_mat_single_term(triangle[0], triangle[1], triangle[2], tracer, b'snr', b'snr', num_samples))
    if printsamples:
        print(samples)
    return np.mean(samples), np.std(samples)

def P2P2_dist(lmin, lmax, stepsize, tracer = b's', printsamples= False):
    ls = np.arange(lmin, lmax + 1, stepsize)

    samples = []
    for l in ls:
        samples.append(
            0.5 * lps_f(l, tracer, tracer)**2 / lps_f_obs(l, tracer, tracer)**2
        )
    if printsamples:
        print(samples)
    return np.mean(samples), np.std(samples)

