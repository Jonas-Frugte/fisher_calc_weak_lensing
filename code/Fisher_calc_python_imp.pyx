import numpy as np
from itertools import *
from cython.parallel import prange, threadid
from libc.string cimport strncmp
import math

from libc.math cimport isnan, isinf
from libc.stdio cimport printf

from data_importer_new cimport lbs_f, lbs_der, lps_f, lps_f_obs, lps_der

cdef int dd = 0

def distinct_l_inverse_check(l1, l2, l3):

    xyz_configs = list(product((b'c', b's'), repeat = 3))
    
    cov_mat = np.zeros((8, 8))
    cov_mat_inv_th = np.zeros((8, 8))
    cov_mat_inv_nu = np.zeros((8, 8))

    for i, j in product(range(8), repeat = 2):
        xyz_config_i = xyz_configs[i]
        xyz_config_j = xyz_configs[j]

        cov_mat[i, j] += lps_f_obs(l1, xyz_config_i[0], xyz_config_j[0]) * lps_f_obs(l2, xyz_config_i[1], xyz_config_j[1]) * lps_f_obs(l3, xyz_config_i[2], xyz_config_j[2])

        cov_mat_inv_th[i, j] += lensing_power_spectrum_inv(l1, xyz_config_i[0], xyz_config_j[0]) * lensing_power_spectrum_inv(l2, xyz_config_i[1], xyz_config_j[1]) * lensing_power_spectrum_inv(l3, xyz_config_i[2], xyz_config_j[2])
    
    print('Condition numbers:')
    print(np.linalg.cond(cov_mat))
    print(np.linalg.cond(cov_mat_inv_th))

    cov_mat_inv_nu = np.linalg.inv(cov_mat)
    print(np.linalg.cond(cov_mat_inv_nu))

    print('Matrices:')
    print(cov_mat)
    print(cov_mat_inv_th)
    print(cov_mat_inv_nu)

    print('Inverse check:')
    print(cov_mat_inv_th @ cov_mat)
    print(cov_mat_inv_nu @ cov_mat)

    pass

def distinct_l_psd_check(lmin, lmax, num_bispec_samples, par1, par2):

    xyz_configs = list(product((b'c', b's'), repeat = 3))
    
    vec_r = np.zeros(8)
    vec_l = np.zeros(8)
    cov_mat = np.zeros((8, 8))
    xyz_configs = list(product((b'c', b's'), repeat = 3))
    result = 0
    
    for l1 in range(lmin, lmax + 1):
        for l2 in range(l1 + 1, lmax + 1):
            for l3 in range(l2 + 1, lmax + 1):
                for config_num in range(len(xyz_configs)):
                    vec_r[config_num] = lbs_der(l1, l2, l3, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples, par2, dd)
                    vec_l[config_num] = lbs_der(l1, l2, l3, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples, par1, dd)
                    
                for i, j in product(range(8), repeat = 2):
                    xyz_config_i = xyz_configs[i]
                    xyz_config_j = xyz_configs[j]

                    cov_mat[i, j] = lps_f_obs(l1, xyz_config_i[0], xyz_config_j[0]) * lps_f_obs(l2, xyz_config_i[1], xyz_config_j[1]) * lps_f_obs(l3, xyz_config_i[2], xyz_config_j[2])
                
                if min(np.linalg.eigvals(cov_mat)) < 0:
                    print('Not psd.')
                    print(f'l config: {(l1, l2, l3)}')
                    print(f'Matrix: \n {cov_mat}')
                    print(f'Eigenvalues: \n {np.linalg.eigvals(cov_mat)}')
                
                result = vec_l.T @ np.linalg.solve(cov_mat, vec_r)

                if result < 0:
                    print(f'Term result is negative: {result}')
                    print(f'Fisher element: {par1}, {par2}')
                    print(f'l config: {(l1, l2, l3)}')
                    # print(f'Matrix: \n {cov_mat}')
                    print(f'Eigenvalues: \n {np.linalg.eigvals(cov_mat)}')

    pass

def condition_check(mat, tolerance = 1e12):
    if np.linalg.cond(mat) > tolerance:
        print(f'Matrix is poorly conditioned, condition number is {np.linalg.cond(mat)}')
    pass

def Fisher_bisp(lmin, lmax, num_bispec_samples, stepsize, par1 = b'snr', par2 = b'snr'):

    result1 = 0
    result2 = 0
    result3 = 0

    # l1 = l2 = l3 sum
    for l in range(lmin, lmax + 1):
        #print((l, l, l))
        xyz_configs = ((b's', b's', b's'), (b'c', b'c', b'c'), (b's', b's', b'c'), (b'c', b'c', b's'))

        vec_l = np.zeros(4)
        vec_r = np.zeros(4)

        # SNR_case
        if par1 == b'snr' and par2 == b'snr':
            for config_num in range(len(xyz_configs)):
                vec_r[config_num] = lbs_f(l, l, l, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples)
                vec_l = vec_r

        # Fisher matrix case
        else:
            for config_num in range(len(xyz_configs)):
                vec_r[config_num] = lbs_der(l, l, l, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples, par2, dd)
                vec_l[config_num] = lbs_der(l, l, l, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples, par1, dd)

        cov_mat = np.zeros((4, 4))
        for i, j in product(range(4), repeat=2):
            xyz_config_i = xyz_configs[i]
            xyz_config_j = xyz_configs[j]
            for x in permutations((0, 1, 2), 3):
                cov_mat[i, j] += lps_f_obs(l, xyz_config_i[0], xyz_config_j[x[0]]) * lps_f_obs(l, xyz_config_i[1], xyz_config_j[x[1]]) * lps_f_obs(l, xyz_config_i[2], xyz_config_j[x[2]])

        condition_check(cov_mat)

        result1 += vec_l.T @ np.linalg.solve(cov_mat, vec_r)

    # l1 = l2 != l3 sum
    for l1, l2 in product(range(lmin, lmax + 1), repeat = 2):
        if l1 != l2:
            #print((l1, l1, l2))
            xyz_configs = ((b's',b's',b's'), (b's',b'c',b's'), (b'c',b'c',b's'), (b's',b's',b'c'), (b's',b'c',b'c'), (b'c',b'c',b'c'))

            vec = np.zeros(6)
            for config_num in range(len(xyz_configs)):
                vec[config_num] = lbs_f(l1, l1, l2, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples)

            vec_l = np.zeros(6)
            vec_r = np.zeros(6)
            
            # SNR_case
            if par1 == b'snr' and par2 == b'snr':
                for config_num in range(len(xyz_configs)):
                    vec_r[config_num] = lbs_f(l1, l1, l2, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples)
                    vec_l = vec_r

            # Fisher matrix case
            else:
                for config_num in range(len(xyz_configs)):
                    vec_r[config_num] = lbs_der(l1, l1, l2, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples, par2, dd)
                    vec_l[config_num] = lbs_der(l1, l1, l2, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples, par1, dd)

            cov_mat = np.zeros((6, 6))
            for i, j in product(range(6), repeat=2):
                xyz_config_i, xyz_config_j = (xyz_configs[i], xyz_configs[j])
                for x in permutations((0, 1), 2):
                    cov_mat[i, j] += lps_f_obs(l1, xyz_config_i[0], xyz_config_j[x[0]]) * lps_f_obs(l1, xyz_config_i[1], xyz_config_j[x[1]]) * lps_f_obs(l2, xyz_config_i[2], xyz_config_j[2])

            condition_check(cov_mat)

            result2 += vec_l.T @ np.linalg.solve(cov_mat, vec_r)

    # l1 < l2 < l3 sum
    # pure python implementation for testing purposes:
    xyz_configs = list(product((b'c', b's'), repeat = 3))
    
    vec_r = np.zeros(8)
    vec_l = np.zeros(8)
    cov_mat = np.zeros((8, 8))

    result3 = 0
    termresult = 0
    num_triangles_counted = 0
    
    for l1 in range(lmin, lmax + 1, stepsize):
        for l2 in range(l1 + 1, lmax + 1, stepsize):
            for l3 in range(l2 + 1, lmax + 1):
                #print((l1, l2, l3))
                num_triangles_counted += 1

                for config_num in range(len(xyz_configs)):
                    if par1 == b'snr' and par2 == b'snr':
                        vec_r[config_num] = lbs_f_obs(l1, l2, l3, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples)
                        vec_l[config_num] = vec_r

                    else:
                        vec_r[config_num] = lbs_der(l1, l2, l3, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples, par2, dd)
                        vec_l[config_num] = lbs_der(l1, l2, l3, xyz_configs[config_num][0], xyz_configs[config_num][1], xyz_configs[config_num][2], num_bispec_samples, par1, dd)
                    
                for i, j in product(range(8), repeat = 2):
                    xyz_config_i = xyz_configs[i]
                    xyz_config_j = xyz_configs[j]

                    cov_mat[i, j] = lps_f_obs(l1, xyz_config_i[0], xyz_config_j[0]) * lps_f_obs(l2, xyz_config_i[1], xyz_config_j[1]) * lps_f_obs(l3, xyz_config_i[2], xyz_config_j[2])
                
                if min(np.linalg.eigvals(cov_mat)) < 0:
                    print('Not psd.')
                    print(f'l config: {(l1, l2, l3)}')
                    print(f'Matrix: \n {cov_mat}')
                    print(f'Eigenvalues: \n {np.linalg.eigvals(cov_mat)}')
                
                termresult = vec_l.T @ np.linalg.solve(cov_mat, vec_r)

                if termresult < 0 and par1 == par2:
                    print(f'Term result is negative: {termresult}')
                    print(f'l config: {(l1, l2, l3)}')
                    print(f'Matrix: \n {cov_mat}')
                
                result3 += termresult
    if lmax > lmin:
        total_num_triangles = math.comb(lmax - lmin + 1, 3)
        print(total_num_triangles / num_triangles_counted)
        result3 *= total_num_triangles / num_triangles_counted

    print(f'l1 = l2 = l3 case contributes: {result1} \n l1 = l2 != l3 case contributes: {result2} \n l1 < l2 < l3 case contributes: {result3}')

    return sum([result1, result2, result3])

cdef int Delta(int l1, int l2, int l3) noexcept nogil:
    if l1 == l2 == l3:
        return 6
    elif l1 == l2 or l1 == l3 or l2 == l3:
        return 2
    else:
        return 1

def Fisher_bisp_single(lmin, lmax, num_bispec_samples, stepsize, par1 = b'snr', par2 = b'snr', source = b'c'):

    result = 0

    for l1 in range(lmin, lmax + 1, stepsize):
        for l2 in range(l1, lmax + 1, stepsize):
            for l3 in range(l2, lmax + 1):
                result += Delta(l1, l2, l3)**(-1) * lbs_der(l1, l2, l3, source, source, source, num_bispec_samples, par1, dd) * lbs_der(l1, l2, l3, source, source, source, num_bispec_samples, par2, dd) / ( lps_f_obs(l1, source, source) * lps_f_obs(l2, source, source) * lps_f_obs(l3, source, source) )
    
    return result * stepsize**2


def Fisher_powersp(lmin, lmax, par1 = b'snr', par2 = b'snr'):

    result1 = 0

    # l1 = l2 sum
    for l in range(lmin, lmax + 1):
        xyz_configs = ((b's', b's'), (b's', b'c'), (b'c', b'c'))

        vec_l = np.zeros(3)
        vec_r = np.zeros(3)

        # SNR_case
        if par1 == b'snr' and par2 == b'snr':
            for config_num in range(len(xyz_configs)):
                vec_r[config_num] = lps_f(l, xyz_configs[config_num][0], xyz_configs[config_num][1])
                vec_l = vec_r
        # Fisher matrix case
        else:
            for config_num in range(len(xyz_configs)):
                vec_r[config_num] = lps_der(l, xyz_configs[config_num][0], xyz_configs[config_num][1], par2, dd)
                vec_l[config_num] = lps_der(l, xyz_configs[config_num][0], xyz_configs[config_num][1], par1, dd)

        cov_mat = np.zeros((3, 3))
        for i, j in product(range(len(xyz_configs)), repeat=2):
            xyz_config_i = xyz_configs[i]
            xyz_config_j = xyz_configs[j]
            for x in permutations((0, 1), 2):
                cov_mat[i, j] += lps_f_obs(l, xyz_config_i[0], xyz_config_j[x[0]]) * lps_f_obs(l, xyz_config_i[1], xyz_config_j[x[1]])

        condition_check(cov_mat)

        result1 += (2 * l + 1) * vec_l.T @ np.linalg.solve(cov_mat, vec_r) # 2 l + 1 factor is new
    
    return result1

def Fisher_powersp_cmb(lmin, lmax, par1 = b'snr', par2 = b'snr'):

    result1 = 0

    # l1 = l2 sum
    for l in range(lmin, lmax + 1):
        xyz_configs = ((b't', b't'), (b't', b'e'), (b'e', b'e'))

        vec_l = np.zeros(3)
        vec_r = np.zeros(3)

        # SNR_case
        if par1 == b'snr' and par2 == b'snr':
            for config_num in range(len(xyz_configs)):
                vec_r[config_num] = lps_f(l, xyz_configs[config_num][0], xyz_configs[config_num][1])
                vec_l = vec_r
        # Fisher matrix case
        else:
            for config_num in range(len(xyz_configs)):
                vec_r[config_num] = lps_der(l, xyz_configs[config_num][0], xyz_configs[config_num][1], par2, dd)
                vec_l[config_num] = lps_der(l, xyz_configs[config_num][0], xyz_configs[config_num][1], par1, dd)

        cov_mat = np.zeros((3, 3))
        for i, j in product(range(len(xyz_configs)), repeat=2):
            xyz_config_i = xyz_configs[i]
            xyz_config_j = xyz_configs[j]
            for x in permutations((0, 1), 2):
                cov_mat[i, j] += lps_f_obs(l, xyz_config_i[0], xyz_config_j[x[0]]) * lps_f_obs(l, xyz_config_i[1], xyz_config_j[x[1]])

        condition_check(cov_mat)

        result1 += (2 * l + 1) * vec_l.T @ np.linalg.solve(cov_mat, vec_r) # 2 l + 1 factor is new
    
    return result1

def Fisher_powersp_single(lmin, lmax, tracer, par1 = b'snr', par2 = b'snr'):

    result1 = 0

    # l1 = l2 sum
    for l in range(lmin, lmax + 1):

        # SNR_case
        if par1 == b'snr' and par2 == b'snr':
            vec_r = lps_f(l, tracer, tracer)
            vec_l = vec_r
        # Fisher matrix case
        else:
            vec_r = lps_der(l, tracer, tracer, par2, dd)
            vec_l = lps_der(l, tracer, tracer, par1, dd)

        cov_mat = (2 / (2 * l + 1)) * lps_f_obs(l, tracer, tracer)**2 # 2l+1 factor is new

        result1 += vec_l * vec_r / cov_mat
    
    return result1

def Fisher_powersp_single_snr(lmin, lmax, tracer):
    # for testing purposes
    result1 = 0
    for l in range(lmin, lmax + 1):
        result1 += lps_f(l, tracer, tracer)**2 / (2 * lps_f_obs(l, tracer, tracer)**2)

    return result1