import numpy as np
import itertools
from cython.parallel import prange, threadid
from libc.string cimport strncmp

from libc.math cimport isnan, isinf
from libc.stdio cimport printf

from data_importer_new cimport lbs_f, lbs_der, lps_f, lps_f_obs, lps_der
from data_importer_new import get_k_max
from bispec_cov_solver cimport index_2d, index_3d, solve_3
from itertools import *

cdef bint pb_correction = False
cdef int num_samples = 100

def set_pb_correction(val):
    global pb_correction
    pb_correction = val
    print(f'Post-Born corrections are now: {pb_correction}')
    pass


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

cdef double Fisher_mat_single_term(int l1, int l2, int l3, int type, char* par1, char* par2, int num_samples) noexcept nogil:
        return lbs_der(l1, l2, l3, type, type, type, num_samples, pb_correction, par1) * lps_f_obs(l1, type, type)**(-1) \
        * lps_f_obs(l2, type, type)**(-1) * lps_f_obs(l3, type, type)**(-1) * lbs_der(l1, l2, l3, type, type, type, num_samples, pb_correction, par2) / Delta(l1, l2, l3)

cpdef double Fisher_mat_single_alpha_beta(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, int type):
    # print(f'Post-Born corrections: {pb_correction}')

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

cpdef double[:, :] Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, int num_cores, int type, pars_list=None):
    # pars_list: optional list/array of parameter numbers (ints). Default: 1..9
    if pars_list is None:
        pars_py = np.array([1,2,3,4,5,6,7,8,9], dtype=np.int32)
    else:
        pars_py = np.array(pars_list, dtype=np.int32)

    cdef int[:] pars_ints = pars_py
    cdef int n_pars = len(pars_py)
    mat = np.zeros((n_pars, n_pars))

    cdef int i, j
    cdef char[3] name_i, name_j

    for i in range(n_pars):
        for j in range(n_pars):
            if i <= j:
                index_to_par(pars_ints[i], name_i)
                index_to_par(pars_ints[j], name_j)
                tmp = Fisher_mat_single_alpha_beta(lmin, lminbin, lmax, triangle_step_size, num_bispec_samples, &name_i[0], &name_j[0], num_cores, type)
                mat[i, j] = tmp
                mat[j, i] = tmp

    return mat


cdef inline int index_xyzalpha(int X, int Y, int Z, int alpha, int n_tracers) noexcept nogil:
    return index_3d(X, Y, Z, n_tracers) + n_tracers**3 * alpha

cdef inline int fish_index(int alpha, int beta, int n_pars) noexcept nogil:
    return alpha * n_pars + beta

cdef void index_to_par(int par_num, char* par_name) noexcept nogil:
    # par_name should be a char array of atleast length 3 
    if par_num == 0:
        par_name[0] = 115  # 's'
        par_name[1] = 110  # 'n'
        par_name[2] = 114  # 'r'
    elif par_num == 1:
        par_name[0] = 72  # 'H'
        par_name[1] = 32  # ' '
        par_name[2] = 32  # ' '
    elif par_num == 2:
        par_name[0] = 111  # 'o'
        par_name[1] = 109  # 'm'
        par_name[2] = 98   # 'b'
    elif par_num == 3:
        par_name[0] = 111  # 'o'
        par_name[1] = 109  # 'm'
        par_name[2] = 99   # 'c'
    elif par_num == 4:
        par_name[0] = 110  # 'n'
        par_name[1] = 32  # ' '
        par_name[2] = 32  # ' '
    elif par_num == 5:
        par_name[0] = 65  # 'A'
        par_name[1] = 32  # ' '
        par_name[2] = 32  # ' '
    elif par_num == 6:
        par_name[0] = 116  # 't'
        par_name[1] = 32  # ' '
        par_name[2] = 32  # ' '
    elif par_num == 7:
        par_name[0] = 109  # 'm'
        par_name[1] = 32  # ' '
        par_name[2] = 32  # ' '
    elif par_num == 8:
        par_name[0] = 119  # 'w'
        par_name[1] = 32  # ' '
        par_name[2] = 32  # ' '
    elif par_num == 9:
        par_name[0] = 108  # 'l' (logT_AGN)
        par_name[1] = 32  # ' '
        par_name[2] = 32  # ' '

cdef int Fisher_mat_full_term(
    int* tracers,
    int l1, 
    int l2, 
    int l3, 
    double* fish_mat, 
    int num_bispec_samples,
    double* tensor, # n_tracer^3 # TODO: I never use this????
    double* lps_mat_1, # n_tracer^2
    double* lps_mat_2, # n_tracer^2
    double* lps_mat_3, # n_tracer^2
    double* bispec_vecs, # n_tracer^3 x 9 (number of cosm pars
    double* bispec_vecs_solved, # n_tracer^3 x 9 (number of cosm pars, this is C^{-1}C^{-1}C^{-1}B)
    double* bispec_vec_temp, # n_tracer^3,
    double* bispec_vec_solved_temp, # n_tracer^3
    double* vecss_solve3, # n_tracer^3
    double* b_solve3, # n_tracer
    double* x_solve3, # n_tracer
    double* x2_solve3, # n_tracer
    double* vecs_solve2, # n_tracer^2
    double* b_solve2, # n_tracer
    double* x_solve2, # n_tracer
    double* vec_solve1, # n_tracer
    double* A_solve, # n_tracer^2
    double* b_solve, # n_tracer
    double* y_solve, # n_tracer
    int* pars_ints, # n_pars
    int n_tracers,
    int n_pars
    ) noexcept nogil:

    # step 1: calculate powerspectra matrices
    cdef int x, y, z, alpha, beta

    x=0
    while x < n_tracers: # populate powerspectrum matrices
        y=0
        while y < n_tracers:
            
            lps_mat_1[index_2d(x, y, n_tracers)] = lps_f_obs(l1, tracers[x], tracers[y])
            lps_mat_2[index_2d(x, y, n_tracers)] = lps_f_obs(l2, tracers[x], tracers[y])
            lps_mat_3[index_2d(x, y, n_tracers)] = lps_f_obs(l3, tracers[x], tracers[y])
            
            y += 1
        x += 1

    # step 2: calculate bispectra vectors of vectors, this will be 5x5x5x9 with flat index [index_3d(X, Y, Z, n_tracers) + n_tracers**3 * alpha],
    # where alpha is the parameter with respect to which you take the derivative
    
    alpha=0
    cdef char[3] alpha_name # char array of length 3 to hold parameter name for lbs_der
    alpha_name[0] = b' '
    alpha_name[1] = b' '
    alpha_name[2] = b' '
    
    while alpha < n_pars:
        index_to_par(pars_ints[alpha], alpha_name)
        x=0
        while x < n_tracers:
            y=0
            while y < n_tracers:
                z = 0
                while z < n_tracers:
                    # TODO: optimize by only calculating for unique combinations of x, y, z
                    #printf('Calculating bispectrum for tracer combination (%d, %d, %d) and parameter %d\n', x, y, z, alpha)
                    bispec_vecs[index_xyzalpha(x, y, z, alpha, n_tracers)] = lbs_der(l1, l2, l3, tracers[x], tracers[y], tracers[z], num_samples, pb_correction, alpha_name)
                    
                    z+=1
                y+=1
            x+=1
        alpha += 1
    


    # step 3: calculate C^{-1}C^{-1}C^{-1}B once for each alpha

    cdef double tmp

    alpha=0
    while alpha < n_pars:
        x=0
        while x < n_tracers:
            y=0
            while y < n_tracers:
                z = 0
                while z < n_tracers:
                    # TODO: optimize by only calculating for unique combinations of x, y, z
                    bispec_vec_temp[index_3d(x, y, z, n_tracers)] = bispec_vecs[index_xyzalpha(x, y, z, alpha, n_tracers)]
                    
                    z+=1
                y+=1
            x+=1

        solve_3(
            lps_mat_1,
            lps_mat_2,
            lps_mat_3,
            bispec_vec_temp,
            bispec_vec_solved_temp,
            vecss_solve3,
            b_solve3,
            x_solve3,
            x2_solve3,
            vecs_solve2,
            b_solve2,
            x_solve2,
            vec_solve1,
            A_solve,
            b_solve,
            y_solve,
            n_tracers
            )

        x=0
        while x < n_tracers:
            y=0
            while y < n_tracers:
                z = 0
                while z < n_tracers:
                    # TODO: optimize by only calculating for unique combinations of x, y, z

                    bispec_vecs_solved[index_xyzalpha(x, y, z, alpha, n_tracers)] = bispec_vec_solved_temp[index_3d(x, y, z, n_tracers)]
                    
                    z+=1
                y+=1
            x+=1
        alpha += 1
    
    

    # step 4: now contract bispectra[index_3d(X, Y, Z, n_tracers) + n_tracers**3 * alpha] * solved_bispectra[index_3d(X, Y, Z, n_tracers) + n_tracers**3 * beta]
    # and add to [alpha, beta] element of existing Fisher matrix that was given

    alpha = 0
    while alpha < n_pars:
        beta = 0
        while beta < n_pars:
            x = 0
            while x < n_tracers:
                y = 0
                while y < n_tracers:
                    z = 0
                    while z < n_tracers:
                        if alpha < beta:
                            tmp = bispec_vecs[index_xyzalpha(x, y, z, alpha, n_tracers)] * bispec_vecs_solved[index_xyzalpha(x, y, z, beta, n_tracers)]
                            fish_mat[fish_index(alpha, beta, n_pars)] += tmp / Delta(l1, l2, l3)
                            fish_mat[fish_index(beta, alpha, n_pars)] += tmp / Delta(l1, l2, l3)

                        if alpha == beta:
                            tmp = bispec_vecs[index_xyzalpha(x, y, z, alpha, n_tracers)] * bispec_vecs_solved[index_xyzalpha(x, y, z, beta, n_tracers)]
                            fish_mat[fish_index(alpha, beta, n_pars)] += tmp / Delta(l1, l2, l3)

                        z+=1
                    y+=1
                x+=1
            beta+=1
        alpha+=1
    


    return 0

def Fisher_mat_full(
    int[:] tracers, 
    int lmin, 
    int lminbin, 
    int lmax, 
    int triangle_step_size, 
    int num_bispec_samples, 
    int num_cores, 
    int n_tracers,
    pars_list=None
    ):
    # print(f'Post-Born corrections: {pb_correction}')

    if lmax > k_max:
        raise Exception('lmax should be less than upper bound for interpolation (k_max)')

    cdef int[:, :] triangles
    cdef double prop_calc
    triangles, prop_calc = generate_triangles(lmin, lminbin, lmax, triangle_step_size, ordermatters = False)

    cdef int num_samples = int(len(triangles))
    cdef int l1, l2, l3
    cdef int index
    cdef int tid


    # Determine parameter list and counts
    if pars_list is None:
        pars_py = np.array([1,2,3,4,5,6,7,8,9], dtype=np.int32)
    else:
        pars_py = np.array(pars_list, dtype=np.int32)
    cdef int[:] pars_ints = pars_py
    cdef int n_pars = int(len(pars_py))

    # Thread-local storage for each thread's partial Fisher matrix
    cdef int nthreads = num_cores
    cdef double[:] local_fish_mats = np.zeros(nthreads * n_pars * n_pars, dtype = np.float64)

    cdef double[:] tensor = np.zeros(nthreads * n_tracers**3, dtype = np.float64) # n_tracer^3
    cdef double[:] lps_mat_1 = np.zeros(nthreads * n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] lps_mat_2 = np.zeros(nthreads * n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] lps_mat_3 = np.zeros(nthreads * n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] bispec_vecs = np.zeros(nthreads * n_tracers**3 * n_pars, dtype = np.float64) # n_tracer^3 x n_pars
    cdef double[:] bispec_vecs_solved = np.zeros(nthreads * n_tracers**3 * n_pars, dtype = np.float64) # n_tracer^3 x n_pars
    cdef double[:] bispec_vec_temp = np.zeros(nthreads * n_tracers**3, dtype = np.float64) # n_tracer^3,
    cdef double[:] bispec_vec_solved_temp = np.zeros(nthreads * n_tracers**3, dtype = np.float64) # n_tracer^3
    cdef double[:] vecss_solve3 = np.zeros(nthreads * n_tracers**3, dtype = np.float64) # n_tracer^3
    cdef double[:] b_solve3 = np.zeros(nthreads * n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] x_solve3 = np.zeros(nthreads * n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] x2_solve3 = np.zeros(nthreads * n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] vecs_solve2 = np.zeros(nthreads * n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] b_solve2 = np.zeros(nthreads * n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] x_solve2 = np.zeros(nthreads * n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] vec_solve1 = np.zeros(nthreads * n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] A_solve = np.zeros(nthreads * n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] b_solve = np.zeros(nthreads * n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] y_solve = np.zeros(nthreads * n_tracers, dtype = np.float64) # n_tracer
    # Parallel loop: each thread accumulates into its own local array
    for index in prange(0, num_samples, schedule='static', num_threads=num_cores, nogil=True):
        l1 = triangles[index, 0]
        l2 = triangles[index, 1]
        l3 = triangles[index, 2]
        tid = threadid()
        # printf('Thread %d processing triangle (%d, %d, %d)\n', tid, l1, l2, l3)
        # Use the thread's own row as the local array
        Fisher_mat_full_term(
            &tracers[0],
            l1, 
            l2, 
            l3, 
            &local_fish_mats[tid * n_pars * n_pars], 
            num_bispec_samples,
            &tensor[tid * n_tracers**3],
            &lps_mat_1[tid * n_tracers**2],   
            &lps_mat_2[tid * n_tracers**2],  
            &lps_mat_3[tid * n_tracers**2],   
            &bispec_vecs[tid * n_tracers**3 * n_pars],   
            &bispec_vecs_solved[tid * n_tracers**3 * n_pars],   
            &bispec_vec_temp[tid * n_tracers**3],  
            &bispec_vec_solved_temp[tid * n_tracers**3],   
            &vecss_solve3[tid * n_tracers**3],   
            &b_solve3[tid * n_tracers],   
            &x_solve3[tid * n_tracers],   
            &x2_solve3[tid * n_tracers],   
            &vecs_solve2[tid * n_tracers**2],
            &b_solve2[tid * n_tracers],   
            &x_solve2[tid * n_tracers],   
            &vec_solve1[tid * n_tracers],   
            &A_solve[tid * n_tracers**2],  
            &b_solve[tid * n_tracers],   
            &y_solve[tid * n_tracers],
            &pars_ints[0],
            n_tracers,
            n_pars
        )


    # Sum over threads
    cdef double[:] fish_mat_flat = np.zeros(n_pars * n_pars, dtype=np.float64)
    cdef int t, i
    for t in range(nthreads):
        for i in range(n_pars * n_pars):
            fish_mat_flat[i] += local_fish_mats[t * n_pars * n_pars + i]

    fish_mat = np.array([[fish_mat_flat[fish_index(alpha, beta, n_pars)] for alpha in range(n_pars)] for beta in range(n_pars)]) * prop_calc
    return fish_mat

def Fisher_mat(
    tracers: tuple,
    lmin: int,
    lminbin: int,
    lmax: int,
    triangle_step_size: int,
    num_bispec_samples: int,
    num_cores: int,
    pars=None
    ):

    n_tracers = len(tracers)
    if n_tracers == 1:
        mat = Fisher_mat_single(lmin, lminbin, lmax, triangle_step_size, num_bispec_samples, num_cores, tracers[0], pars_list=pars)
    else:
        tracers_int = np.array(tracers, dtype=np.int32) # if you don't do this it will pass array of longs instead of ints and give error
        mat = Fisher_mat_full(tracers_int, lmin, lminbin, lmax, triangle_step_size, num_bispec_samples, num_cores, n_tracers, pars_list=pars)
    
    # check for symmetry of Fisher_matrix
    rtol = 1e-05
    if np.allclose(mat, mat.T, rtol=rtol, atol=1e-30):
        return mat
    else:
        print(f'Fisher matrix not symmetric to withing relative tolerance of {rtol}.')
        #print(f'Relative difference is {}') #TODO: this?
        print('returnng symmetrized version of matrix (F + F^T) / 2')
        return (mat + mat.T) / 2

# Below is some code to test Fisher_mat_full_term by comparing to a simple Python implementation. The test was passed on 26 feb 2026

def Fisher_term_test(l1, l2, l3, tracers, pars=None):
    pars_names = [b's', b'H', b'ombh', b'omch', b'n', b'A', b't', b'm', b'w']
    if pars is None:
        pars = [1,2,3,4,5,6,7,8,9]
    num_samples = 100

    bispec_vecs = np.array(
        [[[[lbs_der(l1, l2, l3, x, y, z, num_samples, pb_correction, pars_names[par]) for x in tracers] for y in tracers] for z in tracers] for par in pars]
        )

    lps_mat_1 = np.array([[lps_f_obs(l1, x, y) for x in tracers] for y in tracers])
    lps_mat_2 = np.array([[lps_f_obs(l2, x, y) for x in tracers] for y in tracers])
    lps_mat_3 = np.array([[lps_f_obs(l3, x, y) for x in tracers] for y in tracers])

    fish_mat_py = np.einsum('ixyz,xa,yb,zc,jabc->ij', bispec_vecs, np.linalg.inv(lps_mat_1), np.linalg.inv(lps_mat_2), np.linalg.inv(lps_mat_3), bispec_vecs) / Delta(l1, l2, l3)

    return np.array(fish_mat_py)

def Fisher_term_pyx_wrap(l1, l2, l3, tracers, pars=None):
    cdef int n_tracers = len(tracers)
    cdef int[:] tracers_c = np.array(tracers).astype(np.int32) # ensure tracers are int32 for Cython
    # default parameter list (same as previous behaviour): 1..9
    if pars is None:
        pars_py = np.array([1,2,3,4,5,6,7,8,9], dtype=np.int32)
    else:
        pars_py = np.array(pars, dtype=np.int32)
    cdef int[:] pars_ints = np.array(pars_py, dtype=np.int32)
    cdef int n_pars = len(pars_py)
    cdef double[:] local_fish_mats = np.zeros(n_pars * n_pars, dtype = np.float64)

    cdef double[:] tensor = np.zeros(n_tracers**3, dtype = np.float64) # n_tracer^3
    cdef double[:] lps_mat_1 = np.zeros(n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] lps_mat_2 = np.zeros(n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] lps_mat_3 = np.zeros(n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] bispec_vecs = np.zeros(n_tracers**3 * n_pars, dtype = np.float64) # n_tracer^3 x n_pars
    cdef double[:] bispec_vecs_solved = np.zeros(n_tracers**3 * n_pars, dtype = np.float64) # n_tracer^3 x n_pars
    cdef double[:] bispec_vec_temp = np.zeros(n_tracers**3, dtype = np.float64) # n_tracer^3,
    cdef double[:] bispec_vec_solved_temp = np.zeros(n_tracers**3, dtype = np.float64) # n_tracer^3
    cdef double[:] vecss_solve3 = np.zeros(n_tracers**3, dtype = np.float64) # n_tracer^3
    cdef double[:] b_solve3 = np.zeros(n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] x_solve3 = np.zeros(n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] x2_solve3 = np.zeros(n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] vecs_solve2 = np.zeros(n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] b_solve2 = np.zeros(n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] x_solve2 = np.zeros(n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] vec_solve1 = np.zeros(n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] A_solve = np.zeros(n_tracers**2, dtype = np.float64) # n_tracer^2
    cdef double[:] b_solve = np.zeros(n_tracers, dtype = np.float64) # n_tracer
    cdef double[:] y_solve = np.zeros(n_tracers, dtype = np.float64) # n_tracer

    Fisher_mat_full_term(
        &tracers_c[0],
        l1, 
        l2, 
        l3, 
        &local_fish_mats[0], 
        100,
        &tensor[0],
        &lps_mat_1[0],   
        &lps_mat_2[0],  
        &lps_mat_3[0],   
        &bispec_vecs[0],   
        &bispec_vecs_solved[0],   
        &bispec_vec_temp[0],  
        &bispec_vec_solved_temp[0],   
        &vecss_solve3[0],   
        &b_solve3[0],   
        &x_solve3[0],   
        &x2_solve3[0],   
        &vecs_solve2[0],
        &b_solve2[0],   
        &x_solve2[0],   
        &vec_solve1[0],   
        &A_solve[0],  
        &b_solve[0],   
        &y_solve[0],
        &pars_ints[0],
        n_tracers,
        n_pars
        )
    
    return np.array([[local_fish_mats[fish_index(alpha, beta, n_pars)] for alpha in range(n_pars)] for beta in range(n_pars)])

def test_fisher_terms(l1, l2, l3, tracers, pars=None):
    '''
    Returns True if the Fisher term calculated by the Cython code matches the Fisher term calculated by the Python code to within a relative tolerance of 1e-5, and False otherwise. This is a test to verify that the Cython code is correctly calculating the Fisher term.
    '''
    return np.allclose(Fisher_term_test(l1, l2, l3, tracers, pars), Fisher_term_pyx_wrap(l1, l2, l3, tracers, pars), rtol=1e-5)

def condition_check(mat, l, tolerance = 1e12):
    if np.linalg.cond(mat) > tolerance:
        print(f'Matrix is poorly conditioned, condition number is {np.linalg.cond(mat)}')
        print(f'l = {l}')
        print(mat)
    pass
    
def Fisher_powersp_el(lmin, lmax, tracers_tuple, par1 = b'snr', par2 = b'snr'):
    result = 0

    xyz_configs = list(combinations_with_replacement(tracers_tuple, 2))
    num_configs = len(xyz_configs)

    for l in range(lmin, lmax + 1 - 1): # at l=2000 sometimes gives poorly conditioned matrix, so skip that one
        vec_l = np.zeros(num_configs)
        vec_r = np.zeros(num_configs)

        for config_num in range(len(xyz_configs)): # TODO: check if the SNR case is also working
            vec_r[config_num] = lps_der(l, xyz_configs[config_num][0], xyz_configs[config_num][1], par2)
            vec_l[config_num] = lps_der(l, xyz_configs[config_num][0], xyz_configs[config_num][1], par1)

        cov_mat = np.zeros((num_configs, num_configs))
        for i, j in product(range(len(xyz_configs)), repeat=2):
            xyz_config_i = xyz_configs[i]
            xyz_config_j = xyz_configs[j]
            for x in permutations((0, 1), 2):
                cov_mat[i, j] += lps_f_obs(l, xyz_config_i[0], xyz_config_j[x[0]]) * lps_f_obs(l, xyz_config_i[1], xyz_config_j[x[1]])

        condition_check(cov_mat, l)

        result += (2 * l + 1) * vec_l.T @ np.linalg.solve(cov_mat, vec_r)
    
    return result

def Fisher_powersp(lmin, lmax, tracer, snr_or_constraints = 'snr'):
    pars = [b'H', b'omb', b'omc', b'n', b'A', b't', b'm', b'w', b'l']


    if type(tracer) == str:
        if tracer == 'c':
            tracers_tuple = (0,)
        elif tracer == 's':
            tracers_tuple = (1, 2, 3, 4)
        elif tracer == 'both':
            tracers_tuple = (0, 1, 2, 3, 4)
        elif tracer == 't':
            tracers_tuple = (5,)
        elif tracer == 'e':
            tracers_tuple = (6,)
        elif tracer == 'bothCMB':
            tracers_tuple = (5, 6)
    else:
        tracers_tuple = tracer

    if snr_or_constraints == 'snr':
        return Fisher_powersp_el(lmin, lmax, tracers_tuple, par1 = b'snr', par2 = b'snr')
    elif snr_or_constraints == 'constraints':
        mat = [[Fisher_powersp_el(lmin, lmax, tracers_tuple, par1 = pars[i], par2 = pars[j]) for i in range(9)] for j in range(9)]
        return mat