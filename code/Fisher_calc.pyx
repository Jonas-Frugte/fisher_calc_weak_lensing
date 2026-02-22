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

cpdef double[:, :] Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, int num_cores, int type):
    pars = [b'H', b'omb', b'omc', b'n', b'A', b't', b'm', b'w', b'l']
    mat = np.zeros((len(pars), len(pars)))

    for i in range(len(pars)):
        for j in range(len(pars)):
            if i <= j:
                tmp = Fisher_mat_single_alpha_beta(lmin, lminbin, lmax, triangle_step_size, num_bispec_samples, pars[i], pars[j], num_cores, type)
                mat[i, j] = tmp
                mat[j, i] = tmp

    return mat


# # CAUTION: function doesn't check for invalid type1, 2 input to optimize performance
# cdef double lensing_power_spectrum_inv(int l, char* type1, char* type2) noexcept nogil:
#     cdef double[:, :] lps_mat = np.zeros((5, 5), dtype=np.float64)
#     cdef int i, j
#     for i in range(5): # TODO: check if I can iterate of i and j like this
#         for j in range(5):
#             lps_mat[i, j] = lps_f_obs(l, tracers[i], tracers[j])
    

#     cdef double det = lps_f_obs(l, b'c', b'c') * lps_f_obs(l, b's', b's') - lps_f_obs(l, b'c', b's')**2
#     if type1[0] == b'c' and type2[0] == b'c':
#         return lps_f_obs(l, b's', b's') / det
#     if type1[0] == b's' and type2[0] == b's':
#         return lps_f_obs(l, b'c', b'c') / det
#     else:
#         return -1 * lps_f_obs(l, b'c', b's') / det

# cdef double Fisher_mat_full_term(int l1, int l2, int l3, char* type11, char* type12, char* type13, char* type21, char* type22, char* type23, char* par1, char* par2, int num_samples) noexcept nogil:
#         return lbs_der(l1, l2, l3, type11, type12, type13, num_samples, pb_correction, par1) * lensing_power_spectrum_inv(l1, type11, type21) \
#         * lensing_power_spectrum_inv(l2, type12, type22) * lensing_power_spectrum_inv(l3, type13, type23) * lbs_der(l1, l2, l3, type21, type22, type23, num_samples, pb_correction, par2) / Delta(l1, l2, l3)

cdef inline int index_xyzalpha(int X, int Y, int Z, int alpha, int n_tracers) noexcept nogil:
    return index_3d(X, Y, Z, n_tracers) + 125 * alpha

cdef inline int fish_index(int alpha, int beta) noexcept nogil:
    return alpha * 9 + beta

cdef void index_to_par(int par_num, char* par_name) noexcept nogil:
    # par_name should be a char array of atleast length 3 

    if par_num == 0:
        par_name[0] = b'H'
    elif par_num == 1:
        par_name[0] = b'o'
        par_name[2] = b'b'
    elif par_num == 2:
        par_name[0] = b'o'
        par_name[2] = b'c'
    elif par_num == 3:
        par_name[0] = b'n'
    elif par_num == 4:
        par_name[0] = b'A'
    elif par_num == 5:
        par_name[0] = b't'
    elif par_num == 6:
        par_name[0] = b'm'
    elif par_num == 7:
        par_name[0] = b'w'
    elif par_num == 8:
        par_name[0] = b'l'  # logT_AGN

cdef int Fisher_mat_full_term_new(
    int* tracers,
    int l1, 
    int l2, 
    int l3, 
    double* fish_mat, 
    int num_bispec_samples,
    double* tensor, # n_tracer^3
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
    int n_tracers
    ) noexcept nogil:
    # this function now returns a matrix 

    # step 1: calculate powerspectra matrices
    cdef int x, y, z, alpha, beta

    cdef char[3] alpha_name
    alpha_name[0] = b' '
    alpha_name[1] = b' '
    alpha_name[2] = b' '

    x=0
    while x < n_tracers: # populate powerspectrum matrices
        y=0
        while y < n_tracers:
            
            lps_mat_1[index_2d(x, y, n_tracers)] = lps_f_obs(l1, tracers[x], tracers[y])
            lps_mat_2[index_2d(x, y, n_tracers)] = lps_f_obs(l2, tracers[x], tracers[y])
            lps_mat_3[index_2d(x, y, n_tracers)] = lps_f_obs(l3, tracers[x], tracers[y])
            
            y += 1
        x += 1

    # step 2: calculate bispectra vectors of vectors, this will be 5x5x5x9 with flat index [index_3d(X, Y, Z, n_tracers) + 125 * alpha],
    # where alpha is the parameter with respect to which you take the derivative
    
    alpha=0
    while alpha < 9: #TODO: check that it's actually 9
        x=0
        while x < n_tracers:
            y=0
            while y < n_tracers:
                z = 0
                while z < n_tracers:
                    
                    index_to_par(alpha, &alpha_name[0])
                    bispec_vecs[index_xyzalpha(x, y, z, alpha, n_tracers)] = lbs_der(l1, l2, l3, tracers[x], tracers[y], tracers[z], num_samples, pb_correction, &alpha_name[0])
                    
                    z+=1
                y+=1
            x+=1
        alpha += 1


    # step 3: calculate C^{-1}C^{-1}C^{-1}B once for each alpha

    cdef double tmp

    alpha=0
    while alpha < 9:
        x=0
        while x < n_tracers:
            y=0
            while y < n_tracers:
                z = 0
                while z < n_tracers:
                    
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
                    
                    bispec_vecs_solved[index_xyzalpha(x, y, z, alpha, n_tracers)] = bispec_vec_solved_temp[index_3d(x, y, z, n_tracers)]
                    
                    z+=1
                y+=1
            x+=1
        alpha += 1

    # step 4: now contract bispectra[index_3d(X, Y, Z, n_tracers) + 125 * alpha] * solved_bispectra[index_3d(X, Y, Z, n_tracers) + 125 * beta]
    # and add to [alpha, beta] element of existing Fisher matrix that was given

    alpha = 0
    while alpha < 9:
        beta = 0
        while beta < 9:
            x = 0
            while x < n_tracers:
                y = 0
                while y < n_tracers:
                    z = 0
                    while z < n_tracers:
                        if alpha < beta:
                            tmp = bispec_vecs[index_xyzalpha(x, y, z, alpha, n_tracers)] * bispec_vecs_solved[index_xyzalpha(x, y, z, beta, n_tracers)]
                            fish_mat[fish_index(alpha, beta)] += tmp
                            fish_mat[fish_index(beta, alpha)] += tmp

                        if alpha == beta:
                            tmp = bispec_vecs[index_xyzalpha(x, y, z, alpha, n_tracers)] * bispec_vecs_solved[index_xyzalpha(x, y, z, beta, n_tracers)]
                            fish_mat[fish_index(alpha, beta)] += tmp

                        z+=1
                    y+=1
                x+=1
            beta+=1
        alpha+=1

    return 0

cpdef double[:, :] Fisher_mat_full(int[:] tracers, int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, int num_cores, int n_tracers):
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


    # Thread-local storage for each thread's partial Fisher matrix
    cdef int nthreads = num_cores
    cdef double[:, :] local_fish_mats = np.zeros((nthreads, 81), dtype=np.float64)

    cdef double[:, :] tensor = np.zeros((nthreads, n_tracers**3)) # n_tracer^3
    cdef double[:, :] lps_mat_1 = np.zeros((nthreads, n_tracers**2)) # n_tracer^2
    cdef double[:, :] lps_mat_2 = np.zeros((nthreads, n_tracers**2)) # n_tracer^2
    cdef double[:, :] lps_mat_3 = np.zeros((nthreads, n_tracers**2)) # n_tracer^2
    cdef double[:, :] bispec_vecs = np.zeros((nthreads, n_tracers**3 * 9)) # n_tracer^3 x 9 (number of cosm pars)
    cdef double[:, :] bispec_vecs_solved = np.zeros((nthreads, n_tracers**3 * 9)) # n_tracer^3 x 9 (number of cosm pars)
    cdef double[:, :] bispec_vec_temp = np.zeros((nthreads, n_tracers**3)) # n_tracer^3,
    cdef double[:, :] bispec_vec_solved_temp = np.zeros((nthreads, n_tracers**3)) # n_tracer^3
    cdef double[:, :] vecss_solve3 = np.zeros((nthreads, n_tracers**3)) # n_tracer^3
    cdef double[:, :] b_solve3 = np.zeros((nthreads, n_tracers)) # n_tracer
    cdef double[:, :] x_solve3 = np.zeros((nthreads, n_tracers)) # n_tracer
    cdef double[:, :] x2_solve3 = np.zeros((nthreads, n_tracers)) # n_tracer
    cdef double[:, :] vecs_solve2 = np.zeros((nthreads, n_tracers**2)) # n_tracer^2
    cdef double[:, :] b_solve2 = np.zeros((nthreads, n_tracers)) # n_tracer
    cdef double[:, :] x_solve2 = np.zeros((nthreads, n_tracers)) # n_tracer
    cdef double[:, :] vec_solve1 = np.zeros((nthreads, n_tracers)) # n_tracer
    cdef double[:, :] A_solve = np.zeros((nthreads, n_tracers**2)) # n_tracer^2
    cdef double[:, :] b_solve = np.zeros((nthreads, n_tracers)) # n_tracer
    cdef double[:, :] y_solve = np.zeros((nthreads, n_tracers)) # n_tracer

    # Parallel loop: each thread accumulates into its own local array
    for index in prange(0, num_samples, schedule='static', num_threads=num_cores, nogil=True):
        l1 = triangles[index, 0]
        l2 = triangles[index, 1]
        l3 = triangles[index, 2]
        tid = threadid()

        # Use the thread's own row as the local array
        Fisher_mat_full_term_new(
            &tracers[0],
            l1, 
            l2, 
            l3, 
            &local_fish_mats[tid, 0], 
            num_bispec_samples,
            &tensor[tid, 0],
            &lps_mat_1[tid, 0],   
            &lps_mat_2[tid, 0],  
            &lps_mat_3[tid, 0],   
            &bispec_vecs[tid, 0],   
            &bispec_vecs_solved[tid, 0],   
            &bispec_vec_temp[tid, 0],  
            &bispec_vec_solved_temp[tid, 0],   
            &vecss_solve3[tid, 0],   
            &b_solve3[tid, 0],   
            &x_solve3[tid, 0],   
            &x2_solve3[tid, 0],   
            &vecs_solve2[tid, 0],
            &b_solve2[tid, 0],   
            &x_solve2[tid, 0],   
            &vec_solve1[tid, 0],   
            &A_solve[tid, 0],  
            &b_solve[tid, 0],   
            &y_solve[tid, 0],
            n_tracers   
        )

    # Sum over threads
    cdef double[81] fish_mat_flat = np.zeros(81)
    cdef int t, i
    for t in range(nthreads):
        for i in range(81):
            fish_mat_flat[i] += local_fish_mats[t, i]

    fish_mat = np.array([[fish_mat_flat[fish_index(alpha, beta)] for alpha in range(9)] for beta in range(9)])

    return fish_mat * prop_calc

def Fisher_mat(
    tracers: tuple,
    lmin: int,
    lminbin: int,
    lmax: int,
    triangle_step_size: int,
    num_bispec_samples: int,
    num_cores: int
  ):

    n_tracers = len(tracers)
    if n_tracers == 1:
        return Fisher_mat_single(lmin, lminbin, lmax, triangle_step_size, num_bispec_samples, num_cores, tracers[0])
    else:
        return Fisher_mat_full(tracers, lmin, lminbin, lmax, triangle_step_size, num_bispec_samples, num_cores, n_tracers)