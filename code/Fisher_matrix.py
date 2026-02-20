import numpy as np
import multiprocessing
import Fisher_calc as vis
from itertools import *

lmin = 2
lmax = 2000
stepsizes = [20, 5 * 20, 10 * 20]
#stepsizes = [100, 400, 1000]
num_bispec_samples = 100 # should be 100 or so
num_cores = 64

# Wrapper function for multiprocessing
def fisher_calc_wrapper(tracers):
    # Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores)
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    if tracers == 'c':
        return np.array(vis.Fisher_mat_single(lmin, 0, 200, stepsizes[0], num_bispec_samples, num_cores, b'c')) + np.array(vis.Fisher_mat_single(lmin, 200, 1000, stepsizes[1], num_bispec_samples, num_cores, b'c')) + np.array(vis.Fisher_mat_single(lmin, 1000, lmax, stepsizes[2], num_bispec_samples, num_cores, b'c'))
    if tracers == 's':
        return np.array(vis.Fisher_mat_single(lmin, 0, 200, stepsizes[0], num_bispec_samples, num_cores, b's')) + np.array(vis.Fisher_mat_single(lmin, 200, 1000, stepsizes[1], num_bispec_samples, num_cores, b's')) + np.array(vis.Fisher_mat_single(lmin, 1000, lmax, stepsizes[2], num_bispec_samples, num_cores, b's'))
    if tracers == 'both':
        return np.array(vis.Fisher_mat_full(lmin, 0, 200, stepsizes[0], num_bispec_samples, num_cores)) + np.array(vis.Fisher_mat_full(lmin, 200, 1000, stepsizes[1], num_bispec_samples, num_cores)) + np.array(vis.Fisher_mat_full(lmin, 1000, lmax, stepsizes[2], num_bispec_samples, num_cores))

def main(pb_correction):
    for tracer in ['c', 's', 'both']:
        mat = fisher_calc_wrapper(tracer)
        if pb_correction:
            np.savetxt(f'fisher_matrices/fish_mat_bisp_{tracer}_pb.txt', mat)
        else:
            np.savetxt(f'fisher_matrices/fish_mat_bisp_{tracer}.txt', mat)
        print(mat)

if __name__ == '__main__':
    pb_correction = False
    vis.set_pb_correction(pb_correction)
    main(pb_correction)
    pb_correction = True
    vis.set_pb_correction(pb_correction)
    main(pb_correction)


# for python implementation
# if __name__ == "__main__":
#     num_cores = multiprocessing.cpu_count()
    
#     mat = np.zeros((len(pars), len(pars)))
    
#     # Generate index pairs for the lower triangular part of the matrix (excluding diagonal)
#     index_pairs = [(i, j) for i in range(len(pars)) for j in range(i + 1)]
    
#     with multiprocessing.Pool(processes=num_cores) as pool:
#         results = pool.map(fisher_calc_wrapper, index_pairs)
    
#     # Fill the matrix with computed values
#     for i, j, result in results:
#         mat[i, j] = result
#         mat[j, i] = result  # Symmetric assignment

#     print(mat)
    
#     print("Fisher matrix computed successfully!")