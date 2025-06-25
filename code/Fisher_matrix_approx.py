import numpy as np
import multiprocessing
import Fisher_calc as vis
from itertools import *
lmin = 2
lmax = 2000
stepsizes = [33, 133, 166]

print(f'Stepsizes are: {stepsizes}')

num_bispec_samples = 50
num_cores = 64

pars = [b'H', b'ombh2', b'omch2', b'ns', b'mnu', b'tau', b'As', b'w0']

# Wrapper function for multiprocessing
def fisher_calc_wrapper(args, tracers):
    i, j = args
    # Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores)
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    if tracers == 'c':
        return vis.Fisher_mat_single(lmin, 0, 200, stepsizes[0], num_bispec_samples, pars[i], pars[j], num_cores, b'c') + vis.Fisher_mat_single(lmin, 200, 1000, stepsizes[1], num_bispec_samples, pars[i], pars[j], num_cores, b'c') + vis.Fisher_mat_single(lmin, 1000, lmax, stepsizes[2], num_bispec_samples, pars[i], pars[j], num_cores, b'c')
    if tracers == 's':
        return vis.Fisher_mat_single(lmin, 0, 200, stepsizes[0], num_bispec_samples, pars[i], pars[j], num_cores, b's') + vis.Fisher_mat_single(lmin, 200, 1000, stepsizes[1], num_bispec_samples, pars[i], pars[j], num_cores, b's') + vis.Fisher_mat_single(lmin, 1000, lmax, stepsizes[2], num_bispec_samples, pars[i], pars[j], num_cores, b's')
    if tracers == 'both':
        return vis.Fisher_mat_full(lmin, 0, 200, stepsizes[0], num_bispec_samples, pars[i], pars[j], num_cores) + vis.Fisher_mat_full(lmin, 200, 1000, stepsizes[1], num_bispec_samples, pars[i], pars[j], num_cores) + vis.Fisher_mat_full(lmin, 1000, lmax, stepsizes[2], num_bispec_samples, pars[i], pars[j], num_cores)
    
def fisher_calc_wrapper_planck(args, tracers):
    i, j = args
    # Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores)
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    if tracers == 'c':
        return vis.Fisher_mat_single(8, 0, 400, 8, num_bispec_samples, pars[i], pars[j], num_cores, b'c')

def fisher_calc_wrapper_takada_jain(args, tracers):
    i, j = args
    # Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores)
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    if tracers == 's':
        # lmax should actually be 3000, but interpolation currently doesn't go that high
        return vis.Fisher_mat_single(50, 0, 200, stepsizes[0], num_bispec_samples, pars[i], pars[j], num_cores, b's') + vis.Fisher_mat_single(50, 200, 1000, stepsizes[1], num_bispec_samples, pars[i], pars[j], num_cores, b's') + vis.Fisher_mat_single(50, 1000, 2000, stepsizes[2], num_bispec_samples, pars[i], pars[j], num_cores, b's')
num_pars = len(pars)
print('snoopie2')

def main():
    for tracer in ['both', 's', 'c']:
        mat = np.zeros((num_pars, num_pars))
        counter = 1
        for i in range(num_pars):
            for j in range(i, num_pars):
                result = fisher_calc_wrapper((i, j), tracer)
                mat[i, j] = result
                mat[j, i] = result  # Symmetric assignment
                print(f'{tracer}:', counter, '/', num_pars * (num_pars + 1) / 2)
                counter += 1
        
        np.savetxt(f'fisher_matrices/fish_mat_bisp_{tracer}_pb_exp_fix.txt', mat)
        print(mat)

if __name__ == '__main__':
    main()