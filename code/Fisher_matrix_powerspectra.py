import Fisher_calc_python_imp as vispy

import numpy as np
from itertools import *

lmin = 2
lmax = 2000

pars = [b'H', b'ombh2', b'omch2', b'ns', b'mnu', b'tau', b'As', b'w0']
# Wrapper function for multiprocessing
def fisher_calc_wrapper(args, tracers):
    i, j = args
    if tracers == 'c':
        return vispy.Fisher_powersp_single(lmin, lmax, b'c', par1 = pars[i], par2 = pars[j])
    if tracers == 's':
        return vispy.Fisher_powersp_single(lmin, lmax, b's', par1 = pars[i], par2 = pars[j])
    if tracers == 'both':
        return vispy.Fisher_powersp(lmin, lmax, pars[i], pars[j])

def fisher_calc_wrapper_planck(args, tracers):
    print('Planck simulation')
    i, j = args
    # Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores)
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    if tracers == 'c':
        return vispy.Fisher_powersp_single(8, 400, b'c', par1 = pars[i], par2 = pars[j])
    
def fisher_calc_wrapper_takada_jain(args, tracers):
    i, j = args
    # Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores)
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    if tracers == 's':
        # lmax should actually be 3000, but interpolation currently doesn't go that high
        return vispy.Fisher_powersp_single(50, 2000, b's', par1 = pars[i], par2 = pars[j])
    
def fisher_calc_wrapper_cmb(args, tracers):
    print(lmin)
    i, j = args
    # Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores)
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    if tracers == 't':
        return vispy.Fisher_powersp_single(lmin, lmax, b't', par1 = pars[i], par2 = pars[j])
    if tracers == 'e':
        return vispy.Fisher_powersp_single(lmin, lmax, b'e', par1 = pars[i], par2 = pars[j])
    if tracers == 'both':
        # lmax should actually be 3000, but interpolation currently doesn't go that high
        return vispy.Fisher_powersp_cmb(lmin, lmax, par1 = pars[i], par2 = pars[j])

def main():
    for tracer in ['c', 's', 'both']:
        mat = np.zeros((len(pars), len(pars)))
        for i, j in product(range(len(pars)), repeat = 2):
            result = fisher_calc_wrapper((i, j), tracer)
            mat[i, j] = result
            mat[j, i] = result  # Symmetric assignment

        filepath = f'fisher_matrices/fish_mat_powersp_{tracer}.txt'
        np.savetxt(filepath, mat)
        print(f'Created Fisher matrix at {filepath}')
        print(mat)
    pass

if __name__ == '__main__':
    main()