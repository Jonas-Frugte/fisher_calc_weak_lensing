import Fisher_calc_python_imp as vispy

import numpy as np
from itertools import *

lmin = 2
lmax = 2000

pars = [b'H', b'ombh2', b'omch2', b'ns', b'As', b'tau', b'mnu', b'w0', b'logT_AGN']

def fisher_calc_wrapper_planck(args, tracers):
    print('Planck simulation')
    i, j = args
    # Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores)
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    if tracers == 'c':
        return vispy.Fisher_powersp(8, 400, par1 = pars[i], par2 = pars[j], tracer_name = 'c')
    
def fisher_calc_wrapper_takada_jain(args, tracers):
    i, j = args
    # Fisher_mat_full(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores)
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    if tracers == 's':
        # lmax should actually be 3000, but interpolation currently doesn't go that high
        return vispy.Fisher_powersp(50, 2000, par1 = pars[i], par2 = pars[j], tracer_name = 's')

def main():
    for tracer in ['c', 's', 'both', 't', 'e', 'bothCMB']:
        mat = np.zeros((len(pars), len(pars)))
        for i, j in product(range(len(pars)), repeat = 2):
            result = vispy.Fisher_powersp(lmin, lmax, par1 = pars[i], par2 = pars[j], tracer_name = tracer)
            mat[i, j] = result
            mat[j, i] = result  # Symmetric assignment

        filepath = f'fisher_matrices/fish_mat_powersp_{tracer}.txt'
        np.savetxt(filepath, mat)
        print(f'Created Fisher matrix at {filepath}')
        #print(mat)
    pass

if __name__ == '__main__':
    main()