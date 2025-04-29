import Fisher_calc_python_imp as vispy

import numpy as np
from itertools import *

lmin = 2
lmax = 2000

pars = [b'H', b'ombh2', b'omch2', b'ns', b'mnu', b'As', b'w0']

# Wrapper function for multiprocessing
def fisher_calc_wrapper(args, tracers):
    i, j = args
    if tracers == 'c':
        return vispy.Fisher_powersp_single(lmin, lmax, b'c', par1 = pars[i], par2 = pars[j])
    if tracers == 's':
        return vispy.Fisher_powersp_single(lmin, lmax, b's', par1 = pars[i], par2 = pars[j])
    if tracers == 'both':
        return vispy.Fisher_powersp(lmin, lmax, pars[i], pars[j])

for tracer in ['c', 's', 'both']:
    mat = np.zeros((len(pars), len(pars)))
    for i, j in product(range(len(pars)), repeat = 2):
        result = fisher_calc_wrapper((i, j), tracer)
        mat[i, j] = result
        mat[j, i] = result  # Symmetric assignment

    np.savetxt(f'fisher_matrices/fish_mat_powersp_{tracer}.txt', mat)
    print(mat)