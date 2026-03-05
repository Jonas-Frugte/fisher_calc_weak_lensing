import faulthandler
faulthandler.enable()

import numpy as np
import Fisher_calc as vis
from itertools import *
from datetime import datetime
from data_importer_new import set_noise_types

np.set_printoptions(precision=1, suppress=False)

lmin = 2
lmax = 2000
stepsizes = [60, 200, 400]
#stepsizes = [20 * 5, 5 * 20 * 5, 10 * 20 * 5] # for testing
num_bispec_samples = 100 # should be 100 or so
num_cores = 128


def fisher_calc_wrapper(tracers):

    if tracers == 'c':
        return np.array(vis.Fisher_mat((0,), lmin, 0, 200, stepsizes[0], num_bispec_samples, num_cores)) \
            + np.array(vis.Fisher_mat((0,), lmin, 200, 1000, stepsizes[1], num_bispec_samples, num_cores)) \
            + np.array(vis.Fisher_mat((0,), lmin, 1000, lmax, stepsizes[2], num_bispec_samples, num_cores))

    if tracers == 's':
        return np.array(vis.Fisher_mat((1, 2, 3, 4), lmin, 0, 100, stepsizes[0], num_bispec_samples, num_cores)) \
            + np.array(vis.Fisher_mat((1, 2, 3, 4), lmin, 200, 1000, stepsizes[1], num_bispec_samples, num_cores)) \
            + np.array(vis.Fisher_mat((1, 2, 3, 4), lmin, 1000, lmax, stepsizes[2], num_bispec_samples, num_cores))

    if tracers == 'both':
        return np.array(vis.Fisher_mat((0, 1, 2, 3, 4), lmin, 0, 200, stepsizes[0], num_bispec_samples, num_cores)) \
            + np.array(vis.Fisher_mat((0, 1, 2, 3, 4), lmin, 200, 1000, stepsizes[1], num_bispec_samples, num_cores)) \
            + np.array(vis.Fisher_mat((0, 1, 2, 3, 4), lmin, 1000, lmax, stepsizes[2], num_bispec_samples, num_cores))


def main(pb_correction):
    vis.set_pb_correction(pb_correction)
    for tracer in ['c', 's', 'both']:
        mat = fisher_calc_wrapper(tracer)
        
        if pb_correction:
            filepath = f'fisher_matrices/fish_mat_bisp_{tracer}_pb.txt'
        else:
            filepath = f'fisher_matrices/fish_mat_bisp_{tracer}.txt'
        
        np.savetxt(filepath, mat)
        print(f"Saved Fisher matrix at: {filepath}, \n at {str(datetime.now())}. \n")


if __name__ == '__main__':
    set_noise_types(2, 2) # Stage 4 noise for both CMB and galaxy

    pb_correction = False
    main(pb_correction)

    pb_correction = True
    main(pb_correction)