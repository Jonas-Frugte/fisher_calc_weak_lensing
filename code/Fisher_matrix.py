import numpy as np
import multiprocessing
import Fisher_calc as vis
from itertools import *

lmin = 2
lmax = 2000
stepsizes = [20, 5 * 20, 10 * 20]
num_bispec_samples = 100 # should be 100 or so
num_cores = 64


def fisher_calc_wrapper(tracers):

    if tracers == 'c':
        return np.array(vis.Fisher_mat((0,), lmin, 0, 200, stepsizes[0], num_bispec_samples, num_cores)) \
            + np.array(vis.Fisher_mat((0,), lmin, 200, 1000, stepsizes[1], num_bispec_samples, num_cores)) \
            + np.array(vis.Fisher_mat((0,), lmin, 1000, lmax, stepsizes[2], num_bispec_samples, num_cores))

    if tracers == 's':
        return np.array(vis.Fisher_mat((1, 2, 3, 4), lmin, 0, 200, stepsizes[0], num_bispec_samples, num_cores)) \
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
            np.savetxt(f'fisher_matrices/fish_mat_bisp_{tracer}_pb.txt', mat)
        else:
            np.savetxt(f'fisher_matrices/fish_mat_bisp_{tracer}.txt', mat)
        print(mat)


if __name__ == '__main__':

    pb_correction = False
    main(pb_correction)

    pb_correction = True
    main(pb_correction)