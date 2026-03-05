from Fisher_calc import Fisher_mat_single, set_pb_correction
from data_importer_new import set_noise_types

import numpy as np

# Configuration
num_cores = 64
stepsize = 100
num_samples = 50
print(f'Stepsize: {stepsize}')

# l-bins
ls = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]
lmin = 2

# Stage configurations: (cmb_noise_type, galaxy_noise_type)
# CMB types: 1=Stage 3, 2=Stage 4, 3=Planck
# Galaxy types: 1=Stage 3-like, 2=Stage 4-like, 3=optimistic
stages = {
    'Stage 3': (1, 1),
    'Stage 4': (2, 2)
}

print('Starting computation\n')


def compute_snr2_for_stage(stage_name, cmb_type, galaxy_type, post_born=False):
    """Compute cumulative bispectrum SNR^2 for convergence and galaxy bins.

    Uses `pars_list=[0]` to request the `snr` parameter only (returns 1x1 matrix).
    For tracer types: 0 = convergence, 1..4 = galaxy redshift bins.
    """
    set_noise_types(cmb_type, galaxy_type)
    set_pb_correction(post_born)
    
    pb_str = "with Post-Born" if post_born else "without Post-Born"
    print(f"\n{'='*60}")
    print(f"{stage_name} - {pb_str}")
    print(f"CMB noise type: {cmb_type}, Galaxy noise type: {galaxy_type}")
    print(f"{'='*60}\n")
    
    # convergence
    print("Convergence bispectrum SNR^2:")
    val = Fisher_mat_single(lmin, ls[0], ls[1], stepsize, num_samples, num_cores, 0, pars_list=[0])
    # extract scalar from 1x1 matrix
    data = [float(val[0, 0])]
    for i in range(2, len(ls)):
        v = Fisher_mat_single(lmin, ls[i-1], ls[i], stepsize, num_samples, num_cores, 0, pars_list=[0])
        data.append(data[-1] + float(v[0, 0]))
    print(data)

    # galaxy redshift bins (closest and furthest)
    for galaxy_bin in [1, 4]:
        bin_label = "closest" if galaxy_bin == 1 else "furthest"
        print(f"\nGalaxy bin {galaxy_bin} ({bin_label}) bispectrum SNR^2:")
        val = Fisher_mat_single(lmin, ls[0], ls[1], stepsize, num_samples, num_cores, galaxy_bin, pars_list=[0])
        data = [float(val[0, 0])]
        for i in range(2, len(ls)):
            v = Fisher_mat_single(lmin, ls[i-1], ls[i], stepsize, num_samples, num_cores, galaxy_bin, pars_list=[0])
            data.append(data[-1] + float(v[0, 0]))
        print(data)


# Compute for all stages with and without post-Born
for stage_name, (cmb_type, galaxy_type) in stages.items():
    compute_snr2_for_stage(stage_name, cmb_type, galaxy_type, post_born=False)
    compute_snr2_for_stage(stage_name, cmb_type, galaxy_type, post_born=True)

print("\n" + "="*60)
print("Computation complete!")
print("="*60)