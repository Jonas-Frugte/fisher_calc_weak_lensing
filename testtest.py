# follows quadratic estimator of https://arxiv.org/pdf/astro-ph/0111606

print('boba')

import numpy as np
import multiprocessing
import itertools
lmaxintdef = 1

#############################################################################################
configs = (('T', 'T'), ('T', 'E'), ('T', 'B'), ('E', 'E'), ('E', 'B'), ('B', 'B'))
# configs = (('E', 'E'), ('E', 'B'), ('B', 'B'))
#############################################################################################

num_configs = len(configs)
import itertools

def config_mat(L, sigma, Delta_T, Delta_P, num_samples, lmaxint):
    mat = np.zeros((num_samples, num_samples))
    for i, j in itertools.product(range(num_samples), repeat = 2):
        if i >= j:
            result = np.random.rand()
            mat[i, j] = result
            mat[j, i] = result
        # np.array([[N(L, *config1, *config2, sigma, Delta_T, Delta_P, lmaxint = lmaxint, num_samples = num_samples) for config1 in configs] for config2 in configs])
    return mat

def lps_noise(L, sigma, Delta_T, Delta_P, num_samples, lmaxint = lmaxintdef):
    result = 1 / np.sum(np.linalg.inv(config_mat(L, sigma, Delta_T, Delta_P, num_samples, lmaxint)))
    print(result * L**2)
    return result

if __name__ == '__main__':
    print('kiki')
    import numpy as np
    import os
    #import multiprocessing
    #from functools import partial
    # Create a logarithmically spaced array for L.
    lmin = 10
    lmax = 3000
    lnum = 32
    Ls = np.logspace(np.log10(lmin), np.log10(lmax), lnum)
    np.savetxt(f"cmb_noise_files/ls_{lmin}_{lmax}_{lnum}.txt", Ls)
    
    def create_noise_vals(sigma, Delta_T, Delta_P, num_samples):
        # Use multiprocessing to compute N for each L in parallel.
        settings = [[L, sigma, Delta_T, Delta_P, num_samples] for L in Ls]
        pool = multiprocessing.Pool()
        noise_vals = pool.starmap(lps_noise, settings)
        
        # noise_vals = [lps_noise(L, sigma, Delta_T, Delta_P) for L in Ls]
        
        # Ensure the output directory exists.
        out_dir = "cmb_noise_files"
        os.makedirs(out_dir, exist_ok=True)
        
        # Filename includes sigma, Delta_T, and Delta_P.
        filename = f"{out_dir}/Ns_sigma{sigma}_DeltaT{Delta_T}_DeltaP{Delta_P}.txt"
        np.savetxt(filename, noise_vals)
        
        print(f'\nParameters: sigma={sigma}, Delta_T={Delta_T}, Delta_P={Delta_P}\nComputed Ns:\n', noise_vals)
    
    # Noise configurations: (sigma, Delta_T, Delta_P)
    # noise_configs = [
    # #     (1, 0, 6),  # fsky = 0.5
    # #     (1, 0, 3),  # fsky = 0.05
    #     (3, 0, 1),  # fsky = 0.5
    # ]

    noise_configs = [
        (4, 188, 1.41),
    ]

    num_samples = 451
    # Loop over the noise configurations.
    for sigma, Delta_T, Delta_P in noise_configs:
        create_noise_vals(sigma, Delta_T, Delta_P, num_samples)
