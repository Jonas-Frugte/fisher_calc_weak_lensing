# follows quadratic estimator of https://arxiv.org/pdf/astro-ph/0111606

import cmb_noise_fast as cmbnf
import multiprocessing
import itertools
import numpy as np
import os

#############################################################################################
configs = ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3))#, (3, 3))
# configs = ((2, 2), (2, 3))
# configs = [[2, 3]]
# configs = [[1, 1]]
#############################################################################################

num_configs = len(configs)

def config_mat(L, sigma, Delta_T, Delta_P, num_samples, lmaxint):
    mat = np.zeros((num_configs, num_configs))
    for i, j in itertools.product(range(num_configs), repeat = 2):
        if i >= j:
            result = cmbnf.N(L, *configs[i], *configs[j], sigma, Delta_T, Delta_P, lmaxint, num_samples)
            mat[i, j] = result
            mat[j, i] = result
        # np.array([[N(L, *config1, *config2, sigma, Delta_T, Delta_P, lmaxint = lmaxint, num_samples = num_samples) for config1 in configs] for config2 in configs])
    # print('matrix:', mat, '\n L:', L)
    # print(np.linalg.cond(mat))
    return mat

def lps_noise(L, sigma, Delta_T, Delta_P, num_samples, lmaxint = 10000):
    #result = 1 / np.sum(np.linalg.inv(config_mat(L, sigma, Delta_T, Delta_P, num_samples, lmaxint)))
    result = 1 / sum([1 / cmbnf.N(L, *config, *config, sigma, Delta_T, Delta_P, lmaxint, num_samples) for config in configs])
    print(result * L**2)
    return result

if __name__ == '__main__':
    print('kiki')
    import numpy as np
    import os
    #import multiprocessing
    #from functools import partial
    # Create a logarithmically spaced array for L.
    lmin = 1
    lmax = 3000
    lnum = 64
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
        filename = f"{out_dir}/Ns_sigma{sigma}_DeltaT{Delta_T}_DeltaP{Delta_P}_EB.txt"
        np.savetxt(filename, noise_vals)
        
        print(f'\nParameters: sigma={sigma}, Delta_T={Delta_T}, Delta_P={Delta_P}\nComputed l(l+1)Ns/2pi:\n', np.array(noise_vals) * Ls * (Ls + 1) / 2 / np.pi)

        import matplotlib.pyplot as plt
        plt.loglog(Ls, np.array(noise_vals) * Ls * (Ls + 1) / 2 / np.pi)
        plt.xlim(1, 2000)
        plt.savefig('test.png')
    
    # Noise configurations: (sigma, Delta_T, Delta_P)
    # noise_configs = [
    #     (1, 0, 6),  # fsky = 0.5
    # # #     (1, 0, 3),  # fsky = 0.05
    #     (3, 0, 1),  # fsky = 0.5
    #     #(4, 1, 1.41)
    # ]

    noise_configs = [
        (5, 30, 52),
    ]

    # noise_configs = [
    #     (7, 27, 40 * 1.41),
    # ]

    num_samples = 3001
    # Loop over the noise configurations.
    for sigma, Delta_T, Delta_P in noise_configs:
        create_noise_vals(sigma, Delta_T, Delta_P, num_samples)
