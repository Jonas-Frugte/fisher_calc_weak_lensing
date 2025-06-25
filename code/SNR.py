from Fisher_calc import Fisher_mat_single, Fisher_mat_full
from Fisher_calc_python_imp import Fisher_powersp_single, Fisher_powersp

import numpy as np

num_cores = 64
stepsize = 49
num_samples = 50

# ls = [0] + ls
ls = [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]
lmin = 2

print('Starting computation')

def ps_SNR2():
    # Fisher_powersp_single(lmin, lmax, tracer, par1 = b'snr', par2 = b'snr')
    print("Convergence powerpectrum SNR^2:")
   
    data = [Fisher_powersp_single(lmin, ls[1], b'c')]
    for i in range(2, len(ls)):
        data.append(data[-1] + Fisher_powersp_single(ls[i - 1], ls[i], b'c'))
    print(data)

    print("Shear powerpectrum SNR^2:")

    data = [Fisher_powersp_single(lmin, ls[1], b's')]
    for i in range(2, len(ls)):
        data.append(data[-1] + Fisher_powersp_single(ls[i - 1], ls[i], b's'))
    print(data)

    # print("Shear + conv powerpectrum SNR^2:")

    # data = [Fisher_powersp(50, ls[1])]
    # for i in range(2, len(ls)):
    #     data.append(data[-1] + Fisher_powersp(ls[i - 1], ls[i]))
    # print(data)

    pass

def single_bs_SNR2():
    # Fisher_mat_single(int lmin, int lminbin, int lmax, int triangle_step_size, int num_bispec_samples, char* par1, char* par2, int num_cores, char* type)
    print("Convergence bispectrum SNR^2:")
   
    data = [Fisher_mat_single(lmin, ls[1-1], ls[1], stepsize, num_samples, b'snr', b'snr', num_cores, b'c')]
    for i in range(2, len(ls)):
        data.append(data[-1] + Fisher_mat_single(lmin, ls[i-1], ls[i], stepsize, num_samples, b'snr', b'snr', num_cores, b'c'))
    print(data)

    print("Shear bispectrum SNR^2:")

    data = [Fisher_mat_single(lmin, ls[1-1], ls[1], stepsize, num_samples, b'snr', b'snr', num_cores, b's')]
    for i in range(2, len(ls)):
        data.append(data[-1] + Fisher_mat_single(lmin, ls[i-1], ls[i], stepsize, num_samples, b'snr', b'snr', num_cores, b's'))
    print(data)

    pass

def all_bs_SNR2():
    # Fisher_mat_full(50, ls[1-1], ls[1], stepsize, 100, b'snr', b'snr', num_cores)
    print("Auto and cross bispectra SNR^2:")

    data = [Fisher_mat_full(lmin, ls[1-1], ls[1], stepsize, 100, b'snr', b'snr', num_cores)]
    for i in range(2, len(ls)):
        data.append(data[-1] + Fisher_mat_full(lmin, ls[1-1], ls[1], stepsize, 100, b'snr', b'snr', num_cores))
    print(data)

    pass

#ps_SNR2()
single_bs_SNR2()
# all_bs_SNR2()
