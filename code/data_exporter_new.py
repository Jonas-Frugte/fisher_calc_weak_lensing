import pathlib
import numpy as np
import os
import cosm_setup
import multiprocessing




# the smallest k value that the matter bispectrum will be evaluated at will be 1 / chi_last_scatter, this is of the order 10^-4
# if we thus take k_min = 1e-4 then we should have the entire interpolation range evaluated

import interp_settings
exp_par = interp_settings.exp_par

def data_export(folder_name, cosm_par, lps = True, a_create = True, b_create = True, c_create = True, mps = True, rest = True):
    cosm_data = cosm_setup.lensing_spectra(*cosm_par)

    filepath = '/scratch/p319950/data/' + folder_name
    pathlib.Path(filepath).mkdir(exist_ok=True)

    #ks = np.linspace(exp_par['k_min'], exp_par['k_max'], exp_par['k_num'])
    ks = np.logspace(np.log10(exp_par['k_min']), np.log10(exp_par['k_max']), exp_par['k_num'])
    ks_log_fine = np.logspace(np.log10(exp_par['k_min']), np.log10(exp_par['k_max']), exp_par['k_num'] * 25)
    chis = np.linspace(exp_par['chi_min'], exp_par['chi_max'], exp_par['chi_num'])
    zs = np.linspace(exp_par['z_min'], exp_par['z_max'], exp_par['z_num'])
    zs_fine = np.linspace(exp_par['z_min'], exp_par['z_max'], exp_par['z_num'] * 25)
    chis_log = np.logspace(np.log10(exp_par['chi_min']), np.log10(exp_par['chi_max']), exp_par['chi_num'])



    # order of inputs: H0, ombh2, omch2, ns, redshifts, As, mnu,
    cosm_par_dict = {
            'H': cosm_par[0],
            'ombh2': cosm_par[1],
            'omch2': cosm_par[2],
            'ns': cosm_par[3],
            'As': cosm_par[4],
            'mnu': cosm_par[5],
            'w0' : cosm_par[6]
            }

    if rest:
        np.save(filepath + '/cosm_par', cosm_par)
        print(f'cosm_par created for {folder_name}')

        with open(filepath + '/rho_bar', 'w') as file:
            file.write(str(cosm_data.rho_bar))
        print(f'rho_bar file created to {folder_name}')

        # scale factor as func of chi (1d: chi)
        with open(filepath + '/scale_factor', 'w') as file:
            data = [cosm_data.scale_factor(chi) for chi in chis]
            for num in data:
                file.write(str(num) + '\n')
        print(f'scale_factor created to {folder_name}')

        # window func (convergence) (1d: chi)
        np.save(filepath + '/window_func_c', [cosm_data.window_func(chi, 'convergence') for chi in chis_log])
        print(f'window_func_c created to {folder_name}')

        # window func (shear) (1d: chi)
        np.save(filepath + '/window_func_s', [cosm_data.window_func(chi, 'shear') for chi in chis_log])
        print(f'window_func_s created to {folder_name}')

        # redshift at comoving radial distance (1d: chi)
        np.save(filepath + '/z_at_chi', [cosm_data.results.redshift_at_comoving_radial_distance(chi) for chi in chis])
        print(f'z at chi created to {folder_name}')

    if lps:
        # lensing power spectrum (convergence, convergence) (1d: k)
        np.save(filepath + '/lensing_power_spectrum_cc', [cosm_data.lps(k, ('convergence', 'convergence')) for k in ks])
        print(f'lensing power spectrum cc created to {folder_name}')
        
        # lensing power spectrum (convergence, shear) (1d: k)
        np.save(filepath + '/lensing_power_spectrum_cs', [cosm_data.lps(k, ('convergence', 'shear')) for k in ks])
        print(f'lensing power spectrum cs created to {folder_name}')
        
        # lensing power spectrum (shear, shear) (1d: k)
        np.save(filepath + '/lensing_power_spectrum_ss', [cosm_data.lps(k, ('shear', 'shear')) for k in ks])
        print(f'lensing power spectrum ss created to {folder_name}')

    if mps:
        # matter power spectrum (2d: z, k)
        np.save(filepath + '/matter_power_spectrum', [[cosm_data.mps(k, z) for z in zs_fine] for k in ks_log_fine])
        print(f'matter power spectrum created to {folder_name}')

    if a_create:
        print(f'creating a for {folder_name}')
        # a (2d: k, z)
        np.save(filepath + '/a', [[cosm_data.a(k, z) for z in zs_fine] for k in ks_log_fine])
        print(f'a created to {folder_name}')

    if b_create:
        print(f'creating b for {folder_name}')
        # b (2d: k, z)
        np.save(filepath + '/b', [[cosm_data.b(k, z) for z in zs_fine] for k in ks_log_fine])
        print(f'b created to {folder_name}')

    if c_create:
        print(f'creating c for {folder_name}')
        # c (2d: k, z)
        np.save(filepath + '/c', [[cosm_data.c(k, z) for z in zs_fine] for k in ks_log_fine])
        print(f'c created to {folder_name}')

    print(f'Done with exporting/updating data to {folder_name}.')
    pass
   

# first arguments of spectra class: H0=67.4, ombh2=0.0224, omch2=0.120, ns=0.965, As=2.1e-9, mnu=0.06, w0 = -1., wa = 0.

par_names = ('H', 'ombh2', 'omch2', 'ns', 'As', 'mnu', 'w0')
fiducial_cosm_par = np.array([67.4, 0.0223, 0.119, 0.965, 2.13e-9, 0.06, -1])
# based on values that toshiya told me about
cosm_par_delta = np.array([fiducial_cosm_par[0] * 0.1,
                           fiducial_cosm_par[1] * 0.15,
                           fiducial_cosm_par[2] * 0.05,
                           fiducial_cosm_par[3] * 0.01,
                           fiducial_cosm_par[4] * 0.1,
                           fiducial_cosm_par[5] * 0.05, # mnu h should probably be smaller
                           0.06])

num_pars = len(par_names)

# perturb the delta used to calculate derivative in order to test stability
delta_delta_coeffs = [2, 1, 0, -1, -2]

def delta_delta_coeffs_to_str(coeff):
    if coeff == 2:
        return '2p'
    if coeff == 1:
        return '1p'
    if coeff == 0:
        return '0'
    if coeff == -1:
        return '1m'
    if coeff == -2:
        return '2m'
    
delta_delta = 0.05
    
# data_export(folder_name, cosm_par, lps = True, a_create = True, b_create = True, c_create = True, mps = True, rest = True, exp_par = exp_par)

# allows control over which data to create
which_to_create = [[True, True, True, True, True, True]]

exports = [['data_fiducial', fiducial_cosm_par, *create_settings] for create_settings in which_to_create]


for delta_delta_coeff in delta_delta_coeffs:
    for i in range(num_pars):
        for create_settings in which_to_create:
            perturbed_params = [fiducial_cosm_par[k] + cosm_par_delta[k] * (1 + delta_delta_coeff * delta_delta) * int(k == i) for k in range(num_pars)]
            exports.append(
                [f'data_{par_names[i]}_p_{delta_delta_coeffs_to_str(delta_delta_coeff)}', perturbed_params, *create_settings]
            )
            print(perturbed_params)

            perturbed_params = [fiducial_cosm_par[k] - cosm_par_delta[k] * (1 + delta_delta_coeff * delta_delta) * int(k == i) for k in range(num_pars)]
            exports.append(
                [f'data_{par_names[i]}_m_{delta_delta_coeffs_to_str(delta_delta_coeff)}', perturbed_params, *create_settings]
            )


# had to include this because otherwise there was some issue where somewhere the program was already running in parallel and thus the cores were overused when combined with multiprocessing which led to the program spending most of its time just waiting

os.environ["OMP_NUM_THREADS"] = "1"  # OpenMP threads
os.environ["MKL_NUM_THREADS"] = "1"  # MKL threads
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # NumExpr threads
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # OpenBLAS threads

if __name__ == '__main__':
    print(f'Number cpus available: {os.cpu_count()}')
    print(exports)

    num_cores = int(len(exports))
    print(f"Using {num_cores} cores.")
    pool = multiprocessing.Pool(num_cores)
    pool.starmap(data_export, exports)
    # for export in exports:
    #     data_export(*export)