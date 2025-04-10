import numpy as np
import json
import traceback
import os

print(f'Number cpus available: {os.cpu_count()}')

# the smallest k value that the matter bispectrum will be evaluated at will be 1 / chi_last_scatter, this is of the order 10^-4
# if we thus take k_min = 1e-6 then we should have the entire interpolation range evaluated

exp_par_fine = {
    'k_min' : 1e-4,
    'k_max' : 2000,
    'k_num' : 2 * 300,
    'chi_min' : 1e-8,
    'chi_max' : 14000, # slightly more than chi_last_scatter
    'chi_num' : 2 * 100 * 50,
    'z_min' : 1e-12,
    'z_max' : 1100, # slightly more than z_last_scatter
    'z_num' : 2 * 100
}

# rougher option
exp_par_rough = {
    'k_min' : 1e-4,
    'k_max' : 2000,
    'k_num' : 300,
    'chi_min' : 1e-8,
    'chi_max' : 14000, # slightly more than chi_last_scatter
    'chi_num' : 100 * 50,
    'z_min' : 1e-12,
    'z_max' : 1100, # slightly more than z_last_scatter
    'z_num' : 100
}

def data_export(folder_name, cosm_par, lps = True, a_create = True, b_create = True, c_create = True, mps = True, rest = True, rough = True):
    import cosm_setup
    cosm_data = cosm_setup.lensing_spectra(*cosm_par)

    # 'rough' version has an interpolation error of around 4 percent and is mainly for testing stuff
    if rough:
        exp_par = exp_par_rough
        filepath = '/scratch/p319950/data_rough/' + folder_name

    # 'fine' version has an error of around 1 percent and is for final calculations
    if not rough:
        exp_par = exp_par_fine
        filepath = '/scratch/p319950/data/' + folder_name

    #ks = np.linspace(exp_par['k_min'], exp_par['k_max'], exp_par['k_num'])
    ks = np.logspace(np.log10(exp_par['k_min']), np.log10(exp_par['k_max']), exp_par['k_num'])
    ks_log_fine = np.logspace(np.log10(exp_par['k_min']), np.log10(exp_par['k_max']), exp_par['k_num'] * 25)
    chis = np.linspace(exp_par['chi_min'], exp_par['chi_max'], exp_par['chi_num'])
    zs = np.linspace(exp_par['z_min'], exp_par['z_max'], exp_par['z_num'])
    zs_fine = np.linspace(exp_par['z_min'], exp_par['z_max'], exp_par['z_num'] * 25)
    chis_log = np.logspace(np.log10(exp_par['chi_min']), np.log10(exp_par['chi_max']), exp_par['chi_num'])



    try:
        # order of inputs: H0, ombh2, omch2, ns, redshifts, As, mnu,
        # MAYBE CHANGE SUCH THAT cosm_par IS INPUTTED AS DICT DIRECTLY INSTEAD OF CONVERTING IT HERE AND THEN EXPORTING
        cosm_par_dict = {
                'H': cosm_par[0],
                'ombh2': cosm_par[1],
                'omch2': cosm_par[2],
                'ns': cosm_par[3],
                'As': cosm_par[4],
                'mnu': cosm_par[5]
                }

        if rest:
            with open(filepath + '/cosm_par.json', 'w') as file:
                json.dump(cosm_par_dict, file)
            print(f'cosm_par_dict created for {folder_name}')

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
            with open(filepath + '/window_func_c', 'w') as file:
                data = [cosm_data.window_func(chi, 'convergence') for chi in chis_log]
                for num in data:
                    file.write(str(num) + '\n')
            print(f'window_func_c created to {folder_name}')

            # window func (shear) (1d: chi)
            with open(filepath + '/window_func_s', 'w') as file:
                data = [cosm_data.window_func(chi, 'shear') for chi in chis_log]
                for num in data:
                    # print('Calculated one window func s value')
                    file.write(str(num) + '\n')
            print(f'window_func_s created to {folder_name}')

            # redshift at comoving radial distance (1d: chi)
            with open(filepath + '/z_at_chi', 'w') as file:
                z_at_chi_data = [cosm_data.results.redshift_at_comoving_radial_distance(chi) for chi in chis]
                for num in z_at_chi_data:
                    file.write(str(num) + '\n')
            print(f'z at chi created to {folder_name}')

        if lps:
            # lensing power spectrum (convergence, convergence) (1d: k)
            with open(filepath + '/lensing_power_spectrum_cc', 'w') as file:
                data = [cosm_data.lps(k, ('convergence', 'convergence')) for k in ks]
                for i in range(len(data)):
                    file.write(str(data[i]) + '\n')
            print(f'lensing power spectrum cc created to {folder_name}')
            
            # lensing power spectrum (convergence, shear) (1d: k)
            with open(filepath + '/lensing_power_spectrum_cs', 'w') as file:
                data = [cosm_data.lps(k, ('convergence', 'shear')) for k in ks]
                for num in data:
                    file.write(str(num) + '\n')
            print(f'lensing power spectrum cs created to {folder_name}')
            
            # lensing power spectrum (shear, shear) (1d: k)
            with open(filepath + '/lensing_power_spectrum_ss', 'w') as file:
                data = [cosm_data.lps(k, ('shear', 'shear')) for k in ks]
                for num in data:
                    file.write(str(num) + '\n')
            print(f'lensing power spectrum ss created to {folder_name}')

        if mps:
            # matter power spectrum (2d: z, k)
            with open(filepath + '/matter_power_spectrum', 'w') as file:
                data = [[cosm_data.mps(k, z) for z in zs_fine] for k in ks_log_fine]
                for row in data:
                    file.write(' '.join(map(str, row)) + '\n')
            print(f'matter power spectrum created to {folder_name}')

        if a_create:
            print(f'creating a for {folder_name}')
            # a (2d: k, z)
            with open(filepath + '/a', 'w') as file:
                data = [[cosm_data.a(k, z) for k in ks_log_fine] for z in zs_fine]
                for row in data:
                    file.write(' '.join(map(str, row)) + '\n')
            print(f'a created to {folder_name}')

        if b_create:
            print(f'creating b for {folder_name}')
            # b (2d: k, z)
            with open(filepath + '/b', 'w') as file:
                data = [[cosm_data.b(k, z) for k in ks_log_fine] for z in zs_fine]
                for row in data:
                    file.write(' '.join(map(str, row)) + '\n')
            print(f'b created to {folder_name}')

        if c_create:
            print(f'creating c for {folder_name}')
            # c (2d: k, z)
            with open(filepath + '/c', 'w') as file:
                data = [[cosm_data.c(k, z) for k in ks_log_fine] for z in zs_fine]
                for row in data:
                    file.write(' '.join(map(str, row)) + '\n')
            print(f'c created to {folder_name}')

        print(f'Done with exporting/updating data to {folder_name}.')
    except Exception as e:
        print(e)
        print(type(e).__name__)
        traceback.print_exc()
    pass
   

# first arguments of spectra class: H0=67.4, ombh2=0.0224, omch2=0.120, ns=0.965, As=2.1e-9, mnu=0.06, w0 = -1., wa = 0.

par_names = ('H', 'ombh2', 'omch2', 'ns', 'As', 'mnu', 'w0')
fiducial_cosm_par = (67.4, 0.0223, 0.119, 0.965, 2.13e-9, 0.06, -1)
num_pars = len(par_names)

dx = 0.025 # if you change this here it also needs to be changed at end of data_importer.pyx !!!

dx_coeffs = [2, 1, -1, -2]

def dx_coeffs_to_str(coeff):
    if coeff == 2:
        return '2p'
    if coeff == 1:
        return '1p'
    if coeff == -1:
        return '1m'
    if coeff == -2:
        return '2m'
    else:
        raise Exception
    
# data_export(folder_name, cosm_par, lps = True, a_create = True, b_create = True, c_create = True, mps = True, rest = True, exp_par = exp_par)

# for 1 + 8 x 4 x 3 = 97 core config
# which_to_create = [[True, False, False, False, True, True], [False, True, False, False, False, False], [False, False, True, True, False, False]]
# which_to_create = [[False, True, False, False, False, False], [False, False, True, False, False, False], [False, False, False, True, False, False]]
which_to_create = [[True, False, False, False, True, True]]
# which_to_create = [[False, False, False, True, False, False]]

# for 1 + 8 x 4 x 1 = 33 core config
# which_to_create = [[False, True, False, False, False, False]]

exports = [['data_fiducial', fiducial_cosm_par, *create_settings] for create_settings in which_to_create]


for dx_coeff in dx_coeffs:
    for i in range(num_pars):
        for create_settings in which_to_create:
            perturbed_params = [fiducial_cosm_par[k] * (1 + dx_coeff * dx * int(k == i)) for k in range(num_pars)]
            
            w_cross_safety = 0.0005
            # here we handle the special case of wa which can't be perturbed relative to its fid value because that's 0
            if i == 7:
                perturbed_params[7] = dx_coeff * dx # * int(7 == i)
                # w is not allowed to cross -1, so to prevent this from happening at a = 1 due to floating point errors we add a small offset to w0
                if dx_coeff == 2 or dx_coeff == 1:
                    perturbed_params[6] = -1 + w_cross_safety
                if dx_coeff == -2 or dx_coeff == -1:
                    perturbed_params[6] = -1 - w_cross_safety


            exports.append(
                ['data_' + par_names[i] + dx_coeffs_to_str(dx_coeff), perturbed_params, *create_settings]
            )

print(exports)

# had to include this because otherwise there was some issue where somewhere the program was already running in parallel and thus the cores were overused when combined with multiprocessing which led to the program spending most of its time just waiting

os.environ["OMP_NUM_THREADS"] = "1"  # OpenMP threads
os.environ["MKL_NUM_THREADS"] = "1"  # MKL threads
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # NumExpr threads
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # OpenBLAS threads

if __name__ == '__main__':
    import multiprocessing

    # num_cores = 25
    num_cores = int(len(exports))
    print(f"Using {num_cores} cores.")
    pool = multiprocessing.Pool(num_cores)
    pool.starmap(data_export, exports)