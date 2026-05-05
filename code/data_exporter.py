import pathlib
import numpy as np
import os
import multiprocessing
from itertools import combinations_with_replacement
import cosm_setup




# the smallest k value that the matter bispectrum will be evaluated at will be 1 / chi_last_scatter, this is of the order 10^-4
# if we thus take k_min = 1e-4 then we should have the entire interpolation range evaluated

import interp_settings
exp_par = interp_settings.exp_par

def data_export(folder_name, cosm_par, lps = True, a_create = True, b_create = True, c_create = True, mps = True, rest = True, gal_bins = True):
    cosm_data = cosm_setup.lensing_spectra(*cosm_par)

    filepath = '/scratch/p319950/data_toshiya_like/' + folder_name
    pathlib.Path(filepath).mkdir(exist_ok=True, parents=True)

    ks = np.logspace(np.log10(exp_par['k_min']), np.log10(exp_par['k_max']), exp_par['k_num'])
    ks_log_fine = np.logspace(np.log10(exp_par['k_min']), np.log10(exp_par['k_max']), exp_par['k_num'] * 25)
    chis = np.linspace(exp_par['chi_min'], exp_par['chi_max'], exp_par['chi_num'])
    chis_log_pb = np.logspace(np.log10(exp_par['chi_min']), np.log10(exp_par['chi_max']), exp_par['chi_num_window_func_pb'])
    zs = np.linspace(exp_par['z_min'], exp_par['z_max'], exp_par['z_num'])
    zs_fine = np.linspace(exp_par['z_min'], exp_par['z_max'], exp_par['z_num'] * 25)
    chis_log = np.logspace(np.log10(exp_par['chi_min']), np.log10(exp_par['chi_max']), exp_par['chi_num'])


    # order of inputs: H0, ombh2, omch2, ns, As, tau, mnu, w0, logT_AGN
    cosm_par_dict = {
        'H': cosm_par[0],
        'ombh2': cosm_par[1],
        'omch2': cosm_par[2],
        'ns': cosm_par[3],
        'As': cosm_par[4],
        'tau': cosm_par[5],
        'mnu': cosm_par[6],
        'w0': cosm_par[7],
        'logT_AGN': cosm_par[8]
    }
    
    if rest:
        np.save(filepath + '/cosm_par', cosm_par)
        print(cosm_par)
        print(f'cosm_par created for {folder_name}')

        with open(filepath + '/rho_bar', 'w') as file:
            file.write(str(cosm_data.rho_bar))
        print(f'rho_bar file created to {folder_name}')

        np.save(filepath + '/scale_factor', [cosm_data.scale_factor(chi) for chi in chis])
        print(f'scale_factor created to {folder_name}')

        # window func (convergence) (2d)
        np.save(filepath + '/window_func', [[cosm_data.window_func(chi, tracer) for chi in chis_log] for tracer in ['c', 's1', 's2', 's3', 's4']])
        print(f'window_func created to {folder_name}')

        # redshift at comoving radial distance (1d: chi)
        np.save(filepath + '/z_at_chi', [cosm_data.results.redshift_at_comoving_radial_distance(chi) for chi in chis])
        print(f'z at chi created to {folder_name}')

        lmax = 3000
        cmb_ps = cosm_data.results.get_lensed_scalar_cls(lmax = lmax, raw_cl=True)

        cmb_ps_with_ls = np.zeros((np.shape(cmb_ps)[0], np.shape(cmb_ps)[1] + 1))
        cmb_ps_with_ls[:, 0] = np.arange(lmax + 1)[:]
        cmb_ps_with_ls[:, 1:] = cmb_ps[:, :]

        np.save(filepath + '/cmb_ps_with_ls', cmb_ps_with_ls)
        print(f'cmb_ps_with_ls created to {folder_name}')

    if lps:
        tracers = ('c', 's1', 's2', 's3', 's4')
        tracer_pairs = list(combinations_with_replacement(tracers, 2))
        print(f'Tracer pairs defined in this order: \n {tracer_pairs}')
        lps_data = np.array([[cosm_data.lps(k, tracer_pair) for k in ks_log_fine] for tracer_pair in tracer_pairs])
        np.save(filepath + '/lensing_power_spectrum', lps_data)
        print(f'lensing power spectrum created to {folder_name}, shape: {np.shape(lps_data)}')

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
        np.save(filepath + '/b', [cosm_data.b(k) for k in ks_log_fine])
        print(f'b created to {folder_name}')

    if c_create:
        print(f'creating c for {folder_name}')
        # c (2d: k, z)
        np.save(filepath + '/c', [cosm_data.c(k) for k in ks_log_fine])
        print(f'c created to {folder_name}')

    if gal_bins:
        # galaxy number density for each redshift bin (2d: chi, bin_number)
        np.save(
            filepath + f'/galaxy_density_chi_bins', 
            [[cosm_data.galaxy_density_chi_bin(chi_val, bin_number) for chi_val in chis] for bin_number in range(1, cosm_data.number_bins + 1)]
            )
        print(f'galaxy_density_chi_bins created to {folder_name}')

    print(f'Done with exporting/updating data to {folder_name}.')

par_names = ('H', 'ombh2', 'omch2', 'ns', 'As', 'tau', 'mnu', 'w0', 'logT_AGN')
fiducial_cosm_par = np.array([67.4, 0.0223, 0.119, 0.965, 2.13e-9, 0.063, 0.06, -1, 7.8])

# based on values that toshiya told me about
cosm_par_delta = np.array([
    fiducial_cosm_par[0] * 0.1,
    fiducial_cosm_par[1] * 0.1,
    fiducial_cosm_par[2] * 0.005,
    fiducial_cosm_par[3] * 0.005,
    fiducial_cosm_par[4] * 0.1,
    fiducial_cosm_par[5] * 0.1,
    fiducial_cosm_par[6] * 0.1,
    0.03,
    0.1
])

num_pars = len(par_names)

# Build export list: fiducial + perturbed parameters (±1σ for each parameter)
which_to_create = [[False, True, True, True, False, False, False]]
exports = [['data_fiducial', fiducial_cosm_par, *create_settings] for create_settings in which_to_create]

for i in range(num_pars):
    for create_settings in which_to_create:
        # +1σ perturbation
        perturbed_params = fiducial_cosm_par.copy()
        perturbed_params[i] = fiducial_cosm_par[i] + cosm_par_delta[i]
        exports.append([f'data_{par_names[i]}_p', perturbed_params, *create_settings])

        # -1σ perturbation
        perturbed_params = fiducial_cosm_par.copy()
        perturbed_params[i] = fiducial_cosm_par[i] - cosm_par_delta[i]
        exports.append([f'data_{par_names[i]}_m', perturbed_params, *create_settings])


# Prevent overuse of CPU cores when combined with multiprocessing
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

if __name__ == '__main__':
    print(f'Number cpus available: {os.cpu_count()}')
    print(f'Exports to run: {len(exports)}')

    num_cores = len(exports)
    print(f'Using {num_cores} cores.')
    pool = multiprocessing.Pool(num_cores)
    pool.starmap(data_export, exports)