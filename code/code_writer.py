
code = ''
for cosm_par in ('H', 'ombh2', 'omch2', 'ns', 'As', 'mnu', 'w0'):
    for pm in ['p', 'm']:
        for delta_delta_coeff in ['2m', '1m', '0', '1p', '2p']:
            code += f'''

cdef double[:] cosm_par_{cosm_par}_{pm}_{delta_delta_coeff}
cdef double C_{cosm_par}_{pm}_{delta_delta_coeff}
cdef double[:, :] a_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_{cosm_par}_{pm}_{delta_delta_coeff} = np.zeros(chi_num, dtype=np.float64)

cosm_par_{cosm_par}_{pm}_{delta_delta_coeff}, C_{cosm_par}_{pm}_{delta_delta_coeff}, a_data_{cosm_par}_{pm}_{delta_delta_coeff}, b_data_{cosm_par}_{pm}_{delta_delta_coeff}, c_data_{cosm_par}_{pm}_{delta_delta_coeff}, lps_cc_data_{cosm_par}_{pm}_{delta_delta_coeff}, lps_cs_data_{cosm_par}_{pm}_{delta_delta_coeff}, lps_ss_data_{cosm_par}_{pm}_{delta_delta_coeff}, scale_factor_data_{cosm_par}_{pm}_{delta_delta_coeff}, window_c_data_{cosm_par}_{pm}_{delta_delta_coeff}, window_s_data_{cosm_par}_{pm}_{delta_delta_coeff}, mps_data_{cosm_par}_{pm}_{delta_delta_coeff}, z_at_chi_data_{cosm_par}_{pm}_{delta_delta_coeff} = data_import_func('data_{cosm_par}_{pm}_{delta_delta_coeff}')

cdef double lbs_{cosm_par}_{pm}_{delta_delta_coeff}(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_{cosm_par}_{pm}_{delta_delta_coeff}, a_data_{cosm_par}_{pm}_{delta_delta_coeff}, b_data_{cosm_par}_{pm}_{delta_delta_coeff}, c_data_{cosm_par}_{pm}_{delta_delta_coeff}, scale_factor_data_{cosm_par}_{pm}_{delta_delta_coeff}, window_c_data_{cosm_par}_{pm}_{delta_delta_coeff}, window_s_data_{cosm_par}_{pm}_{delta_delta_coeff}, mps_data_{cosm_par}_{pm}_{delta_delta_coeff}, z_at_chi_data_{cosm_par}_{pm}_{delta_delta_coeff})

cdef double lps_{cosm_par}_{pm}_{delta_delta_coeff}(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_{cosm_par}_{pm}_{delta_delta_coeff}, lps_cs_data_{cosm_par}_{pm}_{delta_delta_coeff}, lps_ss_data_{cosm_par}_{pm}_{delta_delta_coeff})
            '''

with open('import_code.py', "w") as file:
    print(code, file = file)

