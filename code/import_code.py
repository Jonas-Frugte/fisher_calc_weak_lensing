

cdef double[:] cosm_par_H_p
cdef double C_H_p
cdef double[:, :] a_data_H_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_H_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_H_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_H_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_H_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_H_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_H_p, C_H_p, a_data_H_p, b_data_H_p, c_data_H_p, lps_data_H_p, scale_factor_data_H_p, window_data_H_p, mps_data_H_p, z_at_chi_data_H_p, cmbps_H_p, galaxy_density_chi_bins_H_p = data_import_func('data_H_p')

cdef double lbs_H_p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_H_p, a_data_H_p, b_data_H_p, c_data_H_p, scale_factor_data_H_p, window_data_H_p, mps_data_H_p, z_at_chi_data_H_p, galaxy_density_chi_bins_H_p)

cdef double lps_H_p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_H_p, cmbps_H_p)
            

cdef double[:] cosm_par_H_m
cdef double C_H_m
cdef double[:, :] a_data_H_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_H_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_H_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_H_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_H_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_H_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_H_m, C_H_m, a_data_H_m, b_data_H_m, c_data_H_m, lps_data_H_m, scale_factor_data_H_m, window_data_H_m, mps_data_H_m, z_at_chi_data_H_m, cmbps_H_m, galaxy_density_chi_bins_H_m = data_import_func('data_H_m')

cdef double lbs_H_m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_H_m, a_data_H_m, b_data_H_m, c_data_H_m, scale_factor_data_H_m, window_data_H_m, mps_data_H_m, z_at_chi_data_H_m, galaxy_density_chi_bins_H_m)

cdef double lps_H_m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_H_m, cmbps_H_m)
            

cdef double[:] cosm_par_ombh2_p
cdef double C_ombh2_p
cdef double[:, :] a_data_ombh2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_ombh2_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_ombh2_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_ombh2_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_ombh2_p, C_ombh2_p, a_data_ombh2_p, b_data_ombh2_p, c_data_ombh2_p, lps_data_ombh2_p, scale_factor_data_ombh2_p, window_data_ombh2_p, mps_data_ombh2_p, z_at_chi_data_ombh2_p, cmbps_ombh2_p, galaxy_density_chi_bins_ombh2_p = data_import_func('data_ombh2_p')

cdef double lbs_ombh2_p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ombh2_p, a_data_ombh2_p, b_data_ombh2_p, c_data_ombh2_p, scale_factor_data_ombh2_p, window_data_ombh2_p, mps_data_ombh2_p, z_at_chi_data_ombh2_p, galaxy_density_chi_bins_ombh2_p)

cdef double lps_ombh2_p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_ombh2_p, cmbps_ombh2_p)
            

cdef double[:] cosm_par_ombh2_m
cdef double C_ombh2_m
cdef double[:, :] a_data_ombh2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_ombh2_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_ombh2_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_ombh2_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_ombh2_m, C_ombh2_m, a_data_ombh2_m, b_data_ombh2_m, c_data_ombh2_m, lps_data_ombh2_m, scale_factor_data_ombh2_m, window_data_ombh2_m, mps_data_ombh2_m, z_at_chi_data_ombh2_m, cmbps_ombh2_m, galaxy_density_chi_bins_ombh2_m = data_import_func('data_ombh2_m')

cdef double lbs_ombh2_m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ombh2_m, a_data_ombh2_m, b_data_ombh2_m, c_data_ombh2_m, scale_factor_data_ombh2_m, window_data_ombh2_m, mps_data_ombh2_m, z_at_chi_data_ombh2_m, galaxy_density_chi_bins_ombh2_m)

cdef double lps_ombh2_m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_ombh2_m, cmbps_ombh2_m)
            

cdef double[:] cosm_par_omch2_p
cdef double C_omch2_p
cdef double[:, :] a_data_omch2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_omch2_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_omch2_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_omch2_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_omch2_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_omch2_p, C_omch2_p, a_data_omch2_p, b_data_omch2_p, c_data_omch2_p, lps_data_omch2_p, scale_factor_data_omch2_p, window_data_omch2_p, mps_data_omch2_p, z_at_chi_data_omch2_p, cmbps_omch2_p, galaxy_density_chi_bins_omch2_p = data_import_func('data_omch2_p')

cdef double lbs_omch2_p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_omch2_p, a_data_omch2_p, b_data_omch2_p, c_data_omch2_p, scale_factor_data_omch2_p, window_data_omch2_p, mps_data_omch2_p, z_at_chi_data_omch2_p, galaxy_density_chi_bins_omch2_p)

cdef double lps_omch2_p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_omch2_p, cmbps_omch2_p)
            

cdef double[:] cosm_par_omch2_m
cdef double C_omch2_m
cdef double[:, :] a_data_omch2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_omch2_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_omch2_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_omch2_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_omch2_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_omch2_m, C_omch2_m, a_data_omch2_m, b_data_omch2_m, c_data_omch2_m, lps_data_omch2_m, scale_factor_data_omch2_m, window_data_omch2_m, mps_data_omch2_m, z_at_chi_data_omch2_m, cmbps_omch2_m, galaxy_density_chi_bins_omch2_m = data_import_func('data_omch2_m')

cdef double lbs_omch2_m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_omch2_m, a_data_omch2_m, b_data_omch2_m, c_data_omch2_m, scale_factor_data_omch2_m, window_data_omch2_m, mps_data_omch2_m, z_at_chi_data_omch2_m, galaxy_density_chi_bins_omch2_m)

cdef double lps_omch2_m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_omch2_m, cmbps_omch2_m)
            

cdef double[:] cosm_par_ns_p
cdef double C_ns_p
cdef double[:, :] a_data_ns_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_ns_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_ns_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_ns_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_ns_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_ns_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_ns_p, C_ns_p, a_data_ns_p, b_data_ns_p, c_data_ns_p, lps_data_ns_p, scale_factor_data_ns_p, window_data_ns_p, mps_data_ns_p, z_at_chi_data_ns_p, cmbps_ns_p, galaxy_density_chi_bins_ns_p = data_import_func('data_ns_p')

cdef double lbs_ns_p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ns_p, a_data_ns_p, b_data_ns_p, c_data_ns_p, scale_factor_data_ns_p, window_data_ns_p, mps_data_ns_p, z_at_chi_data_ns_p, galaxy_density_chi_bins_ns_p)

cdef double lps_ns_p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_ns_p, cmbps_ns_p)
            

cdef double[:] cosm_par_ns_m
cdef double C_ns_m
cdef double[:, :] a_data_ns_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_ns_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_ns_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_ns_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_ns_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_ns_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_ns_m, C_ns_m, a_data_ns_m, b_data_ns_m, c_data_ns_m, lps_data_ns_m, scale_factor_data_ns_m, window_data_ns_m, mps_data_ns_m, z_at_chi_data_ns_m, cmbps_ns_m, galaxy_density_chi_bins_ns_m = data_import_func('data_ns_m')

cdef double lbs_ns_m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ns_m, a_data_ns_m, b_data_ns_m, c_data_ns_m, scale_factor_data_ns_m, window_data_ns_m, mps_data_ns_m, z_at_chi_data_ns_m, galaxy_density_chi_bins_ns_m)

cdef double lps_ns_m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_ns_m, cmbps_ns_m)
            

cdef double[:] cosm_par_As_p
cdef double C_As_p
cdef double[:, :] a_data_As_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_As_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_As_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_As_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_As_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_As_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_As_p, C_As_p, a_data_As_p, b_data_As_p, c_data_As_p, lps_data_As_p, scale_factor_data_As_p, window_data_As_p, mps_data_As_p, z_at_chi_data_As_p, cmbps_As_p, galaxy_density_chi_bins_As_p = data_import_func('data_As_p')

cdef double lbs_As_p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_As_p, a_data_As_p, b_data_As_p, c_data_As_p, scale_factor_data_As_p, window_data_As_p, mps_data_As_p, z_at_chi_data_As_p, galaxy_density_chi_bins_As_p)

cdef double lps_As_p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_As_p, cmbps_As_p)
            

cdef double[:] cosm_par_As_m
cdef double C_As_m
cdef double[:, :] a_data_As_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_As_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_As_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_As_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_As_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_As_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_As_m, C_As_m, a_data_As_m, b_data_As_m, c_data_As_m, lps_data_As_m, scale_factor_data_As_m, window_data_As_m, mps_data_As_m, z_at_chi_data_As_m, cmbps_As_m, galaxy_density_chi_bins_As_m = data_import_func('data_As_m')

cdef double lbs_As_m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_As_m, a_data_As_m, b_data_As_m, c_data_As_m, scale_factor_data_As_m, window_data_As_m, mps_data_As_m, z_at_chi_data_As_m, galaxy_density_chi_bins_As_m)

cdef double lps_As_m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_As_m, cmbps_As_m)
            

cdef double[:] cosm_par_tau_p
cdef double C_tau_p
cdef double[:, :] a_data_tau_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_tau_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_tau_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_tau_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_tau_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_tau_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_tau_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_tau_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_tau_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_tau_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_tau_p, C_tau_p, a_data_tau_p, b_data_tau_p, c_data_tau_p, lps_data_tau_p, scale_factor_data_tau_p, window_data_tau_p, mps_data_tau_p, z_at_chi_data_tau_p, cmbps_tau_p, galaxy_density_chi_bins_tau_p = data_import_func('data_tau_p')

cdef double lbs_tau_p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_tau_p, a_data_tau_p, b_data_tau_p, c_data_tau_p, scale_factor_data_tau_p, window_data_tau_p, mps_data_tau_p, z_at_chi_data_tau_p, galaxy_density_chi_bins_tau_p)

cdef double lps_tau_p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_tau_p, cmbps_tau_p)
            

cdef double[:] cosm_par_tau_m
cdef double C_tau_m
cdef double[:, :] a_data_tau_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_tau_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_tau_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_tau_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_tau_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_tau_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_tau_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_tau_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_tau_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_tau_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_tau_m, C_tau_m, a_data_tau_m, b_data_tau_m, c_data_tau_m, lps_data_tau_m, scale_factor_data_tau_m, window_data_tau_m, mps_data_tau_m, z_at_chi_data_tau_m, cmbps_tau_m, galaxy_density_chi_bins_tau_m = data_import_func('data_tau_m')

cdef double lbs_tau_m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_tau_m, a_data_tau_m, b_data_tau_m, c_data_tau_m, scale_factor_data_tau_m, window_data_tau_m, mps_data_tau_m, z_at_chi_data_tau_m, galaxy_density_chi_bins_tau_m)

cdef double lps_tau_m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_tau_m, cmbps_tau_m)
            

cdef double[:] cosm_par_mnu_p
cdef double C_mnu_p
cdef double[:, :] a_data_mnu_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_mnu_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_mnu_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_mnu_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_mnu_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_mnu_p, C_mnu_p, a_data_mnu_p, b_data_mnu_p, c_data_mnu_p, lps_data_mnu_p, scale_factor_data_mnu_p, window_data_mnu_p, mps_data_mnu_p, z_at_chi_data_mnu_p, cmbps_mnu_p, galaxy_density_chi_bins_mnu_p = data_import_func('data_mnu_p')

cdef double lbs_mnu_p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_mnu_p, a_data_mnu_p, b_data_mnu_p, c_data_mnu_p, scale_factor_data_mnu_p, window_data_mnu_p, mps_data_mnu_p, z_at_chi_data_mnu_p, galaxy_density_chi_bins_mnu_p)

cdef double lps_mnu_p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_mnu_p, cmbps_mnu_p)
            

cdef double[:] cosm_par_mnu_m
cdef double C_mnu_m
cdef double[:, :] a_data_mnu_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_mnu_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_mnu_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_mnu_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_mnu_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_mnu_m, C_mnu_m, a_data_mnu_m, b_data_mnu_m, c_data_mnu_m, lps_data_mnu_m, scale_factor_data_mnu_m, window_data_mnu_m, mps_data_mnu_m, z_at_chi_data_mnu_m, cmbps_mnu_m, galaxy_density_chi_bins_mnu_m = data_import_func('data_mnu_m')

cdef double lbs_mnu_m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_mnu_m, a_data_mnu_m, b_data_mnu_m, c_data_mnu_m, scale_factor_data_mnu_m, window_data_mnu_m, mps_data_mnu_m, z_at_chi_data_mnu_m, galaxy_density_chi_bins_mnu_m)

cdef double lps_mnu_m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_mnu_m, cmbps_mnu_m)
            

cdef double[:] cosm_par_w0_p
cdef double C_w0_p
cdef double[:, :] a_data_w0_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_w0_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_w0_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_w0_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_w0_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_w0_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_w0_p, C_w0_p, a_data_w0_p, b_data_w0_p, c_data_w0_p, lps_data_w0_p, scale_factor_data_w0_p, window_data_w0_p, mps_data_w0_p, z_at_chi_data_w0_p, cmbps_w0_p, galaxy_density_chi_bins_w0_p = data_import_func('data_w0_p')

cdef double lbs_w0_p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_w0_p, a_data_w0_p, b_data_w0_p, c_data_w0_p, scale_factor_data_w0_p, window_data_w0_p, mps_data_w0_p, z_at_chi_data_w0_p, galaxy_density_chi_bins_w0_p)

cdef double lps_w0_p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_w0_p, cmbps_w0_p)
            

cdef double[:] cosm_par_w0_m
cdef double C_w0_m
cdef double[:, :] a_data_w0_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_w0_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_w0_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_w0_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_w0_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_w0_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_w0_m, C_w0_m, a_data_w0_m, b_data_w0_m, c_data_w0_m, lps_data_w0_m, scale_factor_data_w0_m, window_data_w0_m, mps_data_w0_m, z_at_chi_data_w0_m, cmbps_w0_m, galaxy_density_chi_bins_w0_m = data_import_func('data_w0_m')

cdef double lbs_w0_m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_w0_m, a_data_w0_m, b_data_w0_m, c_data_w0_m, scale_factor_data_w0_m, window_data_w0_m, mps_data_w0_m, z_at_chi_data_w0_m, galaxy_density_chi_bins_w0_m)

cdef double lps_w0_m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_w0_m, cmbps_w0_m)
            

cdef double[:] cosm_par_logT_AGN_p
cdef double C_logT_AGN_p
cdef double[:, :] a_data_logT_AGN_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_logT_AGN_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_logT_AGN_p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_logT_AGN_p = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_logT_AGN_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_logT_AGN_p = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_logT_AGN_p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_logT_AGN_p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_logT_AGN_p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_logT_AGN_p = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_logT_AGN_p, C_logT_AGN_p, a_data_logT_AGN_p, b_data_logT_AGN_p, c_data_logT_AGN_p, lps_data_logT_AGN_p, scale_factor_data_logT_AGN_p, window_data_logT_AGN_p, mps_data_logT_AGN_p, z_at_chi_data_logT_AGN_p, cmbps_logT_AGN_p, galaxy_density_chi_bins_logT_AGN_p = data_import_func('data_logT_AGN_p')

cdef double lbs_logT_AGN_p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_logT_AGN_p, a_data_logT_AGN_p, b_data_logT_AGN_p, c_data_logT_AGN_p, scale_factor_data_logT_AGN_p, window_data_logT_AGN_p, mps_data_logT_AGN_p, z_at_chi_data_logT_AGN_p, galaxy_density_chi_bins_logT_AGN_p)

cdef double lps_logT_AGN_p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_logT_AGN_p, cmbps_logT_AGN_p)
            

cdef double[:] cosm_par_logT_AGN_m
cdef double C_logT_AGN_m
cdef double[:, :] a_data_logT_AGN_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_logT_AGN_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_logT_AGN_m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_data_logT_AGN_m = np.zeros((15, k_num), dtype=np.float64)
cdef double[:] scale_factor_data_logT_AGN_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_data_logT_AGN_m = np.zeros((5, chi_num), dtype=np.float64)
cdef double[:, :] mps_data_logT_AGN_m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_logT_AGN_m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_logT_AGN_m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)
cdef double[:, :] galaxy_density_chi_bins_logT_AGN_m = np.zeros((4, chi_num), dtype=np.float64)

cosm_par_logT_AGN_m, C_logT_AGN_m, a_data_logT_AGN_m, b_data_logT_AGN_m, c_data_logT_AGN_m, lps_data_logT_AGN_m, scale_factor_data_logT_AGN_m, window_data_logT_AGN_m, mps_data_logT_AGN_m, z_at_chi_data_logT_AGN_m, cmbps_logT_AGN_m, galaxy_density_chi_bins_logT_AGN_m = data_import_func('data_logT_AGN_m')

cdef double lbs_logT_AGN_m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_logT_AGN_m, a_data_logT_AGN_m, b_data_logT_AGN_m, c_data_logT_AGN_m, scale_factor_data_logT_AGN_m, window_data_logT_AGN_m, mps_data_logT_AGN_m, z_at_chi_data_logT_AGN_m, galaxy_density_chi_bins_logT_AGN_m)

cdef double lps_logT_AGN_m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_data_logT_AGN_m, cmbps_logT_AGN_m)
            
