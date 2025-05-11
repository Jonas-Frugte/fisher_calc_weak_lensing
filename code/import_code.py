

cdef double[:] cosm_par_H_p_2m
cdef double C_H_p_2m
cdef double[:, :] a_data_H_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_2m, C_H_p_2m, a_data_H_p_2m, b_data_H_p_2m, c_data_H_p_2m, lps_cc_data_H_p_2m, lps_cs_data_H_p_2m, lps_ss_data_H_p_2m, scale_factor_data_H_p_2m, window_c_data_H_p_2m, window_s_data_H_p_2m, mps_data_H_p_2m, z_at_chi_data_H_p_2m, cmbps_H_p_2m = data_import_func('data_H_p_2m')

cdef double lbs_H_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_2m, a_data_H_p_2m, b_data_H_p_2m, c_data_H_p_2m, scale_factor_data_H_p_2m, window_c_data_H_p_2m, window_s_data_H_p_2m, mps_data_H_p_2m, z_at_chi_data_H_p_2m)

cdef double lps_H_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_2m, lps_cs_data_H_p_2m, lps_ss_data_H_p_2m, cmbps_H_p_2m)
            

cdef double[:] cosm_par_H_p_1m
cdef double C_H_p_1m
cdef double[:, :] a_data_H_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_1m, C_H_p_1m, a_data_H_p_1m, b_data_H_p_1m, c_data_H_p_1m, lps_cc_data_H_p_1m, lps_cs_data_H_p_1m, lps_ss_data_H_p_1m, scale_factor_data_H_p_1m, window_c_data_H_p_1m, window_s_data_H_p_1m, mps_data_H_p_1m, z_at_chi_data_H_p_1m, cmbps_H_p_1m = data_import_func('data_H_p_1m')

cdef double lbs_H_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_1m, a_data_H_p_1m, b_data_H_p_1m, c_data_H_p_1m, scale_factor_data_H_p_1m, window_c_data_H_p_1m, window_s_data_H_p_1m, mps_data_H_p_1m, z_at_chi_data_H_p_1m)

cdef double lps_H_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_1m, lps_cs_data_H_p_1m, lps_ss_data_H_p_1m, cmbps_H_p_1m)
            

cdef double[:] cosm_par_H_p_0
cdef double C_H_p_0
cdef double[:, :] a_data_H_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_0, C_H_p_0, a_data_H_p_0, b_data_H_p_0, c_data_H_p_0, lps_cc_data_H_p_0, lps_cs_data_H_p_0, lps_ss_data_H_p_0, scale_factor_data_H_p_0, window_c_data_H_p_0, window_s_data_H_p_0, mps_data_H_p_0, z_at_chi_data_H_p_0, cmbps_H_p_0 = data_import_func('data_H_p_0')

cdef double lbs_H_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_0, a_data_H_p_0, b_data_H_p_0, c_data_H_p_0, scale_factor_data_H_p_0, window_c_data_H_p_0, window_s_data_H_p_0, mps_data_H_p_0, z_at_chi_data_H_p_0)

cdef double lps_H_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_0, lps_cs_data_H_p_0, lps_ss_data_H_p_0, cmbps_H_p_0)
            

cdef double[:] cosm_par_H_p_1p
cdef double C_H_p_1p
cdef double[:, :] a_data_H_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_1p, C_H_p_1p, a_data_H_p_1p, b_data_H_p_1p, c_data_H_p_1p, lps_cc_data_H_p_1p, lps_cs_data_H_p_1p, lps_ss_data_H_p_1p, scale_factor_data_H_p_1p, window_c_data_H_p_1p, window_s_data_H_p_1p, mps_data_H_p_1p, z_at_chi_data_H_p_1p, cmbps_H_p_1p = data_import_func('data_H_p_1p')

cdef double lbs_H_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_1p, a_data_H_p_1p, b_data_H_p_1p, c_data_H_p_1p, scale_factor_data_H_p_1p, window_c_data_H_p_1p, window_s_data_H_p_1p, mps_data_H_p_1p, z_at_chi_data_H_p_1p)

cdef double lps_H_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_1p, lps_cs_data_H_p_1p, lps_ss_data_H_p_1p, cmbps_H_p_1p)
            

cdef double[:] cosm_par_H_p_2p
cdef double C_H_p_2p
cdef double[:, :] a_data_H_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_p_2p, C_H_p_2p, a_data_H_p_2p, b_data_H_p_2p, c_data_H_p_2p, lps_cc_data_H_p_2p, lps_cs_data_H_p_2p, lps_ss_data_H_p_2p, scale_factor_data_H_p_2p, window_c_data_H_p_2p, window_s_data_H_p_2p, mps_data_H_p_2p, z_at_chi_data_H_p_2p, cmbps_H_p_2p = data_import_func('data_H_p_2p')

cdef double lbs_H_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_p_2p, a_data_H_p_2p, b_data_H_p_2p, c_data_H_p_2p, scale_factor_data_H_p_2p, window_c_data_H_p_2p, window_s_data_H_p_2p, mps_data_H_p_2p, z_at_chi_data_H_p_2p)

cdef double lps_H_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_2p, lps_cs_data_H_p_2p, lps_ss_data_H_p_2p, cmbps_H_p_2p)
            

cdef double[:] cosm_par_H_m_2m
cdef double C_H_m_2m
cdef double[:, :] a_data_H_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_2m, C_H_m_2m, a_data_H_m_2m, b_data_H_m_2m, c_data_H_m_2m, lps_cc_data_H_m_2m, lps_cs_data_H_m_2m, lps_ss_data_H_m_2m, scale_factor_data_H_m_2m, window_c_data_H_m_2m, window_s_data_H_m_2m, mps_data_H_m_2m, z_at_chi_data_H_m_2m, cmbps_H_m_2m = data_import_func('data_H_m_2m')

cdef double lbs_H_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_2m, a_data_H_m_2m, b_data_H_m_2m, c_data_H_m_2m, scale_factor_data_H_m_2m, window_c_data_H_m_2m, window_s_data_H_m_2m, mps_data_H_m_2m, z_at_chi_data_H_m_2m)

cdef double lps_H_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_2m, lps_cs_data_H_m_2m, lps_ss_data_H_m_2m, cmbps_H_m_2m)
            

cdef double[:] cosm_par_H_m_1m
cdef double C_H_m_1m
cdef double[:, :] a_data_H_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_1m, C_H_m_1m, a_data_H_m_1m, b_data_H_m_1m, c_data_H_m_1m, lps_cc_data_H_m_1m, lps_cs_data_H_m_1m, lps_ss_data_H_m_1m, scale_factor_data_H_m_1m, window_c_data_H_m_1m, window_s_data_H_m_1m, mps_data_H_m_1m, z_at_chi_data_H_m_1m, cmbps_H_m_1m = data_import_func('data_H_m_1m')

cdef double lbs_H_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_1m, a_data_H_m_1m, b_data_H_m_1m, c_data_H_m_1m, scale_factor_data_H_m_1m, window_c_data_H_m_1m, window_s_data_H_m_1m, mps_data_H_m_1m, z_at_chi_data_H_m_1m)

cdef double lps_H_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_1m, lps_cs_data_H_m_1m, lps_ss_data_H_m_1m, cmbps_H_m_1m)
            

cdef double[:] cosm_par_H_m_0
cdef double C_H_m_0
cdef double[:, :] a_data_H_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_0, C_H_m_0, a_data_H_m_0, b_data_H_m_0, c_data_H_m_0, lps_cc_data_H_m_0, lps_cs_data_H_m_0, lps_ss_data_H_m_0, scale_factor_data_H_m_0, window_c_data_H_m_0, window_s_data_H_m_0, mps_data_H_m_0, z_at_chi_data_H_m_0, cmbps_H_m_0 = data_import_func('data_H_m_0')

cdef double lbs_H_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_0, a_data_H_m_0, b_data_H_m_0, c_data_H_m_0, scale_factor_data_H_m_0, window_c_data_H_m_0, window_s_data_H_m_0, mps_data_H_m_0, z_at_chi_data_H_m_0)

cdef double lps_H_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_0, lps_cs_data_H_m_0, lps_ss_data_H_m_0, cmbps_H_m_0)
            

cdef double[:] cosm_par_H_m_1p
cdef double C_H_m_1p
cdef double[:, :] a_data_H_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_1p, C_H_m_1p, a_data_H_m_1p, b_data_H_m_1p, c_data_H_m_1p, lps_cc_data_H_m_1p, lps_cs_data_H_m_1p, lps_ss_data_H_m_1p, scale_factor_data_H_m_1p, window_c_data_H_m_1p, window_s_data_H_m_1p, mps_data_H_m_1p, z_at_chi_data_H_m_1p, cmbps_H_m_1p = data_import_func('data_H_m_1p')

cdef double lbs_H_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_1p, a_data_H_m_1p, b_data_H_m_1p, c_data_H_m_1p, scale_factor_data_H_m_1p, window_c_data_H_m_1p, window_s_data_H_m_1p, mps_data_H_m_1p, z_at_chi_data_H_m_1p)

cdef double lps_H_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_1p, lps_cs_data_H_m_1p, lps_ss_data_H_m_1p, cmbps_H_m_1p)
            

cdef double[:] cosm_par_H_m_2p
cdef double C_H_m_2p
cdef double[:, :] a_data_H_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_H_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_H_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_H_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_H_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_H_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_H_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_H_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_H_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_H_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_H_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_H_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_H_m_2p, C_H_m_2p, a_data_H_m_2p, b_data_H_m_2p, c_data_H_m_2p, lps_cc_data_H_m_2p, lps_cs_data_H_m_2p, lps_ss_data_H_m_2p, scale_factor_data_H_m_2p, window_c_data_H_m_2p, window_s_data_H_m_2p, mps_data_H_m_2p, z_at_chi_data_H_m_2p, cmbps_H_m_2p = data_import_func('data_H_m_2p')

cdef double lbs_H_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_H_m_2p, a_data_H_m_2p, b_data_H_m_2p, c_data_H_m_2p, scale_factor_data_H_m_2p, window_c_data_H_m_2p, window_s_data_H_m_2p, mps_data_H_m_2p, z_at_chi_data_H_m_2p)

cdef double lps_H_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_2p, lps_cs_data_H_m_2p, lps_ss_data_H_m_2p, cmbps_H_m_2p)
            

cdef double[:] cosm_par_ombh2_p_2m
cdef double C_ombh2_p_2m
cdef double[:, :] a_data_ombh2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_2m, C_ombh2_p_2m, a_data_ombh2_p_2m, b_data_ombh2_p_2m, c_data_ombh2_p_2m, lps_cc_data_ombh2_p_2m, lps_cs_data_ombh2_p_2m, lps_ss_data_ombh2_p_2m, scale_factor_data_ombh2_p_2m, window_c_data_ombh2_p_2m, window_s_data_ombh2_p_2m, mps_data_ombh2_p_2m, z_at_chi_data_ombh2_p_2m, cmbps_ombh2_p_2m = data_import_func('data_ombh2_p_2m')

cdef double lbs_ombh2_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_2m, a_data_ombh2_p_2m, b_data_ombh2_p_2m, c_data_ombh2_p_2m, scale_factor_data_ombh2_p_2m, window_c_data_ombh2_p_2m, window_s_data_ombh2_p_2m, mps_data_ombh2_p_2m, z_at_chi_data_ombh2_p_2m)

cdef double lps_ombh2_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_2m, lps_cs_data_ombh2_p_2m, lps_ss_data_ombh2_p_2m, cmbps_ombh2_p_2m)
            

cdef double[:] cosm_par_ombh2_p_1m
cdef double C_ombh2_p_1m
cdef double[:, :] a_data_ombh2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_1m, C_ombh2_p_1m, a_data_ombh2_p_1m, b_data_ombh2_p_1m, c_data_ombh2_p_1m, lps_cc_data_ombh2_p_1m, lps_cs_data_ombh2_p_1m, lps_ss_data_ombh2_p_1m, scale_factor_data_ombh2_p_1m, window_c_data_ombh2_p_1m, window_s_data_ombh2_p_1m, mps_data_ombh2_p_1m, z_at_chi_data_ombh2_p_1m, cmbps_ombh2_p_1m = data_import_func('data_ombh2_p_1m')

cdef double lbs_ombh2_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_1m, a_data_ombh2_p_1m, b_data_ombh2_p_1m, c_data_ombh2_p_1m, scale_factor_data_ombh2_p_1m, window_c_data_ombh2_p_1m, window_s_data_ombh2_p_1m, mps_data_ombh2_p_1m, z_at_chi_data_ombh2_p_1m)

cdef double lps_ombh2_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_1m, lps_cs_data_ombh2_p_1m, lps_ss_data_ombh2_p_1m, cmbps_ombh2_p_1m)
            

cdef double[:] cosm_par_ombh2_p_0
cdef double C_ombh2_p_0
cdef double[:, :] a_data_ombh2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_0, C_ombh2_p_0, a_data_ombh2_p_0, b_data_ombh2_p_0, c_data_ombh2_p_0, lps_cc_data_ombh2_p_0, lps_cs_data_ombh2_p_0, lps_ss_data_ombh2_p_0, scale_factor_data_ombh2_p_0, window_c_data_ombh2_p_0, window_s_data_ombh2_p_0, mps_data_ombh2_p_0, z_at_chi_data_ombh2_p_0, cmbps_ombh2_p_0 = data_import_func('data_ombh2_p_0')

cdef double lbs_ombh2_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_0, a_data_ombh2_p_0, b_data_ombh2_p_0, c_data_ombh2_p_0, scale_factor_data_ombh2_p_0, window_c_data_ombh2_p_0, window_s_data_ombh2_p_0, mps_data_ombh2_p_0, z_at_chi_data_ombh2_p_0)

cdef double lps_ombh2_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_0, lps_cs_data_ombh2_p_0, lps_ss_data_ombh2_p_0, cmbps_ombh2_p_0)
            

cdef double[:] cosm_par_ombh2_p_1p
cdef double C_ombh2_p_1p
cdef double[:, :] a_data_ombh2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_1p, C_ombh2_p_1p, a_data_ombh2_p_1p, b_data_ombh2_p_1p, c_data_ombh2_p_1p, lps_cc_data_ombh2_p_1p, lps_cs_data_ombh2_p_1p, lps_ss_data_ombh2_p_1p, scale_factor_data_ombh2_p_1p, window_c_data_ombh2_p_1p, window_s_data_ombh2_p_1p, mps_data_ombh2_p_1p, z_at_chi_data_ombh2_p_1p, cmbps_ombh2_p_1p = data_import_func('data_ombh2_p_1p')

cdef double lbs_ombh2_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_1p, a_data_ombh2_p_1p, b_data_ombh2_p_1p, c_data_ombh2_p_1p, scale_factor_data_ombh2_p_1p, window_c_data_ombh2_p_1p, window_s_data_ombh2_p_1p, mps_data_ombh2_p_1p, z_at_chi_data_ombh2_p_1p)

cdef double lps_ombh2_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_1p, lps_cs_data_ombh2_p_1p, lps_ss_data_ombh2_p_1p, cmbps_ombh2_p_1p)
            

cdef double[:] cosm_par_ombh2_p_2p
cdef double C_ombh2_p_2p
cdef double[:, :] a_data_ombh2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_p_2p, C_ombh2_p_2p, a_data_ombh2_p_2p, b_data_ombh2_p_2p, c_data_ombh2_p_2p, lps_cc_data_ombh2_p_2p, lps_cs_data_ombh2_p_2p, lps_ss_data_ombh2_p_2p, scale_factor_data_ombh2_p_2p, window_c_data_ombh2_p_2p, window_s_data_ombh2_p_2p, mps_data_ombh2_p_2p, z_at_chi_data_ombh2_p_2p, cmbps_ombh2_p_2p = data_import_func('data_ombh2_p_2p')

cdef double lbs_ombh2_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_p_2p, a_data_ombh2_p_2p, b_data_ombh2_p_2p, c_data_ombh2_p_2p, scale_factor_data_ombh2_p_2p, window_c_data_ombh2_p_2p, window_s_data_ombh2_p_2p, mps_data_ombh2_p_2p, z_at_chi_data_ombh2_p_2p)

cdef double lps_ombh2_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_2p, lps_cs_data_ombh2_p_2p, lps_ss_data_ombh2_p_2p, cmbps_ombh2_p_2p)
            

cdef double[:] cosm_par_ombh2_m_2m
cdef double C_ombh2_m_2m
cdef double[:, :] a_data_ombh2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_2m, C_ombh2_m_2m, a_data_ombh2_m_2m, b_data_ombh2_m_2m, c_data_ombh2_m_2m, lps_cc_data_ombh2_m_2m, lps_cs_data_ombh2_m_2m, lps_ss_data_ombh2_m_2m, scale_factor_data_ombh2_m_2m, window_c_data_ombh2_m_2m, window_s_data_ombh2_m_2m, mps_data_ombh2_m_2m, z_at_chi_data_ombh2_m_2m, cmbps_ombh2_m_2m = data_import_func('data_ombh2_m_2m')

cdef double lbs_ombh2_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_2m, a_data_ombh2_m_2m, b_data_ombh2_m_2m, c_data_ombh2_m_2m, scale_factor_data_ombh2_m_2m, window_c_data_ombh2_m_2m, window_s_data_ombh2_m_2m, mps_data_ombh2_m_2m, z_at_chi_data_ombh2_m_2m)

cdef double lps_ombh2_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_2m, lps_cs_data_ombh2_m_2m, lps_ss_data_ombh2_m_2m, cmbps_ombh2_m_2m)
            

cdef double[:] cosm_par_ombh2_m_1m
cdef double C_ombh2_m_1m
cdef double[:, :] a_data_ombh2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_1m, C_ombh2_m_1m, a_data_ombh2_m_1m, b_data_ombh2_m_1m, c_data_ombh2_m_1m, lps_cc_data_ombh2_m_1m, lps_cs_data_ombh2_m_1m, lps_ss_data_ombh2_m_1m, scale_factor_data_ombh2_m_1m, window_c_data_ombh2_m_1m, window_s_data_ombh2_m_1m, mps_data_ombh2_m_1m, z_at_chi_data_ombh2_m_1m, cmbps_ombh2_m_1m = data_import_func('data_ombh2_m_1m')

cdef double lbs_ombh2_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_1m, a_data_ombh2_m_1m, b_data_ombh2_m_1m, c_data_ombh2_m_1m, scale_factor_data_ombh2_m_1m, window_c_data_ombh2_m_1m, window_s_data_ombh2_m_1m, mps_data_ombh2_m_1m, z_at_chi_data_ombh2_m_1m)

cdef double lps_ombh2_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_1m, lps_cs_data_ombh2_m_1m, lps_ss_data_ombh2_m_1m, cmbps_ombh2_m_1m)
            

cdef double[:] cosm_par_ombh2_m_0
cdef double C_ombh2_m_0
cdef double[:, :] a_data_ombh2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_0, C_ombh2_m_0, a_data_ombh2_m_0, b_data_ombh2_m_0, c_data_ombh2_m_0, lps_cc_data_ombh2_m_0, lps_cs_data_ombh2_m_0, lps_ss_data_ombh2_m_0, scale_factor_data_ombh2_m_0, window_c_data_ombh2_m_0, window_s_data_ombh2_m_0, mps_data_ombh2_m_0, z_at_chi_data_ombh2_m_0, cmbps_ombh2_m_0 = data_import_func('data_ombh2_m_0')

cdef double lbs_ombh2_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_0, a_data_ombh2_m_0, b_data_ombh2_m_0, c_data_ombh2_m_0, scale_factor_data_ombh2_m_0, window_c_data_ombh2_m_0, window_s_data_ombh2_m_0, mps_data_ombh2_m_0, z_at_chi_data_ombh2_m_0)

cdef double lps_ombh2_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_0, lps_cs_data_ombh2_m_0, lps_ss_data_ombh2_m_0, cmbps_ombh2_m_0)
            

cdef double[:] cosm_par_ombh2_m_1p
cdef double C_ombh2_m_1p
cdef double[:, :] a_data_ombh2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_1p, C_ombh2_m_1p, a_data_ombh2_m_1p, b_data_ombh2_m_1p, c_data_ombh2_m_1p, lps_cc_data_ombh2_m_1p, lps_cs_data_ombh2_m_1p, lps_ss_data_ombh2_m_1p, scale_factor_data_ombh2_m_1p, window_c_data_ombh2_m_1p, window_s_data_ombh2_m_1p, mps_data_ombh2_m_1p, z_at_chi_data_ombh2_m_1p, cmbps_ombh2_m_1p = data_import_func('data_ombh2_m_1p')

cdef double lbs_ombh2_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_1p, a_data_ombh2_m_1p, b_data_ombh2_m_1p, c_data_ombh2_m_1p, scale_factor_data_ombh2_m_1p, window_c_data_ombh2_m_1p, window_s_data_ombh2_m_1p, mps_data_ombh2_m_1p, z_at_chi_data_ombh2_m_1p)

cdef double lps_ombh2_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_1p, lps_cs_data_ombh2_m_1p, lps_ss_data_ombh2_m_1p, cmbps_ombh2_m_1p)
            

cdef double[:] cosm_par_ombh2_m_2p
cdef double C_ombh2_m_2p
cdef double[:, :] a_data_ombh2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ombh2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ombh2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ombh2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ombh2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ombh2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ombh2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ombh2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ombh2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ombh2_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ombh2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ombh2_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ombh2_m_2p, C_ombh2_m_2p, a_data_ombh2_m_2p, b_data_ombh2_m_2p, c_data_ombh2_m_2p, lps_cc_data_ombh2_m_2p, lps_cs_data_ombh2_m_2p, lps_ss_data_ombh2_m_2p, scale_factor_data_ombh2_m_2p, window_c_data_ombh2_m_2p, window_s_data_ombh2_m_2p, mps_data_ombh2_m_2p, z_at_chi_data_ombh2_m_2p, cmbps_ombh2_m_2p = data_import_func('data_ombh2_m_2p')

cdef double lbs_ombh2_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ombh2_m_2p, a_data_ombh2_m_2p, b_data_ombh2_m_2p, c_data_ombh2_m_2p, scale_factor_data_ombh2_m_2p, window_c_data_ombh2_m_2p, window_s_data_ombh2_m_2p, mps_data_ombh2_m_2p, z_at_chi_data_ombh2_m_2p)

cdef double lps_ombh2_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_2p, lps_cs_data_ombh2_m_2p, lps_ss_data_ombh2_m_2p, cmbps_ombh2_m_2p)
            

cdef double[:] cosm_par_omch2_p_2m
cdef double C_omch2_p_2m
cdef double[:, :] a_data_omch2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_2m, C_omch2_p_2m, a_data_omch2_p_2m, b_data_omch2_p_2m, c_data_omch2_p_2m, lps_cc_data_omch2_p_2m, lps_cs_data_omch2_p_2m, lps_ss_data_omch2_p_2m, scale_factor_data_omch2_p_2m, window_c_data_omch2_p_2m, window_s_data_omch2_p_2m, mps_data_omch2_p_2m, z_at_chi_data_omch2_p_2m, cmbps_omch2_p_2m = data_import_func('data_omch2_p_2m')

cdef double lbs_omch2_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_2m, a_data_omch2_p_2m, b_data_omch2_p_2m, c_data_omch2_p_2m, scale_factor_data_omch2_p_2m, window_c_data_omch2_p_2m, window_s_data_omch2_p_2m, mps_data_omch2_p_2m, z_at_chi_data_omch2_p_2m)

cdef double lps_omch2_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_2m, lps_cs_data_omch2_p_2m, lps_ss_data_omch2_p_2m, cmbps_omch2_p_2m)
            

cdef double[:] cosm_par_omch2_p_1m
cdef double C_omch2_p_1m
cdef double[:, :] a_data_omch2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_1m, C_omch2_p_1m, a_data_omch2_p_1m, b_data_omch2_p_1m, c_data_omch2_p_1m, lps_cc_data_omch2_p_1m, lps_cs_data_omch2_p_1m, lps_ss_data_omch2_p_1m, scale_factor_data_omch2_p_1m, window_c_data_omch2_p_1m, window_s_data_omch2_p_1m, mps_data_omch2_p_1m, z_at_chi_data_omch2_p_1m, cmbps_omch2_p_1m = data_import_func('data_omch2_p_1m')

cdef double lbs_omch2_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_1m, a_data_omch2_p_1m, b_data_omch2_p_1m, c_data_omch2_p_1m, scale_factor_data_omch2_p_1m, window_c_data_omch2_p_1m, window_s_data_omch2_p_1m, mps_data_omch2_p_1m, z_at_chi_data_omch2_p_1m)

cdef double lps_omch2_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_1m, lps_cs_data_omch2_p_1m, lps_ss_data_omch2_p_1m, cmbps_omch2_p_1m)
            

cdef double[:] cosm_par_omch2_p_0
cdef double C_omch2_p_0
cdef double[:, :] a_data_omch2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_0, C_omch2_p_0, a_data_omch2_p_0, b_data_omch2_p_0, c_data_omch2_p_0, lps_cc_data_omch2_p_0, lps_cs_data_omch2_p_0, lps_ss_data_omch2_p_0, scale_factor_data_omch2_p_0, window_c_data_omch2_p_0, window_s_data_omch2_p_0, mps_data_omch2_p_0, z_at_chi_data_omch2_p_0, cmbps_omch2_p_0 = data_import_func('data_omch2_p_0')

cdef double lbs_omch2_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_0, a_data_omch2_p_0, b_data_omch2_p_0, c_data_omch2_p_0, scale_factor_data_omch2_p_0, window_c_data_omch2_p_0, window_s_data_omch2_p_0, mps_data_omch2_p_0, z_at_chi_data_omch2_p_0)

cdef double lps_omch2_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_0, lps_cs_data_omch2_p_0, lps_ss_data_omch2_p_0, cmbps_omch2_p_0)
            

cdef double[:] cosm_par_omch2_p_1p
cdef double C_omch2_p_1p
cdef double[:, :] a_data_omch2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_1p, C_omch2_p_1p, a_data_omch2_p_1p, b_data_omch2_p_1p, c_data_omch2_p_1p, lps_cc_data_omch2_p_1p, lps_cs_data_omch2_p_1p, lps_ss_data_omch2_p_1p, scale_factor_data_omch2_p_1p, window_c_data_omch2_p_1p, window_s_data_omch2_p_1p, mps_data_omch2_p_1p, z_at_chi_data_omch2_p_1p, cmbps_omch2_p_1p = data_import_func('data_omch2_p_1p')

cdef double lbs_omch2_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_1p, a_data_omch2_p_1p, b_data_omch2_p_1p, c_data_omch2_p_1p, scale_factor_data_omch2_p_1p, window_c_data_omch2_p_1p, window_s_data_omch2_p_1p, mps_data_omch2_p_1p, z_at_chi_data_omch2_p_1p)

cdef double lps_omch2_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_1p, lps_cs_data_omch2_p_1p, lps_ss_data_omch2_p_1p, cmbps_omch2_p_1p)
            

cdef double[:] cosm_par_omch2_p_2p
cdef double C_omch2_p_2p
cdef double[:, :] a_data_omch2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_p_2p, C_omch2_p_2p, a_data_omch2_p_2p, b_data_omch2_p_2p, c_data_omch2_p_2p, lps_cc_data_omch2_p_2p, lps_cs_data_omch2_p_2p, lps_ss_data_omch2_p_2p, scale_factor_data_omch2_p_2p, window_c_data_omch2_p_2p, window_s_data_omch2_p_2p, mps_data_omch2_p_2p, z_at_chi_data_omch2_p_2p, cmbps_omch2_p_2p = data_import_func('data_omch2_p_2p')

cdef double lbs_omch2_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_p_2p, a_data_omch2_p_2p, b_data_omch2_p_2p, c_data_omch2_p_2p, scale_factor_data_omch2_p_2p, window_c_data_omch2_p_2p, window_s_data_omch2_p_2p, mps_data_omch2_p_2p, z_at_chi_data_omch2_p_2p)

cdef double lps_omch2_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_2p, lps_cs_data_omch2_p_2p, lps_ss_data_omch2_p_2p, cmbps_omch2_p_2p)
            

cdef double[:] cosm_par_omch2_m_2m
cdef double C_omch2_m_2m
cdef double[:, :] a_data_omch2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_2m, C_omch2_m_2m, a_data_omch2_m_2m, b_data_omch2_m_2m, c_data_omch2_m_2m, lps_cc_data_omch2_m_2m, lps_cs_data_omch2_m_2m, lps_ss_data_omch2_m_2m, scale_factor_data_omch2_m_2m, window_c_data_omch2_m_2m, window_s_data_omch2_m_2m, mps_data_omch2_m_2m, z_at_chi_data_omch2_m_2m, cmbps_omch2_m_2m = data_import_func('data_omch2_m_2m')

cdef double lbs_omch2_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_2m, a_data_omch2_m_2m, b_data_omch2_m_2m, c_data_omch2_m_2m, scale_factor_data_omch2_m_2m, window_c_data_omch2_m_2m, window_s_data_omch2_m_2m, mps_data_omch2_m_2m, z_at_chi_data_omch2_m_2m)

cdef double lps_omch2_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_2m, lps_cs_data_omch2_m_2m, lps_ss_data_omch2_m_2m, cmbps_omch2_m_2m)
            

cdef double[:] cosm_par_omch2_m_1m
cdef double C_omch2_m_1m
cdef double[:, :] a_data_omch2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_1m, C_omch2_m_1m, a_data_omch2_m_1m, b_data_omch2_m_1m, c_data_omch2_m_1m, lps_cc_data_omch2_m_1m, lps_cs_data_omch2_m_1m, lps_ss_data_omch2_m_1m, scale_factor_data_omch2_m_1m, window_c_data_omch2_m_1m, window_s_data_omch2_m_1m, mps_data_omch2_m_1m, z_at_chi_data_omch2_m_1m, cmbps_omch2_m_1m = data_import_func('data_omch2_m_1m')

cdef double lbs_omch2_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_1m, a_data_omch2_m_1m, b_data_omch2_m_1m, c_data_omch2_m_1m, scale_factor_data_omch2_m_1m, window_c_data_omch2_m_1m, window_s_data_omch2_m_1m, mps_data_omch2_m_1m, z_at_chi_data_omch2_m_1m)

cdef double lps_omch2_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_1m, lps_cs_data_omch2_m_1m, lps_ss_data_omch2_m_1m, cmbps_omch2_m_1m)
            

cdef double[:] cosm_par_omch2_m_0
cdef double C_omch2_m_0
cdef double[:, :] a_data_omch2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_0, C_omch2_m_0, a_data_omch2_m_0, b_data_omch2_m_0, c_data_omch2_m_0, lps_cc_data_omch2_m_0, lps_cs_data_omch2_m_0, lps_ss_data_omch2_m_0, scale_factor_data_omch2_m_0, window_c_data_omch2_m_0, window_s_data_omch2_m_0, mps_data_omch2_m_0, z_at_chi_data_omch2_m_0, cmbps_omch2_m_0 = data_import_func('data_omch2_m_0')

cdef double lbs_omch2_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_0, a_data_omch2_m_0, b_data_omch2_m_0, c_data_omch2_m_0, scale_factor_data_omch2_m_0, window_c_data_omch2_m_0, window_s_data_omch2_m_0, mps_data_omch2_m_0, z_at_chi_data_omch2_m_0)

cdef double lps_omch2_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_0, lps_cs_data_omch2_m_0, lps_ss_data_omch2_m_0, cmbps_omch2_m_0)
            

cdef double[:] cosm_par_omch2_m_1p
cdef double C_omch2_m_1p
cdef double[:, :] a_data_omch2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_1p, C_omch2_m_1p, a_data_omch2_m_1p, b_data_omch2_m_1p, c_data_omch2_m_1p, lps_cc_data_omch2_m_1p, lps_cs_data_omch2_m_1p, lps_ss_data_omch2_m_1p, scale_factor_data_omch2_m_1p, window_c_data_omch2_m_1p, window_s_data_omch2_m_1p, mps_data_omch2_m_1p, z_at_chi_data_omch2_m_1p, cmbps_omch2_m_1p = data_import_func('data_omch2_m_1p')

cdef double lbs_omch2_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_1p, a_data_omch2_m_1p, b_data_omch2_m_1p, c_data_omch2_m_1p, scale_factor_data_omch2_m_1p, window_c_data_omch2_m_1p, window_s_data_omch2_m_1p, mps_data_omch2_m_1p, z_at_chi_data_omch2_m_1p)

cdef double lps_omch2_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_1p, lps_cs_data_omch2_m_1p, lps_ss_data_omch2_m_1p, cmbps_omch2_m_1p)
            

cdef double[:] cosm_par_omch2_m_2p
cdef double C_omch2_m_2p
cdef double[:, :] a_data_omch2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_omch2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_omch2_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_omch2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_omch2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_omch2_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_omch2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_omch2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_omch2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_omch2_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_omch2_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_omch2_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_omch2_m_2p, C_omch2_m_2p, a_data_omch2_m_2p, b_data_omch2_m_2p, c_data_omch2_m_2p, lps_cc_data_omch2_m_2p, lps_cs_data_omch2_m_2p, lps_ss_data_omch2_m_2p, scale_factor_data_omch2_m_2p, window_c_data_omch2_m_2p, window_s_data_omch2_m_2p, mps_data_omch2_m_2p, z_at_chi_data_omch2_m_2p, cmbps_omch2_m_2p = data_import_func('data_omch2_m_2p')

cdef double lbs_omch2_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_omch2_m_2p, a_data_omch2_m_2p, b_data_omch2_m_2p, c_data_omch2_m_2p, scale_factor_data_omch2_m_2p, window_c_data_omch2_m_2p, window_s_data_omch2_m_2p, mps_data_omch2_m_2p, z_at_chi_data_omch2_m_2p)

cdef double lps_omch2_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_2p, lps_cs_data_omch2_m_2p, lps_ss_data_omch2_m_2p, cmbps_omch2_m_2p)
            

cdef double[:] cosm_par_ns_p_2m
cdef double C_ns_p_2m
cdef double[:, :] a_data_ns_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_2m, C_ns_p_2m, a_data_ns_p_2m, b_data_ns_p_2m, c_data_ns_p_2m, lps_cc_data_ns_p_2m, lps_cs_data_ns_p_2m, lps_ss_data_ns_p_2m, scale_factor_data_ns_p_2m, window_c_data_ns_p_2m, window_s_data_ns_p_2m, mps_data_ns_p_2m, z_at_chi_data_ns_p_2m, cmbps_ns_p_2m = data_import_func('data_ns_p_2m')

cdef double lbs_ns_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_2m, a_data_ns_p_2m, b_data_ns_p_2m, c_data_ns_p_2m, scale_factor_data_ns_p_2m, window_c_data_ns_p_2m, window_s_data_ns_p_2m, mps_data_ns_p_2m, z_at_chi_data_ns_p_2m)

cdef double lps_ns_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_2m, lps_cs_data_ns_p_2m, lps_ss_data_ns_p_2m, cmbps_ns_p_2m)
            

cdef double[:] cosm_par_ns_p_1m
cdef double C_ns_p_1m
cdef double[:, :] a_data_ns_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_1m, C_ns_p_1m, a_data_ns_p_1m, b_data_ns_p_1m, c_data_ns_p_1m, lps_cc_data_ns_p_1m, lps_cs_data_ns_p_1m, lps_ss_data_ns_p_1m, scale_factor_data_ns_p_1m, window_c_data_ns_p_1m, window_s_data_ns_p_1m, mps_data_ns_p_1m, z_at_chi_data_ns_p_1m, cmbps_ns_p_1m = data_import_func('data_ns_p_1m')

cdef double lbs_ns_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_1m, a_data_ns_p_1m, b_data_ns_p_1m, c_data_ns_p_1m, scale_factor_data_ns_p_1m, window_c_data_ns_p_1m, window_s_data_ns_p_1m, mps_data_ns_p_1m, z_at_chi_data_ns_p_1m)

cdef double lps_ns_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_1m, lps_cs_data_ns_p_1m, lps_ss_data_ns_p_1m, cmbps_ns_p_1m)
            

cdef double[:] cosm_par_ns_p_0
cdef double C_ns_p_0
cdef double[:, :] a_data_ns_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_0, C_ns_p_0, a_data_ns_p_0, b_data_ns_p_0, c_data_ns_p_0, lps_cc_data_ns_p_0, lps_cs_data_ns_p_0, lps_ss_data_ns_p_0, scale_factor_data_ns_p_0, window_c_data_ns_p_0, window_s_data_ns_p_0, mps_data_ns_p_0, z_at_chi_data_ns_p_0, cmbps_ns_p_0 = data_import_func('data_ns_p_0')

cdef double lbs_ns_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_0, a_data_ns_p_0, b_data_ns_p_0, c_data_ns_p_0, scale_factor_data_ns_p_0, window_c_data_ns_p_0, window_s_data_ns_p_0, mps_data_ns_p_0, z_at_chi_data_ns_p_0)

cdef double lps_ns_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_0, lps_cs_data_ns_p_0, lps_ss_data_ns_p_0, cmbps_ns_p_0)
            

cdef double[:] cosm_par_ns_p_1p
cdef double C_ns_p_1p
cdef double[:, :] a_data_ns_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_1p, C_ns_p_1p, a_data_ns_p_1p, b_data_ns_p_1p, c_data_ns_p_1p, lps_cc_data_ns_p_1p, lps_cs_data_ns_p_1p, lps_ss_data_ns_p_1p, scale_factor_data_ns_p_1p, window_c_data_ns_p_1p, window_s_data_ns_p_1p, mps_data_ns_p_1p, z_at_chi_data_ns_p_1p, cmbps_ns_p_1p = data_import_func('data_ns_p_1p')

cdef double lbs_ns_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_1p, a_data_ns_p_1p, b_data_ns_p_1p, c_data_ns_p_1p, scale_factor_data_ns_p_1p, window_c_data_ns_p_1p, window_s_data_ns_p_1p, mps_data_ns_p_1p, z_at_chi_data_ns_p_1p)

cdef double lps_ns_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_1p, lps_cs_data_ns_p_1p, lps_ss_data_ns_p_1p, cmbps_ns_p_1p)
            

cdef double[:] cosm_par_ns_p_2p
cdef double C_ns_p_2p
cdef double[:, :] a_data_ns_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_p_2p, C_ns_p_2p, a_data_ns_p_2p, b_data_ns_p_2p, c_data_ns_p_2p, lps_cc_data_ns_p_2p, lps_cs_data_ns_p_2p, lps_ss_data_ns_p_2p, scale_factor_data_ns_p_2p, window_c_data_ns_p_2p, window_s_data_ns_p_2p, mps_data_ns_p_2p, z_at_chi_data_ns_p_2p, cmbps_ns_p_2p = data_import_func('data_ns_p_2p')

cdef double lbs_ns_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_p_2p, a_data_ns_p_2p, b_data_ns_p_2p, c_data_ns_p_2p, scale_factor_data_ns_p_2p, window_c_data_ns_p_2p, window_s_data_ns_p_2p, mps_data_ns_p_2p, z_at_chi_data_ns_p_2p)

cdef double lps_ns_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_2p, lps_cs_data_ns_p_2p, lps_ss_data_ns_p_2p, cmbps_ns_p_2p)
            

cdef double[:] cosm_par_ns_m_2m
cdef double C_ns_m_2m
cdef double[:, :] a_data_ns_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_2m, C_ns_m_2m, a_data_ns_m_2m, b_data_ns_m_2m, c_data_ns_m_2m, lps_cc_data_ns_m_2m, lps_cs_data_ns_m_2m, lps_ss_data_ns_m_2m, scale_factor_data_ns_m_2m, window_c_data_ns_m_2m, window_s_data_ns_m_2m, mps_data_ns_m_2m, z_at_chi_data_ns_m_2m, cmbps_ns_m_2m = data_import_func('data_ns_m_2m')

cdef double lbs_ns_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_2m, a_data_ns_m_2m, b_data_ns_m_2m, c_data_ns_m_2m, scale_factor_data_ns_m_2m, window_c_data_ns_m_2m, window_s_data_ns_m_2m, mps_data_ns_m_2m, z_at_chi_data_ns_m_2m)

cdef double lps_ns_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_2m, lps_cs_data_ns_m_2m, lps_ss_data_ns_m_2m, cmbps_ns_m_2m)
            

cdef double[:] cosm_par_ns_m_1m
cdef double C_ns_m_1m
cdef double[:, :] a_data_ns_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_1m, C_ns_m_1m, a_data_ns_m_1m, b_data_ns_m_1m, c_data_ns_m_1m, lps_cc_data_ns_m_1m, lps_cs_data_ns_m_1m, lps_ss_data_ns_m_1m, scale_factor_data_ns_m_1m, window_c_data_ns_m_1m, window_s_data_ns_m_1m, mps_data_ns_m_1m, z_at_chi_data_ns_m_1m, cmbps_ns_m_1m = data_import_func('data_ns_m_1m')

cdef double lbs_ns_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_1m, a_data_ns_m_1m, b_data_ns_m_1m, c_data_ns_m_1m, scale_factor_data_ns_m_1m, window_c_data_ns_m_1m, window_s_data_ns_m_1m, mps_data_ns_m_1m, z_at_chi_data_ns_m_1m)

cdef double lps_ns_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_1m, lps_cs_data_ns_m_1m, lps_ss_data_ns_m_1m, cmbps_ns_m_1m)
            

cdef double[:] cosm_par_ns_m_0
cdef double C_ns_m_0
cdef double[:, :] a_data_ns_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_0, C_ns_m_0, a_data_ns_m_0, b_data_ns_m_0, c_data_ns_m_0, lps_cc_data_ns_m_0, lps_cs_data_ns_m_0, lps_ss_data_ns_m_0, scale_factor_data_ns_m_0, window_c_data_ns_m_0, window_s_data_ns_m_0, mps_data_ns_m_0, z_at_chi_data_ns_m_0, cmbps_ns_m_0 = data_import_func('data_ns_m_0')

cdef double lbs_ns_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_0, a_data_ns_m_0, b_data_ns_m_0, c_data_ns_m_0, scale_factor_data_ns_m_0, window_c_data_ns_m_0, window_s_data_ns_m_0, mps_data_ns_m_0, z_at_chi_data_ns_m_0)

cdef double lps_ns_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_0, lps_cs_data_ns_m_0, lps_ss_data_ns_m_0, cmbps_ns_m_0)
            

cdef double[:] cosm_par_ns_m_1p
cdef double C_ns_m_1p
cdef double[:, :] a_data_ns_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_1p, C_ns_m_1p, a_data_ns_m_1p, b_data_ns_m_1p, c_data_ns_m_1p, lps_cc_data_ns_m_1p, lps_cs_data_ns_m_1p, lps_ss_data_ns_m_1p, scale_factor_data_ns_m_1p, window_c_data_ns_m_1p, window_s_data_ns_m_1p, mps_data_ns_m_1p, z_at_chi_data_ns_m_1p, cmbps_ns_m_1p = data_import_func('data_ns_m_1p')

cdef double lbs_ns_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_1p, a_data_ns_m_1p, b_data_ns_m_1p, c_data_ns_m_1p, scale_factor_data_ns_m_1p, window_c_data_ns_m_1p, window_s_data_ns_m_1p, mps_data_ns_m_1p, z_at_chi_data_ns_m_1p)

cdef double lps_ns_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_1p, lps_cs_data_ns_m_1p, lps_ss_data_ns_m_1p, cmbps_ns_m_1p)
            

cdef double[:] cosm_par_ns_m_2p
cdef double C_ns_m_2p
cdef double[:, :] a_data_ns_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_ns_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_ns_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_ns_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_ns_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_ns_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_ns_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_ns_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_ns_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_ns_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_ns_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_ns_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_ns_m_2p, C_ns_m_2p, a_data_ns_m_2p, b_data_ns_m_2p, c_data_ns_m_2p, lps_cc_data_ns_m_2p, lps_cs_data_ns_m_2p, lps_ss_data_ns_m_2p, scale_factor_data_ns_m_2p, window_c_data_ns_m_2p, window_s_data_ns_m_2p, mps_data_ns_m_2p, z_at_chi_data_ns_m_2p, cmbps_ns_m_2p = data_import_func('data_ns_m_2p')

cdef double lbs_ns_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_ns_m_2p, a_data_ns_m_2p, b_data_ns_m_2p, c_data_ns_m_2p, scale_factor_data_ns_m_2p, window_c_data_ns_m_2p, window_s_data_ns_m_2p, mps_data_ns_m_2p, z_at_chi_data_ns_m_2p)

cdef double lps_ns_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_2p, lps_cs_data_ns_m_2p, lps_ss_data_ns_m_2p, cmbps_ns_m_2p)
            

cdef double[:] cosm_par_As_p_2m
cdef double C_As_p_2m
cdef double[:, :] a_data_As_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_2m, C_As_p_2m, a_data_As_p_2m, b_data_As_p_2m, c_data_As_p_2m, lps_cc_data_As_p_2m, lps_cs_data_As_p_2m, lps_ss_data_As_p_2m, scale_factor_data_As_p_2m, window_c_data_As_p_2m, window_s_data_As_p_2m, mps_data_As_p_2m, z_at_chi_data_As_p_2m, cmbps_As_p_2m = data_import_func('data_As_p_2m')

cdef double lbs_As_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_2m, a_data_As_p_2m, b_data_As_p_2m, c_data_As_p_2m, scale_factor_data_As_p_2m, window_c_data_As_p_2m, window_s_data_As_p_2m, mps_data_As_p_2m, z_at_chi_data_As_p_2m)

cdef double lps_As_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_2m, lps_cs_data_As_p_2m, lps_ss_data_As_p_2m, cmbps_As_p_2m)
            

cdef double[:] cosm_par_As_p_1m
cdef double C_As_p_1m
cdef double[:, :] a_data_As_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_1m, C_As_p_1m, a_data_As_p_1m, b_data_As_p_1m, c_data_As_p_1m, lps_cc_data_As_p_1m, lps_cs_data_As_p_1m, lps_ss_data_As_p_1m, scale_factor_data_As_p_1m, window_c_data_As_p_1m, window_s_data_As_p_1m, mps_data_As_p_1m, z_at_chi_data_As_p_1m, cmbps_As_p_1m = data_import_func('data_As_p_1m')

cdef double lbs_As_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_1m, a_data_As_p_1m, b_data_As_p_1m, c_data_As_p_1m, scale_factor_data_As_p_1m, window_c_data_As_p_1m, window_s_data_As_p_1m, mps_data_As_p_1m, z_at_chi_data_As_p_1m)

cdef double lps_As_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_1m, lps_cs_data_As_p_1m, lps_ss_data_As_p_1m, cmbps_As_p_1m)
            

cdef double[:] cosm_par_As_p_0
cdef double C_As_p_0
cdef double[:, :] a_data_As_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_0, C_As_p_0, a_data_As_p_0, b_data_As_p_0, c_data_As_p_0, lps_cc_data_As_p_0, lps_cs_data_As_p_0, lps_ss_data_As_p_0, scale_factor_data_As_p_0, window_c_data_As_p_0, window_s_data_As_p_0, mps_data_As_p_0, z_at_chi_data_As_p_0, cmbps_As_p_0 = data_import_func('data_As_p_0')

cdef double lbs_As_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_0, a_data_As_p_0, b_data_As_p_0, c_data_As_p_0, scale_factor_data_As_p_0, window_c_data_As_p_0, window_s_data_As_p_0, mps_data_As_p_0, z_at_chi_data_As_p_0)

cdef double lps_As_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_0, lps_cs_data_As_p_0, lps_ss_data_As_p_0, cmbps_As_p_0)
            

cdef double[:] cosm_par_As_p_1p
cdef double C_As_p_1p
cdef double[:, :] a_data_As_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_1p, C_As_p_1p, a_data_As_p_1p, b_data_As_p_1p, c_data_As_p_1p, lps_cc_data_As_p_1p, lps_cs_data_As_p_1p, lps_ss_data_As_p_1p, scale_factor_data_As_p_1p, window_c_data_As_p_1p, window_s_data_As_p_1p, mps_data_As_p_1p, z_at_chi_data_As_p_1p, cmbps_As_p_1p = data_import_func('data_As_p_1p')

cdef double lbs_As_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_1p, a_data_As_p_1p, b_data_As_p_1p, c_data_As_p_1p, scale_factor_data_As_p_1p, window_c_data_As_p_1p, window_s_data_As_p_1p, mps_data_As_p_1p, z_at_chi_data_As_p_1p)

cdef double lps_As_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_1p, lps_cs_data_As_p_1p, lps_ss_data_As_p_1p, cmbps_As_p_1p)
            

cdef double[:] cosm_par_As_p_2p
cdef double C_As_p_2p
cdef double[:, :] a_data_As_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_p_2p, C_As_p_2p, a_data_As_p_2p, b_data_As_p_2p, c_data_As_p_2p, lps_cc_data_As_p_2p, lps_cs_data_As_p_2p, lps_ss_data_As_p_2p, scale_factor_data_As_p_2p, window_c_data_As_p_2p, window_s_data_As_p_2p, mps_data_As_p_2p, z_at_chi_data_As_p_2p, cmbps_As_p_2p = data_import_func('data_As_p_2p')

cdef double lbs_As_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_p_2p, a_data_As_p_2p, b_data_As_p_2p, c_data_As_p_2p, scale_factor_data_As_p_2p, window_c_data_As_p_2p, window_s_data_As_p_2p, mps_data_As_p_2p, z_at_chi_data_As_p_2p)

cdef double lps_As_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_2p, lps_cs_data_As_p_2p, lps_ss_data_As_p_2p, cmbps_As_p_2p)
            

cdef double[:] cosm_par_As_m_2m
cdef double C_As_m_2m
cdef double[:, :] a_data_As_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_2m, C_As_m_2m, a_data_As_m_2m, b_data_As_m_2m, c_data_As_m_2m, lps_cc_data_As_m_2m, lps_cs_data_As_m_2m, lps_ss_data_As_m_2m, scale_factor_data_As_m_2m, window_c_data_As_m_2m, window_s_data_As_m_2m, mps_data_As_m_2m, z_at_chi_data_As_m_2m, cmbps_As_m_2m = data_import_func('data_As_m_2m')

cdef double lbs_As_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_2m, a_data_As_m_2m, b_data_As_m_2m, c_data_As_m_2m, scale_factor_data_As_m_2m, window_c_data_As_m_2m, window_s_data_As_m_2m, mps_data_As_m_2m, z_at_chi_data_As_m_2m)

cdef double lps_As_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_2m, lps_cs_data_As_m_2m, lps_ss_data_As_m_2m, cmbps_As_m_2m)
            

cdef double[:] cosm_par_As_m_1m
cdef double C_As_m_1m
cdef double[:, :] a_data_As_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_1m, C_As_m_1m, a_data_As_m_1m, b_data_As_m_1m, c_data_As_m_1m, lps_cc_data_As_m_1m, lps_cs_data_As_m_1m, lps_ss_data_As_m_1m, scale_factor_data_As_m_1m, window_c_data_As_m_1m, window_s_data_As_m_1m, mps_data_As_m_1m, z_at_chi_data_As_m_1m, cmbps_As_m_1m = data_import_func('data_As_m_1m')

cdef double lbs_As_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_1m, a_data_As_m_1m, b_data_As_m_1m, c_data_As_m_1m, scale_factor_data_As_m_1m, window_c_data_As_m_1m, window_s_data_As_m_1m, mps_data_As_m_1m, z_at_chi_data_As_m_1m)

cdef double lps_As_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_1m, lps_cs_data_As_m_1m, lps_ss_data_As_m_1m, cmbps_As_m_1m)
            

cdef double[:] cosm_par_As_m_0
cdef double C_As_m_0
cdef double[:, :] a_data_As_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_0, C_As_m_0, a_data_As_m_0, b_data_As_m_0, c_data_As_m_0, lps_cc_data_As_m_0, lps_cs_data_As_m_0, lps_ss_data_As_m_0, scale_factor_data_As_m_0, window_c_data_As_m_0, window_s_data_As_m_0, mps_data_As_m_0, z_at_chi_data_As_m_0, cmbps_As_m_0 = data_import_func('data_As_m_0')

cdef double lbs_As_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_0, a_data_As_m_0, b_data_As_m_0, c_data_As_m_0, scale_factor_data_As_m_0, window_c_data_As_m_0, window_s_data_As_m_0, mps_data_As_m_0, z_at_chi_data_As_m_0)

cdef double lps_As_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_0, lps_cs_data_As_m_0, lps_ss_data_As_m_0, cmbps_As_m_0)
            

cdef double[:] cosm_par_As_m_1p
cdef double C_As_m_1p
cdef double[:, :] a_data_As_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_1p, C_As_m_1p, a_data_As_m_1p, b_data_As_m_1p, c_data_As_m_1p, lps_cc_data_As_m_1p, lps_cs_data_As_m_1p, lps_ss_data_As_m_1p, scale_factor_data_As_m_1p, window_c_data_As_m_1p, window_s_data_As_m_1p, mps_data_As_m_1p, z_at_chi_data_As_m_1p, cmbps_As_m_1p = data_import_func('data_As_m_1p')

cdef double lbs_As_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_1p, a_data_As_m_1p, b_data_As_m_1p, c_data_As_m_1p, scale_factor_data_As_m_1p, window_c_data_As_m_1p, window_s_data_As_m_1p, mps_data_As_m_1p, z_at_chi_data_As_m_1p)

cdef double lps_As_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_1p, lps_cs_data_As_m_1p, lps_ss_data_As_m_1p, cmbps_As_m_1p)
            

cdef double[:] cosm_par_As_m_2p
cdef double C_As_m_2p
cdef double[:, :] a_data_As_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_As_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_As_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_As_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_As_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_As_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_As_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_As_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_As_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_As_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_As_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_As_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_As_m_2p, C_As_m_2p, a_data_As_m_2p, b_data_As_m_2p, c_data_As_m_2p, lps_cc_data_As_m_2p, lps_cs_data_As_m_2p, lps_ss_data_As_m_2p, scale_factor_data_As_m_2p, window_c_data_As_m_2p, window_s_data_As_m_2p, mps_data_As_m_2p, z_at_chi_data_As_m_2p, cmbps_As_m_2p = data_import_func('data_As_m_2p')

cdef double lbs_As_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_As_m_2p, a_data_As_m_2p, b_data_As_m_2p, c_data_As_m_2p, scale_factor_data_As_m_2p, window_c_data_As_m_2p, window_s_data_As_m_2p, mps_data_As_m_2p, z_at_chi_data_As_m_2p)

cdef double lps_As_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_2p, lps_cs_data_As_m_2p, lps_ss_data_As_m_2p, cmbps_As_m_2p)
            

cdef double[:] cosm_par_mnu_p_2m
cdef double C_mnu_p_2m
cdef double[:, :] a_data_mnu_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_2m, C_mnu_p_2m, a_data_mnu_p_2m, b_data_mnu_p_2m, c_data_mnu_p_2m, lps_cc_data_mnu_p_2m, lps_cs_data_mnu_p_2m, lps_ss_data_mnu_p_2m, scale_factor_data_mnu_p_2m, window_c_data_mnu_p_2m, window_s_data_mnu_p_2m, mps_data_mnu_p_2m, z_at_chi_data_mnu_p_2m, cmbps_mnu_p_2m = data_import_func('data_mnu_p_2m')

cdef double lbs_mnu_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_2m, a_data_mnu_p_2m, b_data_mnu_p_2m, c_data_mnu_p_2m, scale_factor_data_mnu_p_2m, window_c_data_mnu_p_2m, window_s_data_mnu_p_2m, mps_data_mnu_p_2m, z_at_chi_data_mnu_p_2m)

cdef double lps_mnu_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_2m, lps_cs_data_mnu_p_2m, lps_ss_data_mnu_p_2m, cmbps_mnu_p_2m)
            

cdef double[:] cosm_par_mnu_p_1m
cdef double C_mnu_p_1m
cdef double[:, :] a_data_mnu_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_1m, C_mnu_p_1m, a_data_mnu_p_1m, b_data_mnu_p_1m, c_data_mnu_p_1m, lps_cc_data_mnu_p_1m, lps_cs_data_mnu_p_1m, lps_ss_data_mnu_p_1m, scale_factor_data_mnu_p_1m, window_c_data_mnu_p_1m, window_s_data_mnu_p_1m, mps_data_mnu_p_1m, z_at_chi_data_mnu_p_1m, cmbps_mnu_p_1m = data_import_func('data_mnu_p_1m')

cdef double lbs_mnu_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_1m, a_data_mnu_p_1m, b_data_mnu_p_1m, c_data_mnu_p_1m, scale_factor_data_mnu_p_1m, window_c_data_mnu_p_1m, window_s_data_mnu_p_1m, mps_data_mnu_p_1m, z_at_chi_data_mnu_p_1m)

cdef double lps_mnu_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_1m, lps_cs_data_mnu_p_1m, lps_ss_data_mnu_p_1m, cmbps_mnu_p_1m)
            

cdef double[:] cosm_par_mnu_p_0
cdef double C_mnu_p_0
cdef double[:, :] a_data_mnu_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_0, C_mnu_p_0, a_data_mnu_p_0, b_data_mnu_p_0, c_data_mnu_p_0, lps_cc_data_mnu_p_0, lps_cs_data_mnu_p_0, lps_ss_data_mnu_p_0, scale_factor_data_mnu_p_0, window_c_data_mnu_p_0, window_s_data_mnu_p_0, mps_data_mnu_p_0, z_at_chi_data_mnu_p_0, cmbps_mnu_p_0 = data_import_func('data_mnu_p_0')

cdef double lbs_mnu_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_0, a_data_mnu_p_0, b_data_mnu_p_0, c_data_mnu_p_0, scale_factor_data_mnu_p_0, window_c_data_mnu_p_0, window_s_data_mnu_p_0, mps_data_mnu_p_0, z_at_chi_data_mnu_p_0)

cdef double lps_mnu_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_0, lps_cs_data_mnu_p_0, lps_ss_data_mnu_p_0, cmbps_mnu_p_0)
            

cdef double[:] cosm_par_mnu_p_1p
cdef double C_mnu_p_1p
cdef double[:, :] a_data_mnu_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_1p, C_mnu_p_1p, a_data_mnu_p_1p, b_data_mnu_p_1p, c_data_mnu_p_1p, lps_cc_data_mnu_p_1p, lps_cs_data_mnu_p_1p, lps_ss_data_mnu_p_1p, scale_factor_data_mnu_p_1p, window_c_data_mnu_p_1p, window_s_data_mnu_p_1p, mps_data_mnu_p_1p, z_at_chi_data_mnu_p_1p, cmbps_mnu_p_1p = data_import_func('data_mnu_p_1p')

cdef double lbs_mnu_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_1p, a_data_mnu_p_1p, b_data_mnu_p_1p, c_data_mnu_p_1p, scale_factor_data_mnu_p_1p, window_c_data_mnu_p_1p, window_s_data_mnu_p_1p, mps_data_mnu_p_1p, z_at_chi_data_mnu_p_1p)

cdef double lps_mnu_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_1p, lps_cs_data_mnu_p_1p, lps_ss_data_mnu_p_1p, cmbps_mnu_p_1p)
            

cdef double[:] cosm_par_mnu_p_2p
cdef double C_mnu_p_2p
cdef double[:, :] a_data_mnu_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_p_2p, C_mnu_p_2p, a_data_mnu_p_2p, b_data_mnu_p_2p, c_data_mnu_p_2p, lps_cc_data_mnu_p_2p, lps_cs_data_mnu_p_2p, lps_ss_data_mnu_p_2p, scale_factor_data_mnu_p_2p, window_c_data_mnu_p_2p, window_s_data_mnu_p_2p, mps_data_mnu_p_2p, z_at_chi_data_mnu_p_2p, cmbps_mnu_p_2p = data_import_func('data_mnu_p_2p')

cdef double lbs_mnu_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_p_2p, a_data_mnu_p_2p, b_data_mnu_p_2p, c_data_mnu_p_2p, scale_factor_data_mnu_p_2p, window_c_data_mnu_p_2p, window_s_data_mnu_p_2p, mps_data_mnu_p_2p, z_at_chi_data_mnu_p_2p)

cdef double lps_mnu_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_2p, lps_cs_data_mnu_p_2p, lps_ss_data_mnu_p_2p, cmbps_mnu_p_2p)
            

cdef double[:] cosm_par_mnu_m_2m
cdef double C_mnu_m_2m
cdef double[:, :] a_data_mnu_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_2m, C_mnu_m_2m, a_data_mnu_m_2m, b_data_mnu_m_2m, c_data_mnu_m_2m, lps_cc_data_mnu_m_2m, lps_cs_data_mnu_m_2m, lps_ss_data_mnu_m_2m, scale_factor_data_mnu_m_2m, window_c_data_mnu_m_2m, window_s_data_mnu_m_2m, mps_data_mnu_m_2m, z_at_chi_data_mnu_m_2m, cmbps_mnu_m_2m = data_import_func('data_mnu_m_2m')

cdef double lbs_mnu_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_2m, a_data_mnu_m_2m, b_data_mnu_m_2m, c_data_mnu_m_2m, scale_factor_data_mnu_m_2m, window_c_data_mnu_m_2m, window_s_data_mnu_m_2m, mps_data_mnu_m_2m, z_at_chi_data_mnu_m_2m)

cdef double lps_mnu_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_2m, lps_cs_data_mnu_m_2m, lps_ss_data_mnu_m_2m, cmbps_mnu_m_2m)
            

cdef double[:] cosm_par_mnu_m_1m
cdef double C_mnu_m_1m
cdef double[:, :] a_data_mnu_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_1m, C_mnu_m_1m, a_data_mnu_m_1m, b_data_mnu_m_1m, c_data_mnu_m_1m, lps_cc_data_mnu_m_1m, lps_cs_data_mnu_m_1m, lps_ss_data_mnu_m_1m, scale_factor_data_mnu_m_1m, window_c_data_mnu_m_1m, window_s_data_mnu_m_1m, mps_data_mnu_m_1m, z_at_chi_data_mnu_m_1m, cmbps_mnu_m_1m = data_import_func('data_mnu_m_1m')

cdef double lbs_mnu_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_1m, a_data_mnu_m_1m, b_data_mnu_m_1m, c_data_mnu_m_1m, scale_factor_data_mnu_m_1m, window_c_data_mnu_m_1m, window_s_data_mnu_m_1m, mps_data_mnu_m_1m, z_at_chi_data_mnu_m_1m)

cdef double lps_mnu_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_1m, lps_cs_data_mnu_m_1m, lps_ss_data_mnu_m_1m, cmbps_mnu_m_1m)
            

cdef double[:] cosm_par_mnu_m_0
cdef double C_mnu_m_0
cdef double[:, :] a_data_mnu_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_0, C_mnu_m_0, a_data_mnu_m_0, b_data_mnu_m_0, c_data_mnu_m_0, lps_cc_data_mnu_m_0, lps_cs_data_mnu_m_0, lps_ss_data_mnu_m_0, scale_factor_data_mnu_m_0, window_c_data_mnu_m_0, window_s_data_mnu_m_0, mps_data_mnu_m_0, z_at_chi_data_mnu_m_0, cmbps_mnu_m_0 = data_import_func('data_mnu_m_0')

cdef double lbs_mnu_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_0, a_data_mnu_m_0, b_data_mnu_m_0, c_data_mnu_m_0, scale_factor_data_mnu_m_0, window_c_data_mnu_m_0, window_s_data_mnu_m_0, mps_data_mnu_m_0, z_at_chi_data_mnu_m_0)

cdef double lps_mnu_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_0, lps_cs_data_mnu_m_0, lps_ss_data_mnu_m_0, cmbps_mnu_m_0)
            

cdef double[:] cosm_par_mnu_m_1p
cdef double C_mnu_m_1p
cdef double[:, :] a_data_mnu_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_1p, C_mnu_m_1p, a_data_mnu_m_1p, b_data_mnu_m_1p, c_data_mnu_m_1p, lps_cc_data_mnu_m_1p, lps_cs_data_mnu_m_1p, lps_ss_data_mnu_m_1p, scale_factor_data_mnu_m_1p, window_c_data_mnu_m_1p, window_s_data_mnu_m_1p, mps_data_mnu_m_1p, z_at_chi_data_mnu_m_1p, cmbps_mnu_m_1p = data_import_func('data_mnu_m_1p')

cdef double lbs_mnu_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_1p, a_data_mnu_m_1p, b_data_mnu_m_1p, c_data_mnu_m_1p, scale_factor_data_mnu_m_1p, window_c_data_mnu_m_1p, window_s_data_mnu_m_1p, mps_data_mnu_m_1p, z_at_chi_data_mnu_m_1p)

cdef double lps_mnu_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_1p, lps_cs_data_mnu_m_1p, lps_ss_data_mnu_m_1p, cmbps_mnu_m_1p)
            

cdef double[:] cosm_par_mnu_m_2p
cdef double C_mnu_m_2p
cdef double[:, :] a_data_mnu_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_mnu_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_mnu_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_mnu_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_mnu_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_mnu_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_mnu_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_mnu_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_mnu_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_mnu_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_mnu_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_mnu_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_mnu_m_2p, C_mnu_m_2p, a_data_mnu_m_2p, b_data_mnu_m_2p, c_data_mnu_m_2p, lps_cc_data_mnu_m_2p, lps_cs_data_mnu_m_2p, lps_ss_data_mnu_m_2p, scale_factor_data_mnu_m_2p, window_c_data_mnu_m_2p, window_s_data_mnu_m_2p, mps_data_mnu_m_2p, z_at_chi_data_mnu_m_2p, cmbps_mnu_m_2p = data_import_func('data_mnu_m_2p')

cdef double lbs_mnu_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_mnu_m_2p, a_data_mnu_m_2p, b_data_mnu_m_2p, c_data_mnu_m_2p, scale_factor_data_mnu_m_2p, window_c_data_mnu_m_2p, window_s_data_mnu_m_2p, mps_data_mnu_m_2p, z_at_chi_data_mnu_m_2p)

cdef double lps_mnu_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_2p, lps_cs_data_mnu_m_2p, lps_ss_data_mnu_m_2p, cmbps_mnu_m_2p)
            

cdef double[:] cosm_par_w0_p_2m
cdef double C_w0_p_2m
cdef double[:, :] a_data_w0_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_2m, C_w0_p_2m, a_data_w0_p_2m, b_data_w0_p_2m, c_data_w0_p_2m, lps_cc_data_w0_p_2m, lps_cs_data_w0_p_2m, lps_ss_data_w0_p_2m, scale_factor_data_w0_p_2m, window_c_data_w0_p_2m, window_s_data_w0_p_2m, mps_data_w0_p_2m, z_at_chi_data_w0_p_2m, cmbps_w0_p_2m = data_import_func('data_w0_p_2m')

cdef double lbs_w0_p_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_2m, a_data_w0_p_2m, b_data_w0_p_2m, c_data_w0_p_2m, scale_factor_data_w0_p_2m, window_c_data_w0_p_2m, window_s_data_w0_p_2m, mps_data_w0_p_2m, z_at_chi_data_w0_p_2m)

cdef double lps_w0_p_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_2m, lps_cs_data_w0_p_2m, lps_ss_data_w0_p_2m, cmbps_w0_p_2m)
            

cdef double[:] cosm_par_w0_p_1m
cdef double C_w0_p_1m
cdef double[:, :] a_data_w0_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_1m, C_w0_p_1m, a_data_w0_p_1m, b_data_w0_p_1m, c_data_w0_p_1m, lps_cc_data_w0_p_1m, lps_cs_data_w0_p_1m, lps_ss_data_w0_p_1m, scale_factor_data_w0_p_1m, window_c_data_w0_p_1m, window_s_data_w0_p_1m, mps_data_w0_p_1m, z_at_chi_data_w0_p_1m, cmbps_w0_p_1m = data_import_func('data_w0_p_1m')

cdef double lbs_w0_p_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_1m, a_data_w0_p_1m, b_data_w0_p_1m, c_data_w0_p_1m, scale_factor_data_w0_p_1m, window_c_data_w0_p_1m, window_s_data_w0_p_1m, mps_data_w0_p_1m, z_at_chi_data_w0_p_1m)

cdef double lps_w0_p_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_1m, lps_cs_data_w0_p_1m, lps_ss_data_w0_p_1m, cmbps_w0_p_1m)
            

cdef double[:] cosm_par_w0_p_0
cdef double C_w0_p_0
cdef double[:, :] a_data_w0_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_0, C_w0_p_0, a_data_w0_p_0, b_data_w0_p_0, c_data_w0_p_0, lps_cc_data_w0_p_0, lps_cs_data_w0_p_0, lps_ss_data_w0_p_0, scale_factor_data_w0_p_0, window_c_data_w0_p_0, window_s_data_w0_p_0, mps_data_w0_p_0, z_at_chi_data_w0_p_0, cmbps_w0_p_0 = data_import_func('data_w0_p_0')

cdef double lbs_w0_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_0, a_data_w0_p_0, b_data_w0_p_0, c_data_w0_p_0, scale_factor_data_w0_p_0, window_c_data_w0_p_0, window_s_data_w0_p_0, mps_data_w0_p_0, z_at_chi_data_w0_p_0)

cdef double lps_w0_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_0, lps_cs_data_w0_p_0, lps_ss_data_w0_p_0, cmbps_w0_p_0)
            

cdef double[:] cosm_par_w0_p_1p
cdef double C_w0_p_1p
cdef double[:, :] a_data_w0_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_1p, C_w0_p_1p, a_data_w0_p_1p, b_data_w0_p_1p, c_data_w0_p_1p, lps_cc_data_w0_p_1p, lps_cs_data_w0_p_1p, lps_ss_data_w0_p_1p, scale_factor_data_w0_p_1p, window_c_data_w0_p_1p, window_s_data_w0_p_1p, mps_data_w0_p_1p, z_at_chi_data_w0_p_1p, cmbps_w0_p_1p = data_import_func('data_w0_p_1p')

cdef double lbs_w0_p_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_1p, a_data_w0_p_1p, b_data_w0_p_1p, c_data_w0_p_1p, scale_factor_data_w0_p_1p, window_c_data_w0_p_1p, window_s_data_w0_p_1p, mps_data_w0_p_1p, z_at_chi_data_w0_p_1p)

cdef double lps_w0_p_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_1p, lps_cs_data_w0_p_1p, lps_ss_data_w0_p_1p, cmbps_w0_p_1p)
            

cdef double[:] cosm_par_w0_p_2p
cdef double C_w0_p_2p
cdef double[:, :] a_data_w0_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_p_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_p_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_p_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_p_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_p_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_p_2p, C_w0_p_2p, a_data_w0_p_2p, b_data_w0_p_2p, c_data_w0_p_2p, lps_cc_data_w0_p_2p, lps_cs_data_w0_p_2p, lps_ss_data_w0_p_2p, scale_factor_data_w0_p_2p, window_c_data_w0_p_2p, window_s_data_w0_p_2p, mps_data_w0_p_2p, z_at_chi_data_w0_p_2p, cmbps_w0_p_2p = data_import_func('data_w0_p_2p')

cdef double lbs_w0_p_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_p_2p, a_data_w0_p_2p, b_data_w0_p_2p, c_data_w0_p_2p, scale_factor_data_w0_p_2p, window_c_data_w0_p_2p, window_s_data_w0_p_2p, mps_data_w0_p_2p, z_at_chi_data_w0_p_2p)

cdef double lps_w0_p_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_2p, lps_cs_data_w0_p_2p, lps_ss_data_w0_p_2p, cmbps_w0_p_2p)
            

cdef double[:] cosm_par_w0_m_2m
cdef double C_w0_m_2m
cdef double[:, :] a_data_w0_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_2m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_2m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_2m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_2m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_2m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_2m, C_w0_m_2m, a_data_w0_m_2m, b_data_w0_m_2m, c_data_w0_m_2m, lps_cc_data_w0_m_2m, lps_cs_data_w0_m_2m, lps_ss_data_w0_m_2m, scale_factor_data_w0_m_2m, window_c_data_w0_m_2m, window_s_data_w0_m_2m, mps_data_w0_m_2m, z_at_chi_data_w0_m_2m, cmbps_w0_m_2m = data_import_func('data_w0_m_2m')

cdef double lbs_w0_m_2m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_2m, a_data_w0_m_2m, b_data_w0_m_2m, c_data_w0_m_2m, scale_factor_data_w0_m_2m, window_c_data_w0_m_2m, window_s_data_w0_m_2m, mps_data_w0_m_2m, z_at_chi_data_w0_m_2m)

cdef double lps_w0_m_2m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_2m, lps_cs_data_w0_m_2m, lps_ss_data_w0_m_2m, cmbps_w0_m_2m)
            

cdef double[:] cosm_par_w0_m_1m
cdef double C_w0_m_1m
cdef double[:, :] a_data_w0_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_1m = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_1m = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_1m = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_1m = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_1m = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_1m, C_w0_m_1m, a_data_w0_m_1m, b_data_w0_m_1m, c_data_w0_m_1m, lps_cc_data_w0_m_1m, lps_cs_data_w0_m_1m, lps_ss_data_w0_m_1m, scale_factor_data_w0_m_1m, window_c_data_w0_m_1m, window_s_data_w0_m_1m, mps_data_w0_m_1m, z_at_chi_data_w0_m_1m, cmbps_w0_m_1m = data_import_func('data_w0_m_1m')

cdef double lbs_w0_m_1m(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_1m, a_data_w0_m_1m, b_data_w0_m_1m, c_data_w0_m_1m, scale_factor_data_w0_m_1m, window_c_data_w0_m_1m, window_s_data_w0_m_1m, mps_data_w0_m_1m, z_at_chi_data_w0_m_1m)

cdef double lps_w0_m_1m(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_1m, lps_cs_data_w0_m_1m, lps_ss_data_w0_m_1m, cmbps_w0_m_1m)
            

cdef double[:] cosm_par_w0_m_0
cdef double C_w0_m_0
cdef double[:, :] a_data_w0_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_0, C_w0_m_0, a_data_w0_m_0, b_data_w0_m_0, c_data_w0_m_0, lps_cc_data_w0_m_0, lps_cs_data_w0_m_0, lps_ss_data_w0_m_0, scale_factor_data_w0_m_0, window_c_data_w0_m_0, window_s_data_w0_m_0, mps_data_w0_m_0, z_at_chi_data_w0_m_0, cmbps_w0_m_0 = data_import_func('data_w0_m_0')

cdef double lbs_w0_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_0, a_data_w0_m_0, b_data_w0_m_0, c_data_w0_m_0, scale_factor_data_w0_m_0, window_c_data_w0_m_0, window_s_data_w0_m_0, mps_data_w0_m_0, z_at_chi_data_w0_m_0)

cdef double lps_w0_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_0, lps_cs_data_w0_m_0, lps_ss_data_w0_m_0, cmbps_w0_m_0)
            

cdef double[:] cosm_par_w0_m_1p
cdef double C_w0_m_1p
cdef double[:, :] a_data_w0_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_1p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_1p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_1p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_1p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_1p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_1p, C_w0_m_1p, a_data_w0_m_1p, b_data_w0_m_1p, c_data_w0_m_1p, lps_cc_data_w0_m_1p, lps_cs_data_w0_m_1p, lps_ss_data_w0_m_1p, scale_factor_data_w0_m_1p, window_c_data_w0_m_1p, window_s_data_w0_m_1p, mps_data_w0_m_1p, z_at_chi_data_w0_m_1p, cmbps_w0_m_1p = data_import_func('data_w0_m_1p')

cdef double lbs_w0_m_1p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_1p, a_data_w0_m_1p, b_data_w0_m_1p, c_data_w0_m_1p, scale_factor_data_w0_m_1p, window_c_data_w0_m_1p, window_s_data_w0_m_1p, mps_data_w0_m_1p, z_at_chi_data_w0_m_1p)

cdef double lps_w0_m_1p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_1p, lps_cs_data_w0_m_1p, lps_ss_data_w0_m_1p, cmbps_w0_m_1p)
            

cdef double[:] cosm_par_w0_m_2p
cdef double C_w0_m_2p
cdef double[:, :] a_data_w0_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_w0_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_w0_m_2p = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_w0_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_w0_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_w0_m_2p = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_w0_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_w0_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_w0_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_w0_m_2p = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_w0_m_2p = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_w0_m_2p = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_w0_m_2p, C_w0_m_2p, a_data_w0_m_2p, b_data_w0_m_2p, c_data_w0_m_2p, lps_cc_data_w0_m_2p, lps_cs_data_w0_m_2p, lps_ss_data_w0_m_2p, scale_factor_data_w0_m_2p, window_c_data_w0_m_2p, window_s_data_w0_m_2p, mps_data_w0_m_2p, z_at_chi_data_w0_m_2p, cmbps_w0_m_2p = data_import_func('data_w0_m_2p')

cdef double lbs_w0_m_2p(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, C_w0_m_2p, a_data_w0_m_2p, b_data_w0_m_2p, c_data_w0_m_2p, scale_factor_data_w0_m_2p, window_c_data_w0_m_2p, window_s_data_w0_m_2p, mps_data_w0_m_2p, z_at_chi_data_w0_m_2p)

cdef double lps_w0_m_2p(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_2p, lps_cs_data_w0_m_2p, lps_ss_data_w0_m_2p, cmbps_w0_m_2p)
            
