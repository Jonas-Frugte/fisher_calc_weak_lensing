

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

cdef double lbs_H_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_H_p_0, a_data_H_p_0, b_data_H_p_0, c_data_H_p_0, scale_factor_data_H_p_0, window_c_data_H_p_0, window_s_data_H_p_0, mps_data_H_p_0, z_at_chi_data_H_p_0)

cdef double lps_H_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_p_0, lps_cs_data_H_p_0, lps_ss_data_H_p_0, cmbps_H_p_0)
            

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

cdef double lbs_H_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_H_m_0, a_data_H_m_0, b_data_H_m_0, c_data_H_m_0, scale_factor_data_H_m_0, window_c_data_H_m_0, window_s_data_H_m_0, mps_data_H_m_0, z_at_chi_data_H_m_0)

cdef double lps_H_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_H_m_0, lps_cs_data_H_m_0, lps_ss_data_H_m_0, cmbps_H_m_0)
            

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

cdef double lbs_ombh2_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ombh2_p_0, a_data_ombh2_p_0, b_data_ombh2_p_0, c_data_ombh2_p_0, scale_factor_data_ombh2_p_0, window_c_data_ombh2_p_0, window_s_data_ombh2_p_0, mps_data_ombh2_p_0, z_at_chi_data_ombh2_p_0)

cdef double lps_ombh2_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_p_0, lps_cs_data_ombh2_p_0, lps_ss_data_ombh2_p_0, cmbps_ombh2_p_0)
            

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

cdef double lbs_ombh2_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ombh2_m_0, a_data_ombh2_m_0, b_data_ombh2_m_0, c_data_ombh2_m_0, scale_factor_data_ombh2_m_0, window_c_data_ombh2_m_0, window_s_data_ombh2_m_0, mps_data_ombh2_m_0, z_at_chi_data_ombh2_m_0)

cdef double lps_ombh2_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ombh2_m_0, lps_cs_data_ombh2_m_0, lps_ss_data_ombh2_m_0, cmbps_ombh2_m_0)
            

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

cdef double lbs_omch2_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_omch2_p_0, a_data_omch2_p_0, b_data_omch2_p_0, c_data_omch2_p_0, scale_factor_data_omch2_p_0, window_c_data_omch2_p_0, window_s_data_omch2_p_0, mps_data_omch2_p_0, z_at_chi_data_omch2_p_0)

cdef double lps_omch2_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_p_0, lps_cs_data_omch2_p_0, lps_ss_data_omch2_p_0, cmbps_omch2_p_0)
            

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

cdef double lbs_omch2_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_omch2_m_0, a_data_omch2_m_0, b_data_omch2_m_0, c_data_omch2_m_0, scale_factor_data_omch2_m_0, window_c_data_omch2_m_0, window_s_data_omch2_m_0, mps_data_omch2_m_0, z_at_chi_data_omch2_m_0)

cdef double lps_omch2_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_omch2_m_0, lps_cs_data_omch2_m_0, lps_ss_data_omch2_m_0, cmbps_omch2_m_0)
            

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

cdef double lbs_ns_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ns_p_0, a_data_ns_p_0, b_data_ns_p_0, c_data_ns_p_0, scale_factor_data_ns_p_0, window_c_data_ns_p_0, window_s_data_ns_p_0, mps_data_ns_p_0, z_at_chi_data_ns_p_0)

cdef double lps_ns_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_p_0, lps_cs_data_ns_p_0, lps_ss_data_ns_p_0, cmbps_ns_p_0)
            

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

cdef double lbs_ns_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_ns_m_0, a_data_ns_m_0, b_data_ns_m_0, c_data_ns_m_0, scale_factor_data_ns_m_0, window_c_data_ns_m_0, window_s_data_ns_m_0, mps_data_ns_m_0, z_at_chi_data_ns_m_0)

cdef double lps_ns_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_ns_m_0, lps_cs_data_ns_m_0, lps_ss_data_ns_m_0, cmbps_ns_m_0)
            

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

cdef double lbs_As_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_As_p_0, a_data_As_p_0, b_data_As_p_0, c_data_As_p_0, scale_factor_data_As_p_0, window_c_data_As_p_0, window_s_data_As_p_0, mps_data_As_p_0, z_at_chi_data_As_p_0)

cdef double lps_As_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_p_0, lps_cs_data_As_p_0, lps_ss_data_As_p_0, cmbps_As_p_0)
            

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

cdef double lbs_As_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_As_m_0, a_data_As_m_0, b_data_As_m_0, c_data_As_m_0, scale_factor_data_As_m_0, window_c_data_As_m_0, window_s_data_As_m_0, mps_data_As_m_0, z_at_chi_data_As_m_0)

cdef double lps_As_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_As_m_0, lps_cs_data_As_m_0, lps_ss_data_As_m_0, cmbps_As_m_0)
            

cdef double[:] cosm_par_tau_p_0
cdef double C_tau_p_0
cdef double[:, :] a_data_tau_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_tau_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_tau_p_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_tau_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_tau_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_tau_p_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_tau_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_tau_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_tau_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_tau_p_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_tau_p_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_tau_p_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_tau_p_0, C_tau_p_0, a_data_tau_p_0, b_data_tau_p_0, c_data_tau_p_0, lps_cc_data_tau_p_0, lps_cs_data_tau_p_0, lps_ss_data_tau_p_0, scale_factor_data_tau_p_0, window_c_data_tau_p_0, window_s_data_tau_p_0, mps_data_tau_p_0, z_at_chi_data_tau_p_0, cmbps_tau_p_0 = data_import_func('data_tau_p_0')

cdef double lbs_tau_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_tau_p_0, a_data_tau_p_0, b_data_tau_p_0, c_data_tau_p_0, scale_factor_data_tau_p_0, window_c_data_tau_p_0, window_s_data_tau_p_0, mps_data_tau_p_0, z_at_chi_data_tau_p_0)

cdef double lps_tau_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_tau_p_0, lps_cs_data_tau_p_0, lps_ss_data_tau_p_0, cmbps_tau_p_0)
            

cdef double[:] cosm_par_tau_m_0
cdef double C_tau_m_0
cdef double[:, :] a_data_tau_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] b_data_tau_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:, :] c_data_tau_m_0 = np.zeros((z_num, k_num), dtype=np.float64)
cdef double[:] lps_cc_data_tau_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_cs_data_tau_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] lps_ss_data_tau_m_0 = np.zeros(k_num, dtype=np.float64)
cdef double[:] scale_factor_data_tau_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_c_data_tau_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:] window_s_data_tau_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] mps_data_tau_m_0 = np.zeros((k_num_fine, z_num_fine), dtype=np.float64)
cdef double[:] z_at_chi_data_tau_m_0 = np.zeros(chi_num, dtype=np.float64)
cdef double[:, :] cmbps_tau_m_0 = np.zeros((lmax_cmbps+1, 5), dtype=np.float64)

cosm_par_tau_m_0, C_tau_m_0, a_data_tau_m_0, b_data_tau_m_0, c_data_tau_m_0, lps_cc_data_tau_m_0, lps_cs_data_tau_m_0, lps_ss_data_tau_m_0, scale_factor_data_tau_m_0, window_c_data_tau_m_0, window_s_data_tau_m_0, mps_data_tau_m_0, z_at_chi_data_tau_m_0, cmbps_tau_m_0 = data_import_func('data_tau_m_0')

cdef double lbs_tau_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_tau_m_0, a_data_tau_m_0, b_data_tau_m_0, c_data_tau_m_0, scale_factor_data_tau_m_0, window_c_data_tau_m_0, window_s_data_tau_m_0, mps_data_tau_m_0, z_at_chi_data_tau_m_0)

cdef double lps_tau_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_tau_m_0, lps_cs_data_tau_m_0, lps_ss_data_tau_m_0, cmbps_tau_m_0)
            

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

cdef double lbs_mnu_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_mnu_p_0, a_data_mnu_p_0, b_data_mnu_p_0, c_data_mnu_p_0, scale_factor_data_mnu_p_0, window_c_data_mnu_p_0, window_s_data_mnu_p_0, mps_data_mnu_p_0, z_at_chi_data_mnu_p_0)

cdef double lps_mnu_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_p_0, lps_cs_data_mnu_p_0, lps_ss_data_mnu_p_0, cmbps_mnu_p_0)
            

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

cdef double lbs_mnu_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_mnu_m_0, a_data_mnu_m_0, b_data_mnu_m_0, c_data_mnu_m_0, scale_factor_data_mnu_m_0, window_c_data_mnu_m_0, window_s_data_mnu_m_0, mps_data_mnu_m_0, z_at_chi_data_mnu_m_0)

cdef double lps_mnu_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_mnu_m_0, lps_cs_data_mnu_m_0, lps_ss_data_mnu_m_0, cmbps_mnu_m_0)
            

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

cdef double lbs_w0_p_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_w0_p_0, a_data_w0_p_0, b_data_w0_p_0, c_data_w0_p_0, scale_factor_data_w0_p_0, window_c_data_w0_p_0, window_s_data_w0_p_0, mps_data_w0_p_0, z_at_chi_data_w0_p_0)

cdef double lps_w0_p_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_p_0, lps_cs_data_w0_p_0, lps_ss_data_w0_p_0, cmbps_w0_p_0)
            

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

cdef double lbs_w0_m_0(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil:
    return lbs(k1, k2, k3, type1, type2, type3, num_samples, pb_correction, C_w0_m_0, a_data_w0_m_0, b_data_w0_m_0, c_data_w0_m_0, scale_factor_data_w0_m_0, window_c_data_w0_m_0, window_s_data_w0_m_0, mps_data_w0_m_0, z_at_chi_data_w0_m_0)

cdef double lps_w0_m_0(int l, char* type1, char* type2) noexcept nogil:
    return lps(l, type1, type2, lps_cc_data_w0_m_0, lps_cs_data_w0_m_0, lps_ss_data_w0_m_0, cmbps_w0_m_0)
            
