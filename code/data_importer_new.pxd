# cdef double lensing_bi_spectrum_f(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil
cdef double lps_f(int l, char* type1, char* type2) noexcept nogil
cdef double lps_f_obs(int l, char* type1, char* type2) noexcept nogil
cdef double lbs_f(int l_1, int l_2, int l_3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction) noexcept nogil

cdef double lbs_der(int k1, int k2, int k3, char* type1, char* type2, char* type3, int num_samples, bint pb_correction, char* par, int delta_delta) noexcept nogil
cdef double lps_der(int k, char* type1, char* type2, char* par, int delta_delta) noexcept nogil

# cpdef double a_test(double l, double z)

# cpdef double b_test(double l, double z)

# cpdef double c_test(double l, double z)

# cpdef double a_test(double l, double z)
# cpdef double b_test(double l, double z)
# cpdef double c_test(double l, double z)

#cdef double get_k_max()

# cdef double lensing_bi_spectrum_full_sky(int l_1, int l_2, int l_3, char* type1, char* type2, char* type3, int num_samples) noexcept nogil
# cdef double lbs_der_full_sky(int l_1, int l_2, int l_3, char* type1, char* type2, char* type3, int num_samples, char* par) noexcept nogil