cdef inline int index_3d_5(int X, int Y, int Z) noexcept nogil

cdef inline int index_2d_5(int X, int Y) noexcept nogil

cdef int soe_solver_5(double* matrix, double* vec, double* sol) noexcept nogil

cdef int solve_1(double* CXp, int Yp, int Zp, double* bispec_vec, double* out5) noexcept nogil

cdef int solve_2(double* CXp, double* CYp, int Zp, double* bispec_vec, double* out25) noexcept nogil

cdef int solve_3(double* CXp, double* CYp, double* CZp, double* bispec_vec, double* out125) noexcept nogil
