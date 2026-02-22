cdef int index_3d(int X, int Y, int Z, int n_tracers) noexcept nogil

cdef int index_2d(int X, int Y, int n_tracers) noexcept nogil

cdef int solve_3(
    double* CXp, # n_tracer^2
    double* CYp, # n_tracer^2
    double* CZp, # n_tracer^2
    double* bispec_vec, # n_tracer^3
    double* outn3, # n_tracer^3
    double* vecss, # n_tracer^3
    double* b, # n_tracer
    double* x, # n_tracer
    double* x2, # n_tracer^2
    double* vecs_solve2, # n_tracer^2
    double* b_solve2, # n_tracer
    double* x_solve2, # n_tracer
    double* vec_solve1, # n_tracer
    double* A_solve, # n_tracer^2
    double* b_solve, # n_tracer
    double* y_solve, # n_tracer
    int n_tracers
    ) noexcept nogil
