# header for cython_interpolation.pyx

from libc.math cimport log10

# Declare cpdef functions
cpdef double linear_interp(double x, double xmin, double xmax, int n_points, double[:] ys) noexcept nogil
cpdef double logspace_linear_interp(double x, double xmin, double xmax, int n_points, double[:] ys) noexcept nogil
cpdef double linear_interp2d(double x, double y, double xmin, double xmax, int n_points_x, double ymin, double ymax, int n_points_y, double[:, :] zs) noexcept nogil
cpdef double logspace_linear_interp2d(double x, double y, double xmin, double xmax, int n_points_x, double ymin, double ymax, int n_points_y, double[:, :] zs) noexcept nogil