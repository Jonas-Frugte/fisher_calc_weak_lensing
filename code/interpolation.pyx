from libc.math cimport log10, isnan, isinf
from libc.stdio cimport printf

cpdef double linear_interp(double x, double xmin, double xmax, int n_points, double[:] ys) noexcept nogil:
    # xs should be of the form np.linspace(xmin, xmax, n_points)
    # xs[i] = xmin + i * (xmax - xmin) / (n_points - 1)
    if x > xmax or x < xmin:
        return 0 # fill value outside of interp range

    cdef double i_f = (x - xmin) * (n_points - 1) / (xmax - xmin)
    cdef int i = int(i_f)
    cdef double dx = x - (xmin + i * (xmax - xmin) / (n_points - 1))
    cdef double dydx = (ys[i + 1] - ys[i]) / ((xmax - xmin) / (n_points - 1))


    if isnan(ys[i] + dx * dydx) or isinf(ys[i] + dx * dydx):
        printf('nonono')
    return ys[i] + dx * dydx

cpdef double logspace_linear_interp(double x, double xmin, double xmax, int n_points, double[:] ys) noexcept nogil:
    # xs should be of the form np.logspace(np.log10(xmin), np.log10(xmax), n_points)
    # xs[i] = 10 ** (log10(xmin) + i * (log10(xmax) - log10(xmin)) / (n_points - 1))
    if x > xmax or x < xmin:
        return 0 # fill value outside of interp range

    cdef double xmin_log = log10(xmin)
    cdef double xmax_log = log10(xmax)

    cdef double i_f = (log10(x) - xmin_log) * (n_points - 1) / (xmax_log - xmin_log)
    cdef int i = int(i_f)
    cdef double dx = x - 10 ** (xmin_log + i * (xmax_log - xmin_log) / (n_points - 1))
    cdef double dydx = (ys[i + 1] - ys[i]) / (10 ** (xmin_log + (i + 1) * (xmax_log - xmin_log) / (n_points - 1)) - 10 ** (xmin_log + i * (xmax_log - xmin_log) / (n_points - 1)))

    # if i >= ys.shape[0]:
    #     printf("oh noooo (%d and %d)", i, ys.shape[0])

    if isnan(ys[i] + dx * dydx) or isinf(ys[i] + dx * dydx):
        printf('NaN value in logspace_linear_interp')

    return ys[i] + dx * dydx

cdef double two_pt_interp(double x, double xmin, double xmax, double ymin, double ymax) noexcept nogil:
    cdef double dx = x - xmin
    cdef double dydx = (ymax - ymin) / (xmax - xmin)
    return ymin + dx * dydx

cdef inline int linear_index(double x, double xmin, double xmax, int n_points) noexcept nogil:
    i_f = (x - xmin) * (n_points - 1) / (xmax - xmin)
    return int(i_f)

cdef inline int logspace_index(double x, double xmin, double xmax, int n_points) noexcept nogil:
    i_f = (log10(x) - log10(xmin)) * (n_points - 1) / (log10(xmax) - log10(xmin))
    return int(i_f)

cdef inline double linear_i_to_x(int i, double xmin, double xmax, int n_points) noexcept nogil:
    return xmin + i * (xmax - xmin) / (n_points - 1)

cdef inline double logspace_i_to_x(int i, double xmin, double xmax, int n_points) noexcept nogil:
    return 10 ** (log10(xmin) + i * (log10(xmax) - log10(xmin)) / (n_points - 1))

cpdef double linear_interp2d(double x, double y, double xmin, double xmax, int n_points_x, double ymin, double ymax, int n_points_y, double[:, :] zs) noexcept nogil:
    # associates x[i], y[j] to z[i, j]
    # x, y are of the form np.linspace(xmin, xmax, n_points_x), ...
    if x < xmin or x > xmax or y < ymin or y > ymax:
        return 0
    cdef int i_x = linear_index(x, xmin, xmax, n_points_x)
    cdef int i_y = linear_index(y, ymin, ymax, n_points_y)
    return two_pt_interp(
        y,
        linear_i_to_x(i_y, ymin, ymax, n_points_y),
        linear_i_to_x(i_y + 1, ymin, ymax, n_points_y),
        two_pt_interp(x, linear_i_to_x(i_x, xmin, xmax, n_points_x), linear_i_to_x(i_x + 1, xmin, xmax, n_points_x), zs[i_x, i_y], zs[i_x + 1, i_y]),
        two_pt_interp(x, linear_i_to_x(i_x, xmin, xmax, n_points_x), linear_i_to_x(i_x + 1, xmin, xmax, n_points_x), zs[i_x, i_y + 1], zs[i_x + 1, i_y + 1])      
    )

cpdef double logspace_linear_interp2d(double x, double y, double xmin, double xmax, int n_points_x, double ymin, double ymax, int n_points_y, double[:, :] zs) noexcept nogil:
    # associates x[i], y[j] to z[i, j]
    # x of the form np.logspace(log10(xmin), log10(xmax), n_points_x), np.linspace(ymax, ymin, n_points_y)
  
   
    if x < xmin or x > xmax or y < ymin or y > ymax:
        return 0
    cdef int i_x = logspace_index(x, xmin, xmax, n_points_x)
    cdef int i_y = linear_index(y, ymin, ymax, n_points_y)
    
    # if i_x + 1 >= zs.shape[0] or i_y + 1 >= zs.shape[1] or i_x < 0 or i_y < 0:
    #     printf("just shat myself")
    #     printf("shape of arrays: (%d, %d)", zs.shape[0], zs.shape[1])
    #     printf("indexes: %d and %d", i_x + 1, i_y + 1)

    cdef float result = two_pt_interp(
        y,
        linear_i_to_x(i_y, ymin, ymax, n_points_y),
        linear_i_to_x(i_y + 1, ymin, ymax, n_points_y),
        two_pt_interp(x, logspace_i_to_x(i_x, xmin, xmax, n_points_x), logspace_i_to_x(i_x + 1, xmin, xmax, n_points_x), zs[i_x, i_y], zs[i_x + 1, i_y]),
        two_pt_interp(x, logspace_i_to_x(i_x, xmin, xmax, n_points_x), logspace_i_to_x(i_x + 1, xmin, xmax, n_points_x), zs[i_x, i_y + 1], zs[i_x + 1, i_y + 1])      
    )

    if isnan(result) or isinf(result):
        printf('just shat myself')
    else:
        return result
