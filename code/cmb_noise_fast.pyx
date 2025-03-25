# follows quadratic estimator of https://arxiv.org/pdf/astro-ph/0111606

import numpy as np
from scipy.integrate import quad, dblquad, simpson
from scipy.interpolate import interp1d
import multiprocessing
cimport interpolation as interp
from libc.math cimport sqrt, atan2, cos, sin, log, exp, abs

# Get C_l for TT, EE, BB, TE in DIMENSIONLESS UNITS!! up to some lmax
cdef double [:, :] cls_ul = np.loadtxt('cmb_ps_unlensed.txt')
cdef double [:, :] cls_l = np.loadtxt('cmb_ps_lensed.txt')

# cls_l[:,0] = TT, cls_l[:,1] = EE, cls_l[:,2] = BB, cls_l[:,3] = TE
cdef double lmin = 0
cdef double lmax = 10000
cdef int lnum = 10001

cdef double fTT_ul(double l) noexcept nogil:
    return interp.linear_interp(l, lmin, lmax, lnum, cls_ul[:, 0])
cdef double fEE_ul(double l) noexcept nogil:
    return interp.linear_interp(l, lmin, lmax, lnum, cls_ul[:, 1])
cdef double fBB_ul(double l) noexcept nogil:
    return interp.linear_interp(l, lmin, lmax, lnum, cls_ul[:, 2])
cdef double fTE_ul(double l) noexcept nogil:
    return interp.linear_interp(l, lmin, lmax, lnum, cls_ul[:, 3])

cdef double fTT_l(double l) noexcept nogil:
    return interp.linear_interp(l, lmin, lmax, lnum, cls_l[:, 0])
cdef double fEE_l(double l) noexcept nogil:
    return interp.linear_interp(l, lmin, lmax, lnum, cls_l[:, 1])
cdef double fBB_l(double l) noexcept nogil:
    return interp.linear_interp(l, lmin, lmax, lnum, cls_l[:, 2])
cdef double fTE_l(double l) noexcept nogil:
    return interp.linear_interp(l, lmin, lmax, lnum, cls_l[:, 3])

# for types:
# 1: T
# 2: E
# 3: B

cdef double cmbps_ul(double l, int type1, int type2) noexcept nogil:
    # Returns C_l in μK^2 for the given combination (type1, type2).
    # type1, type2 can be 'T', 'E', or 'B'.
    if type1 == 1 and type2 == 1:
        return fTT_ul(l)
    elif type1 == 2 and type2 == 2:
        return fEE_ul(l)
    elif type1 == 3 and type2 == 3:
        return fBB_ul(l)
    elif (type1 == 1 and type2 == 2) or (type1 == 2 and type2 == 1):
        return fTE_ul(l)
    else:
        # For non-standard combos like TB or EB (normally zero in standard ΛCDM),
        # just return zero unless you have them computed.
        return 0.0

cdef double cmbps_l(double l, int type1, int type2) noexcept nogil:
    # Returns C_l in μK^2 for the given combination (type1, type2).
    # type1, type2 can be 'T', 'E', or 'B'.
    if type1 == 1 and type2 == 1:
        return fTT_l(l)
    elif type1 == 2 and type2 == 2:
        return fEE_l(l)
    elif type1 == 3 and type2 == 3:
        return fBB_l(l)
    elif (type1 == 1 and type2 == 2) or (type1 == 2 and type2 == 1):
        return fTE_l(l)
    else:
        # For non-standard combos like TB or EB (normally zero in standard ΛCDM),
        # just return zero unless you have them computed.
        return 0.0

cdef double arcmintorad = np.pi / 10800 # convert arcminutes to radians

cdef double cmbps_noise(double l, int type1, int type2, double sigma, double Delta_T, double Delta_P) noexcept nogil:
    cdef double noise = 0

    # units to input:
    # sigma: arcmin
    # Delta_T, Delta_P: microKelvin arcmin

    # constant time!
    cdef double Tcmb = 2.728e6 # in micro kelvin
    sigma *= arcmintorad # in radians
    Delta_T *= arcmintorad
    Delta_P *= arcmintorad
    
    if type1 == 1 and type2 == 1:
        noise = (Delta_T / Tcmb)**2 * exp(l * (l + 1) * sigma**2 / (8 * log(2)))
    if (type1 == 2 and type2 == 2) or (type1 == 3 and type2 == 3):
        noise = (Delta_P / Tcmb)**2 * exp(l * (l + 1) * sigma**2 / (8 * log(2)))

    return noise

cdef double cmbps_obs(double l, int type1, int type2, double sigma, double Delta_T, double Delta_P) noexcept nogil:
    # observed powerspectrum (with noise)
    return cmbps_l(l, type1, type2) + cmbps_noise(l, type1, type2, sigma, Delta_T, Delta_P)

cdef double f(double l1x, double l1y, double l2x, double l2y, int type1, int type2) noexcept nogil:
    # angle is defined as: angle of l1 vec - angle of l2 vec
    # L = l1 + l2

    cdef double l1 = sqrt(l1x**2 + l1y**2)
    cdef double l2 = sqrt(l2x**2 + l2y**2)
    cdef double Ll1 = l1**2 + l2x*l1x + l2y*l1y
    cdef double Ll2 = l1x*l2x + l1y*l2y + l2**2
    cdef double angle12 = atan2(l1y, l1x) - atan2(l2y, l2x)

    if type1 == 1 and type2 == 1:
        return cmbps_ul(l1, 1, 1) * Ll1 + cmbps_ul(l2, 1, 1) * Ll2
    if (type1 == 1 and type2 == 2) or (type1 == 2 and type2 == 1):
        # used to be cmbps_ul(l1, 1, 2) * Ll1 * cos(angle12) + cmbps_ul(l2, 1, 2) * Ll2
        return cmbps_ul(l1, 1, 2) * Ll1 * cos(2 * angle12) + cmbps_ul(l2, 1, 2) * Ll2
    if (type1 == 1 and type2 == 3) or (type1 == 3 and type2 == 1):
        return cmbps_ul(l1, 1, 2) * sin(2 * angle12) * Ll1
    if type1 == 2 and type2 == 2:
        return (cmbps_ul(l1, 2, 2) * Ll1 + cmbps_ul(l2, 2, 2) * Ll2 ) * cos(2 * angle12)
    if (type1 == 2 and type2 == 3) or (type1 == 3 and type2 == 2):
        return (cmbps_ul(l1, 2, 2) * Ll1 - cmbps_ul(l2, 3, 3) * Ll2 ) * sin(2 * angle12)
    if type1 == 3 and type2 == 3:
        return (cmbps_ul(l1, 3, 3) * Ll1 + cmbps_ul(l2, 3, 3) * Ll2 ) * cos(2 * angle12)
    else:
        pass

cdef double F(double l1x, double l1y, double l2x, double l2y, int type1, int type2, double sigma, double Delta_T, double Delta_P) noexcept nogil:
    cdef double l1 = sqrt(l1x**2 + l1y**2)
    cdef double l2 = sqrt(l2x**2 + l2y**2)
    cdef double result, numerator, denominator
    # angle is defined as: angle of l1 vec - angle of l2 vec
    if type1 == type2:
        numerator = f(l1x, l1y, l2x, l2y, type1, type2)
        denominator = 2 * cmbps_obs(l1, type1, type1, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type1, sigma, Delta_T, Delta_P)
    elif (type1 == 1 and type2 == 3) or (type1 == 2 and type2 == 3):
        numerator = f(l1x, l1y, l2x, l2y, type1, type2)
        denominator = cmbps_obs(l1, type1, type1, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type2, type2, sigma, Delta_T, Delta_P)
    else:
        numerator = cmbps_obs(l1, type2, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type1, sigma, Delta_T, Delta_P) * f(l1x, l1y, l2x, l2y, type1, type2) - cmbps_obs(l1, type1, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type2, sigma, Delta_T, Delta_P) * f(l2x, l2y, l1x, l1y, type1, type2)

        denominator = cmbps_obs(l1, type1, type1, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type2, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l1, type2, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type1, sigma, Delta_T, Delta_P) - (cmbps_obs(l1, type1, type2, sigma, Delta_T, Delta_P) * cmbps_obs(l2, type1, type2, sigma, Delta_T, Delta_P))**2

    result = numerator / denominator
            
    return result

cdef double simp(double[:] data, int n, double xrange) noexcept nogil:
    """
    Compute the integral of 'data' from 0..n-1 using Simpson's rule,
    assuming the total range length is 'xrange', so dx = xrange / (n - 1).
    """
    cdef double dx = xrange / (n - 1)
    cdef double total = data[0] + data[n-1]
    cdef int i

    for i in range(1, n - 1):
        if (i & 1) == 1:   # i odd
            total += 4.0 * data[i]
        else:              # i even
            total += 2.0 * data[i]
    
    return dx * total / 3.0

cdef double simp2d(double[:, :] data, int n, double xrange) noexcept nogil:
    # integral for square with length xrange

    cdef double dx = xrange / (n - 1)
    cdef double total = simp(data[0, :], n, xrange)  + simp(data[n-1, :], n, xrange) 
    cdef int i

    for i in range(1, n - 1):
        if (i & 1) == 1:   # i odd
            total += 4.0 * simp(data[i, :], n, xrange) 
        else:              # i even
            total += 2.0 * simp(data[i, :], n, xrange) 
    
    return dx * total / 3.0

cpdef double simp2dp(double[:, :] data, int n, double xrange) noexcept nogil:
    # integral for square with length xrange

    cdef double dx = xrange / (n - 1)
    cdef double total = simp(data[0, :], n, xrange)  + simp(data[n-1, :], n, xrange) 
    cdef int i

    for i in range(1, n - 1):
        if (i & 1) == 1:   # i odd
            total += 4.0 * simp(data[i, :], n, xrange) 
        else:              # i even
            total += 2.0 * simp(data[i, :], n, xrange) 
    
    return dx * total / 3.0

cdef double A_integrand(double l1x, double l1y, double L, int type1, int type2, double sigma, double Delta_T, double Delta_P) noexcept nogil:
    # # Avoid potential issues at the origin:
    # if np.sqrt(l1x**2 + l1y**2) == 0:
    #     return 0
    # l2 = L - l1, we take L to be in x direction wlog
    cdef double l2x = L - l1x 
    cdef double l2y = 0. - l1y
    cdef double l1 = sqrt(l1x**2 + l1y**2)
    cdef double l2 = sqrt(l2x**2 + l2y**2)
    if l1 < 2 or l2 < 2 or l1 > 10000 or l2 > 10000:
        return 0.

    else:
        return f(l1x, l1y, l2x, l2y, type1, type2) * F(l1x, l1y, l2x, l2y, type1, type2, sigma, Delta_T, Delta_P)

def A_integrand_p(l1x, l1y, L, type1, type2, sigma, Delta_T, Delta_P):
    return A_integrand(l1x, l1y, L, type1, type2, sigma, Delta_T, Delta_P)

cdef double A(double L, int type1, int type2, double sigma, double Delta_T, double Delta_P, double lmaxint, int num_samples):
    
    cdef double[:, :] int_data = np.zeros((num_samples, num_samples))
    cdef double[:] ls = np.linspace(-1*lmaxint, lmaxint, num_samples)
    cdef int i, j

    for i in range(num_samples):
        for j in range(num_samples):
            int_data[i, j] = A_integrand(ls[i], ls[j], L, type1, type2, sigma, Delta_T, Delta_P)

    return L**2 * (2 * np.pi)**2 / simp2d(int_data, num_samples, 2 * lmaxint)

cdef double N_integrand(double l1x, double l1y, double L, int type11, int type12, int type21, int type22, double sigma, double Delta_T, double Delta_P) noexcept nogil:
    cdef double l2x = L - l1x 
    cdef double l2y = 0. - l1y

    cdef double l1 = sqrt(l1x**2 + l1y**2)
    cdef double l2 = sqrt(l2x**2 + l2y**2)

    cdef double term1, term2

    if l1 < 2 or l2 < 2 or l1 > 10000 or l2 > 10000:
        return 0.
    else:
        term1 = (F(l1x, l1y, l2x, l2y, type11, type12, sigma, Delta_T, Delta_P) *
                F(l1x, l1y, l2x, l2y, type21, type22, sigma, Delta_T, Delta_P) *
                cmbps_obs(l1, type11, type12, sigma, Delta_T, Delta_P) *
                cmbps_obs(l2, type21, type22, sigma, Delta_T, Delta_P))
        term2 = (F(l1x, l1y, l2x, l2y, type11, type12, sigma, Delta_T, Delta_P) * 
                F(l2x, l2y, l1x, l1y, type21, type22, sigma, Delta_T, Delta_P) *
                cmbps_obs(l1, type11, type22, sigma, Delta_T, Delta_P) *
                cmbps_obs(l2, type21, type12, sigma, Delta_T, Delta_P))
        
        return term1 + term2

cpdef N(double L, int type11, int type12, int type21, int type22, double sigma, double Delta_T, double Delta_P, double lmaxint, int num_samples):
    cdef double[:, :] int_data = np.zeros((num_samples, num_samples))
    cdef double[:] ls = np.linspace(-1*lmaxint, lmaxint, num_samples)
    cdef int i, j

    # print(L, type11, type12, type21, type22)
    if type11 == type21 and type12 == type22:
        return A(L, type11, type12, sigma, Delta_T, Delta_P, lmaxint, num_samples)

    else:
        for i in range(num_samples):
            for j in range(num_samples):
                int_data[i, j] = N_integrand(ls[i], ls[j], L, type11, type12, type21, type22, sigma, Delta_T, Delta_P)
        
        return L**(-2) * A(L, type11, type12, sigma, Delta_T, Delta_P, lmaxint, num_samples) * A(L, type21, type22, sigma, Delta_T, Delta_P, lmaxint, num_samples) * (2 * np.pi)**(-2) * simp2d(int_data, num_samples, 2 * lmaxint)