import numpy as np
import fisherstuff

visbpp_f = fisherstuff.visb_f # + fisherstuff.visp_f

s8_ders = np.array([2.76945054e-03, -6.50570238e+00,  4.30234556e+00,  3.06272435e-01, 1.93167857e+08, -2.15681667e-01, -2.01175433e-01])

s8_const = np.sum(np.outer(s8_ders, s8_ders) * visbpp_f)
print(s8_const)

print((visbpp_f[1, 1] + 2 * visbpp_f[1, 2] + visbpp_f[2, 2])**(-1/2))