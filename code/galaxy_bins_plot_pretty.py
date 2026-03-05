import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import scipy.optimize

# -----------------------
# galaxy distribution
# -----------------------

z0 = 0.64
beta = 1.5

def n_of_z(z):
    return z**2 * np.exp(-(z/z0)**beta)

# normalize
norm = scipy.integrate.quad(n_of_z, 0, 5)[0]
def pz(z):
    return n_of_z(z)/norm


# -----------------------
# cumulative distribution
# -----------------------

def cdf(z):
    return scipy.integrate.quad(pz, 0, z)[0]


# -----------------------
# bin edges (quantiles)
# -----------------------

N_bins = 4
bin_edges = [0]

for i in range(1, N_bins):
    q = i/N_bins
    edge = scipy.optimize.brentq(lambda z: cdf(z)-q, 0, 4)
    bin_edges.append(edge)

bin_edges.append(3.0)
bin_edges = np.array(bin_edges)


# -----------------------
# photo-z error
# -----------------------

def sigma_z(z):
    return 0.05*(1+z)


def gaussian(z, zp):
    sig = sigma_z(zp)
    return np.exp(-0.5*((z-zp)/sig)**2)/(sig*np.sqrt(2*np.pi))


# -----------------------
# smeared bin distributions
# -----------------------

z_grid = np.linspace(0,3,400)

bin_curves = []

for i in range(N_bins):

    def integrand(zp, z):
        if bin_edges[i] <= zp <= bin_edges[i+1]:
            return pz(zp)*gaussian(z, zp)
        return 0

    curve = []

    for z in z_grid:
        val = scipy.integrate.quad(integrand,
                                   bin_edges[i],
                                   bin_edges[i+1],
                                   args=(z,))[0]
        curve.append(val)

    bin_curves.append(np.array(curve))


# -----------------------
# plotting
# -----------------------

plt.figure(figsize=(7,4))

# total distribution
pz_vals = pz(z_grid)

plt.fill_between(z_grid, pz_vals, color="lightgray")
plt.plot(z_grid, pz_vals, color="black", lw=2)

# dashed bin distributions
for curve in bin_curves:
    plt.plot(z_grid, curve, "--", color="black")

# vertical bin edges (only up to curve)
for edge in bin_edges[1:-1]:
    height = pz(edge)
    plt.plot([edge, edge], [0, height], color="black")

plt.text(0.4, 0.3, str(1), ha="center", va="center", fontsize=14)
plt.text(0.75, 0.5, str(2), ha="center", va="center", fontsize=14)
plt.text(1.05, 0.45, str(3), ha="center", va="center", fontsize=14)
plt.text(1.6, 0.18, str(4), ha="center", va="center", fontsize=14)

    

plt.xlabel("z")
plt.ylabel(r"$p_{\rm gal}(z)$")
plt.xlim(0,3)
plt.ylim(0, 1.1*pz(0.64))

plt.tight_layout()
plt.savefig("/home3/p319950/ResearchProject/fisher_calc_weak_lensing/paper/figures/galaxy_bins.pdf")