import cosm_setup as cs
import numpy as np
import matplotlib.pyplot as plt
spectra = cs.lensing_spectra()

# chi**(-6) * scale_factor(chi, scale_factor_data)**(-2) * ((chi - chi_source)  / (chi * chi_source))**2 * matter_power_spectrum(l / chi, z_at_chi(chi, z_at_chi_data), mps_data)

l = 100
chi_source = 12000
chis = np.linspace(0.0000001, 12000, 10000000)

def integrand(chi):
    return chi**(2) * spectra.scale_factor(chi)**(-2) * ((chi - chi_source)  / (chi * chi_source))**2 * spectra.mps(l / chi, spectra.z_at_chi(chi))

vals = [integrand(chi) for chi in chis]

plt.plot(chis, vals)
plt.savefig('test1.png')

plt.semilogx(chis, vals)
plt.savefig('test2.png')

plt.loglog(chis, vals)
plt.savefig('test3.png')
