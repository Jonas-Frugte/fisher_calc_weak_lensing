import cosm_setup as cs
import matplotlib.pyplot as plt
import numpy as np

tau_vals = [0.063*0.85, 0.063 * 0.9, 0.063 * 0.95, 0.063, 0.063 * 1.05, 0.063 * 1.1, 0.063 * 1.15]
As_fid = 2.13e-9
As_vals = [As_fid*0.85, As_fid*0.9, As_fid*0.95, As_fid, As_fid * 1.05, As_fid * 1.1, As_fid * 1.15]
ls = np.arange(2000+1)

fig, axs = plt.subplots()

spectra_fid = cs.lensing_spectra()
cls_fid = spectra_fid.results.get_lensed_scalar_cls(lmax=2000, raw_cl=True)

def tau_plot():
    for tau in tau_vals:
        spectra = cs.lensing_spectra(tau = tau)
        cls = spectra.results.get_lensed_scalar_cls(lmax=2000, raw_cl=True)
        axs.semilogx(ls[2:], cls[2:,0] / cls_fid[2:,0], label = f'TT, tau = {tau}')

    fig.legend()
    fig.savefig('taus_perturbed.png')
    pass

def As_plot():
    for As in As_vals:
        spectra = cs.lensing_spectra(As = As)
        cls = spectra.results.get_lensed_scalar_cls(lmax=2000, raw_cl=True)
        axs.semilogx(ls[2:], cls[2:,0] / cls_fid[2:,0], label = f'TT, As = {As}')

    fig.legend()
    fig.savefig('As_perturbed.png')
    pass
