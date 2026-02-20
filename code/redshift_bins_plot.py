import numpy as np
from cosm_setup import lensing_spectra

import matplotlib.pyplot as plt

if __name__ == '__main__':
    # Define chi values (comoving distance)
    chi = np.linspace(0, 8000, 1000)  # Adjust range as needed
    # Create instance of lensing_spectra class
    cosm = lensing_spectra()
    
    # Plot galaxy number density for each redshift bin
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for bin_number in range(1, cosm.number_bins + 1):
        density = [cosm.galaxy_density_chi_bin_fast(chi_val, bin_number) for chi_val in chi]
        ax.plot(chi, density, label=f'Bin {bin_number}')
    
    ax.set_xlabel('Comoving Distance χ (Mpc)')
    ax.set_ylabel('Galaxy Number Density')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('redshift_bins_fast.png')
    plt.show()