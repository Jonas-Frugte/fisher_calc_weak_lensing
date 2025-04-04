
# Fisher Calculations for Weak Lensing

This code calculates weak lensing powerspectra, bispectra, SNR's, and parameter constraints. It uses CAMB, for initial simulations and Fisher matrix formalism for SNR and parameter constraints. For further information on the formulas used please refer to "Parameter Forecasts From CMB Lensing and Galaxy Lensing Power- and Bispectra" by myself and Prof. Daan Meerburg (in preperation).

## Important files within code folder:
-	cosm_setup.py: Runs cosmological simulation and calculates power and bispectra for CMB and galaxy weak lensing, as well as other related functions. Written in pure Python and thus slow. Meant to be used for testing purposes and to check interpolation accuracy.
-	data_exporter.py: Exports folders with data used to recreate lensing power and bispectra through interpolation. Bispectra is 3 dimensional and thus not viable to interpolate directly, we instead interpolate all functions it is based on that are 2 or 1 dimensional. A new folder is created for each configuration of $\Lambda$CDM parameters. Multiple configurations are required to be able to differentiate the lensing spectra with respect to these parameters.
-	data_importer.pyx: Imports data earlier created and defines lensing spectra in same way as cosm_setup.py. This version is orders of magnitude faster and is used for the Fisher matrix calculations.
-	interpolation.pyx: interpolation methods written for 1D and 2D linearly and log-spaced grids. Significantly faster than using e.g. SciPy's interpolation methods.
-	Fisher_calc.pyx / Fisher_calc_python_imp.pyx / Fisher_matrix_powerspectra.py / Fisher_matrix.py / SNR.py: contain Fisher matrix calculations, to be streamlined later.
-	cmb_noise.py / cmb_noise_fast.pyx / cmb_noise_fast_runner.py: calculates CMB lensing reconstruction noise according to "Mass Reconstruction with CMB Polarization" by Hu & Okamoto, 2002.
-	setup.py: compiles .pyx files, can be run using `python setup.py build_ext --inplace`
