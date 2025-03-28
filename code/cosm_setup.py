print(f'now running setup.py as {__name__}')

import camb
from camb import model
from scipy import integrate
import math
import numpy as np
import scipy
import matplotlib.pyplot as plt
import scipy
import scipy.integrate

# np.seterr(all='warn')

#======================================

matter_bispec_fit_params = (0.484, 3.740, -0.849, 0.392, 1.013, -0.575, 0.128, -0.722, -0.926)

class lensing_spectra:
  '''
  Calculates relevant features of a cosmology given input parameters and matter bispectrum fitting constants
  '''
  def __init__(self, H0=67.4, ombh2=0.0224, omch2=0.120, ns=0.965, As=2.1e-9, mnu=0.06, w0=-1, wa=0, redshifts=np.linspace(1100, 0, 200), matter_bispec_fit_params = matter_bispec_fit_params, fiducial_k_nls = False):
    '''
    Initialize instance of class, set results and some methods/attributes, default inputs are for fiducial cosmology

    Parameters:
    H0 : float
      Hubble constant
    ombh2 : float
      Baryon density
    omch2 : float
      Cold dark matter density
    ns : float
     spectral index
    redshift : array
      array of redshift values used in obtaining results
    As : float
      scalar amplitude
    mnu : float
      sum of neutrino masses
    matter_bispec_fit_params : iterable of length 9
      defines (a_1, ..., a_9) as used in matter bispec fitting function paper
    '''
    self.H0 = H0          # Hubble constant
    self.h = H0 * 0.01    # Unitless Hubble constant
    self.ombh2 = ombh2    # Baryon density
    self.omch2 = omch2    # Cold dark matter density
    self.ns = ns          # Spectral index
    self.redshifts = redshifts  # Array of redshift values
    self.As = As          # Scalar amplitude
    self.mnu = mnu        # Sum of neutrino masses
    self.w0 = w0
    self.wa = wa
    self.matter_bispec_fit_params = matter_bispec_fit_params # params of matter bispectrum fitting function as (a_1, ..., a_9)

    self.lps_count = 0

    self.window_func_normalization = 200. # see window functions


    self.obtain_results()

    self.z_at_chi = self.results.redshift_at_comoving_radial_distance

    # power spectrum in units of Mpc^3, k in units of Mpc^-1, with input (z, k) on RHS 
    mps_swapped = self.results.get_matter_power_interpolator(nonlinear=True, hubble_units=False, k_hunit=False).P
    self.mps = lambda k, z : mps_swapped(z, k)
    linear_mps_swapped = self.results.get_matter_power_interpolator(nonlinear=False, hubble_units=False, k_hunit=False).P
    self.linear_mps = lambda k, z : linear_mps_swapped(z, k)

    self.lightspeed_kms = 299792 # speed of light in km / s
    self.omega_m = (self.ombh2 + self.omch2) / self.h**2 # matter density parameter
    self.rho_bar =  3 * self.omega_m * H0**2 / (2 * self.lightspeed_kms**2) # this is in natural units, used to be called self.C

    self.z_last_scatter = self.results.get_derived_params()['zstar']
    self.chi_last_scatter = self.results.comoving_radial_distance(self.z_last_scatter)

    # galaxy density dist as a function of redshift
    f = lambda z : float(z)**2
    z0 = 0.64
    beta = 1.5
    normalization = scipy.integrate.quad(lambda z : f(z) * np.exp(float(-(z / z0)**beta)), 1e-8, self.results.get_derived_params()['zstar'])[0]

    self.galaxy_density = lambda z : f(z) * np.exp(-(max(1e-12, z / z0))**beta) / normalization # normalizes integral of func to 1

    # source redshift and comoving dist of galaxies using Luna's convention
    self.z_galaxy_source = 5
    self.chi_galaxy_source = self.results.comoving_radial_distance(self.z_galaxy_source)

    k_nls_fiducial = [0.16800002240000297, 0.2560000341333379, 0.40400005386667387, 0.6560000874666784, 1.108000147733353, 1.9480002597333679, 3.57200047626673, 6.792000905600121, 13.384001784533572, 27.320003642667153, 58.656007820801044, 4235.396564719542, 6730.844897445986, 7977.9530637270755, 9005.469200729227, 10325.777376770317, 11815.08157534421, 13027.237736965031, 13962.933861724516, 14637.293951639194, 15099.270013236002, 15421.03005613734, 15688.574091809878, 15973.122129749618, 16283.550171140023, 16617.07821561043, 16970.694262759236, 17341.130312150708, 17724.786363304847, 18117.738415698455, 18515.666468755528, 18913.874521849935, 19307.262574301676, 19690.338625378485, 20059.042674539025, 20413.55872180783, 20754.882767317704, 21084.162811221708, 21402.68685369158, 21711.894894919318, 22013.346935112924, 22308.722974496395, 22599.831013310803, 22888.587051811606, 23177.011090268144, 23466.487128864952, 23756.875167583356, 24047.831206377494, 24338.995245199367, 24629.99128399884, 24920.42332272311, 25209.887361318313, 25497.963399728455, 25784.211437894857, 26068.183475757796, 26349.4155132554, 26627.63955035194, 26902.931587057545, 27175.42762339035, 27445.275659370087, 27712.635695018093, 27977.675730356765, 28240.575765410103, 28501.523800203173, 28760.731834764243, 29018.403869120517, 29274.76790330239, 29529.959937327993, 29783.963971195197]
    zs_used_fiducial = [0.0, 0.5050505050505051, 1.0101010101010102, 1.5151515151515151, 2.0202020202020203, 2.5252525252525255, 3.0303030303030303, 3.5353535353535355, 4.040404040404041, 4.545454545454546, 5.050505050505051, 5.555555555555556, 6.0606060606060606, 6.565656565656566, 7.070707070707071, 7.575757575757576, 8.080808080808081, 8.585858585858587, 9.090909090909092, 9.595959595959597, 10.101010101010102, 10.606060606060607, 11.111111111111112, 11.616161616161618, 12.121212121212121, 12.626262626262626, 13.131313131313131, 13.636363636363637, 14.141414141414142, 14.646464646464647, 15.151515151515152, 15.656565656565657, 16.161616161616163, 16.666666666666668, 17.171717171717173, 17.67676767676768, 18.181818181818183, 18.68686868686869, 19.191919191919194, 19.6969696969697, 20.202020202020204, 20.70707070707071, 21.212121212121215, 21.71717171717172, 22.222222222222225, 22.72727272727273, 23.232323232323235, 23.73737373737374, 24.242424242424242, 24.747474747474747, 25.252525252525253, 25.757575757575758, 26.262626262626263, 26.767676767676768, 27.272727272727273, 27.77777777777778, 28.282828282828284, 28.78787878787879, 29.292929292929294, 29.7979797979798, 30.303030303030305, 30.80808080808081, 31.313131313131315, 31.81818181818182, 32.323232323232325, 32.82828282828283, 33.333333333333336, 33.83838383838384, 34.343434343434346]

    if fiducial_k_nls:
      # saves about 20 minutes or recalculating k_nls every time for fiducial case
      print("Took fiducial values for k_nls")
      self.zs_for_k_nls_used_max = zs_used_fiducial[-1]
      self.k_nls = k_nls_fiducial
      self.zs_for_k_nls_used = zs_used_fiducial

    else:
      print('Generating k_nl vals')
      self.zs_for_k_nls = np.linspace(0, 300, 1000)
      self.k_nls, self.zs_for_k_nls_used = self._k_nl_exact_nonfid(self.zs_for_k_nls, 30000, k_nls_fiducial, zs_used_fiducial)
      self.zs_for_k_nls_used_max = self.zs_for_k_nls_used[-1]
      print('Done generating k_nl vals')
      # print('k_nls : ', [float(k) for k in self.k_nls])
      # print('zs_used : ', self.zs_for_k_nls_used)

  def obtain_results(self):
    '''
    Sets cosmological results
    '''

    print('Obtaining cosmological results.')

    pars = camb.CAMBparams()
    lmax = 10000
    pars.set_cosmology(H0=self.H0, ombh2=self.ombh2, omch2=self.omch2, mnu=self.mnu, omk=0, tau=0.05)
    pars.InitPower.set_params(As=self.As, ns=self.ns, r=0)
    pars.set_for_lmax(lmax, lens_potential_accuracy=1)
    pars.set_matter_power(redshifts=self.redshifts, kmax=10000, nonlinear=True)
    pars.set_dark_energy(w=self.w0, wa=self.wa)

    self.results = camb.get_results(pars)

    print('Cosmological results created.')
    
    pass

  def lps_cc_CAMB(self, lmax=10000):
    """
    Calculate the lensing potential power spectrum (C_phi_phi) using CAMB for the CMB (convergence) case.

    Parameters:
    -----------
    lmax : int
        Maximum multipole for the calculation.

    Returns:
    --------
    ell : ndarray
        Array of multipole moments.
    C_phi_phi : ndarray
        Array of lensing potential power spectrum values.
    """
    # Get the lensing potential power spectrum
    # if raw_cl = false then returns L^2(L+1)^2 / (2pi) times the result
    lens_potential_cls = self.results.get_lens_potential_cls(lmax=lmax, raw_cl=True)
    ell = np.arange(lmax + 1)
    C_phi_phi = lens_potential_cls[:, 0]

    return ell, C_phi_phi
  
  def _k_nl(self, z):
    # used quadratic interpolation for finer k_nl values but this sometimes gives negative k_nl values resulting in errors
    return scipy.interpolate.interp1d(self.zs_for_k_nls_used, self.k_nls, kind = 'linear', bounds_error=False, fill_value = 'extrapolate')(z)
  
  def ns_eff(self, k, z):
    '''
    Parameters:
    k : float
      wavenumber
    z : float
      redshift

    Returns:
    ns_eff : float
      effective spectral index defined as d ln P^lin_m(k) / d ln k
    '''
    return scipy.misc.derivative(
      lambda ln_k : np.log(self.linear_mps(np.exp(ln_k), z)),
      np.log(k),
      dx = 1e-4
    )

  def _sigma_8(self, z):
    '''
    Parameters:
    z : float

    Returns:
    sigma_8 : float
    '''
    # np.interp only accepts increase x values, hence why we need to invert the matrices, self.redshifts is in order of decreasing z values
    return np.interp(z, self.redshifts[::-1], self.results.get_sigma8()[::-1])
  
  def _k_nl_exact(self, zs):
    errors = []
    results = []
    zs_used = []
    for z in zs:
      # if z > 15:
      #   raise Exception('does not work for z > 15')
      error_tol = 0.8
      kmax = 30000
      num_samples = int(kmax / 0.004) # denominator is stepsize in ls, used to be 0.001 which didn't give any issues

      ls = np.linspace(0, kmax, num_samples)
      vals = np.array([np.abs(1 - l**3 * self.linear_mps(l, z) / (2 * np.pi**2)) for l in ls])
      vals_min_index = np.argmin(vals)
      error = vals[vals_min_index]
      result = ls[vals_min_index]

      if vals_min_index == num_samples - 1:
        z_k_nl_max = zs_used[-1]
        print(f"Max z reached. \n k_nls = {results}, zs = {zs_used}")
        return results, errors, zs_used, z_k_nl_max
      
      zs_used.append(z)

      if error > error_tol:
        raise Exception(f'did not reach error tolerance, error is {error}, z = {z}, best val is {result}')
      
      errors.append(error)
      results.append(result)
      # print(f'k created: {result}')
    
    return results, errors, zs_used, z_k_nl_max
      
    

  def _k_nl_exact_nonfid(self, zs, kmax, k_nls_fiducial, zs_used_fiducial):
    num_samples = 10 ** 4
    error_tol = 0.8
    results = []
    zs_used = []

    for z in zs:

      k_start = scipy.interpolate.interp1d(k_nls_fiducial, zs_used_fiducial, kind = 'linear', bounds_error=False, fill_value = 'extrapolate')(z)
      
      l_upper = k_start + 25
      l_lower = max(k_start - 25, 1e-8)
      ls = np.linspace(l_lower, l_upper, num_samples)
      vals = np.array([np.abs(1 - (4 * np.pi * l**3 * self.linear_mps(l, z))) for l in ls])
      vals_min_index = np.argmin(vals)
      error = vals[vals_min_index]
      result = ls[vals_min_index]
      
      #vals_min_index = 10**4 - 1
      vals_lower_min_index = 0
      vals_upper_min_index = 0
      if error > error_tol:
        iterations = 0
        # while error > error_tol or vals_lower_min_index == 0 or vals_lower_min_index == num_samples - 1 or vals_upper_min_index == 0 or vals_upper_min_index == num_samples - 1:
        while error > error_tol:
          if iterations > 100:
            raise Exception(f'more than 50 iterations required for z = {z}, l_lower = {l_lower}, l_upper = {l_upper}, error = {error}, result = {result}')
          iterations += 1
          ls_lower = np.linspace(max(l_lower - 50, 1e-8), l_lower, num_samples)
          ls_upper = np.linspace(max(l_upper, 1e-8), l_upper + 50, num_samples)

          vals_lower = np.array([np.abs(1 - (4 * np.pi * l**3 * self.linear_mps(l, z))) for l in ls_lower])
          vals_upper = np.array([np.abs(1 - (4 * np.pi * l**3 * self.linear_mps(l, z))) for l in ls_upper])

          vals_lower_min_index = np.argmin(vals_lower)
          vals_upper_min_index = np.argmin(vals_upper)
          
          error_lower = vals_lower[vals_lower_min_index]
          error_upper = vals_upper[vals_upper_min_index]

          error = min(error_lower, error_upper)
          if error_lower > error_upper:
            result = ls_upper[vals_upper_min_index]
          else:
            result = ls_lower[vals_lower_min_index]

          if l_lower > 1e-8 + 50:
            l_lower -= 50
          l_upper += 50

    # print('error achieved:',  error)
    # print('k_nl:', result, '\n')
      if result > kmax:
        break
      else:
        zs_used.append(z)
        results.append(result)
    
    return results, zs_used

      
  
  
  #   #return np.interp(z, np.linspace(1e-8, self.z_last_scatter, 50), self.k_nls, left = self.k_nls[0], right = self.k_nls[-1])

  def _q(self, k, z):
    '''
    used in matter bispec fitting func

    Parameters:
    k : float
      wavenumber
    z : float
      redshift
    
    Returns:
    float
    '''
    return k / self._k_nl(z)

  def _Q_3(self, n):
    '''
    used in matter bispec fitting func

    Parameters:
    n : float

    Returns:
    float
    '''
    return (4 - 2**n) / (1 + 2**(n + 1))
  
  def scale_factor(self, chi):
     '''
     Parameters
     chi : float
      comoving radial distance

    Returns:
    float
      scale factor
     '''
     return 1 / (1 + self.results.redshift_at_comoving_radial_distance(chi))
  
  def lps(self, l, types): # used to be "lensing_power_spectrum"
    '''
    Parameters:
    k : float
      wavenumber
    types : iterable of 'shear' or 'convergence' length 2
      denotes which lensing power spectrum to calculate
    
    Returns:
      float
    '''
    integrand = lambda chi : chi**2 * self.scale_factor(chi)**(-2) * self.window_func(chi, types[0]) * self.window_func(chi, types[1]) * self.mps(l/chi, self.results.redshift_at_comoving_radial_distance(chi))
    result =  9 * (1/l)**4 * self.omega_m**2 * self.H0**4 * self.lightspeed_kms**(-4) * scipy.integrate.quad(integrand, 1e-5, self.chi_last_scatter, epsabs=1e-3, epsrel=1e-3, limit = 50)[0]
    
    # print(self.lps_count)

    self.lps_count += 1

    return result

  def a(self, k, z):
    '''
    used in matter bisepc fitting func

    Parameters:
    k : float
      wavenumber
    z : float
      redshift

    Returns:
    float
    '''
    a_1 = self.matter_bispec_fit_params[0]
    a_2 = self.matter_bispec_fit_params[1]
    a_6 = self.matter_bispec_fit_params[5]
    if self._Q_3(self.ns_eff(k, z)) < 0:
      print(self._Q_3(self.ns_eff(k, z)), k, z)
    numerator = 1 + self._sigma_8(z)**a_6 * np.sqrt(0.7 * self._Q_3(self.ns_eff(k, z))) * (self._q(k, z) * a_1)**(self.ns_eff(k, z) + a_2)
    # except Exception as e:
    #   print(e)
    denominator = 1 + (self._q(k, z) * a_1)**(self.ns_eff(k, z) + a_2)
    
    return numerator / denominator
  
  def b(self, k, z):
    '''
    used in matter bisepc fitting func

    Parameters:
    k : float
      wavenumber
    z : float
      redshift

    Returns:
    float
    '''
    a_3 = self.matter_bispec_fit_params[2]
    a_7 = self.matter_bispec_fit_params[6]
    a_8 = self.matter_bispec_fit_params[7]

    numerator = 1 + 0.2 * a_3 * (self.ns_eff(k, z) + 3) * (self._q(k, z) * a_7) ** (self.ns_eff(k, z) + 3 + a_8)
    denominator = 1 + (self._q(k, z) * a_7) ** (self.ns_eff(k, z) + 3.5 + a_8)

    return numerator / denominator
  
  def c(self, k, z):
    '''
    used in matter bisepc fitting func

    Parameters:
    k : float
      wavenumber
    z : float
      redshift

    Returns:
    float
    '''
    a_4 = self.matter_bispec_fit_params[3]
    a_5 = self.matter_bispec_fit_params[4]
    a_9 = self.matter_bispec_fit_params[8]
    
    numerator = 1 + (4.5 * a_4 / (1.5 + (self.ns_eff(k, z) + 3)**4)) * (self._q(k, z) * a_5)**(self.ns_eff(k, z) + 3 + a_9)
    denominator = 1 + (self._q(k, z) * a_5)**(self.ns_eff(k, z) + 3.5 + a_9)

    return numerator / denominator
  
  def law_cosines(self, x, y, z):
    # gives cosine of angle between vector x and y, where we know the magnitudes of x, y, z and that x + y + z = 0 vector
    return -1 * (x**2 + y**2 - z**2) / (2 * x * y)
  
  def F_2(self, k1, k2, k3, z):
    '''
    effective kernel for matter bispec fitting func

    Parameters:
    k_1, k_2 : floats
      wavenumbers
    z : float
      redshift

    Returns:
    float
    '''
    term_1 = (5. / 7.) * self.a(k1, z) * self.a(k2, z)
    term_2 = 0.5 * self.law_cosines(k1, k2, k3) * (k1 / k2 + k2 / k1) * self.b(k1, z) * self.b(k2, z)
    term_3 = (2. / 7.) * self.law_cosines(k1, k2, k3)**2 * self.c(k1, z) * self.c(k2, z)
    return term_1 + term_2 + term_3
  
  def F_2_tree(self, k1, k2, k3, z):
    '''
    tree level kernel

    Parameters:
    k_1, k_2 : floats
      wavenumbers
    z : float
      redshift

    Returns:
    float
    '''
    # set a = b = c = 1
    term_1 = (5. / 7.)
    term_2 = 0.5 * self.law_cosines(k1, k2, k3) * (k1 / k2 + k2 / k1)
    term_3 = (2. / 7.) * self.law_cosines(k1, k2, k3)**2
    return term_1 + term_2 + term_3
  
  def mbs(self, k_1, k_2, k_3, z):
    '''
    Parameters:
    k_1, k_2, k_3 : floats
      wavenumbers
    z : float
      redshift

    Returns:
    float
    '''
    term = lambda p_1, p_2, p_3 : 2 * self.F_2(p_1, p_2, p_3, z) * self.mps(p_1, z) * self.mps(p_2, z)
    return term(k_1, k_2, k_3) + term(k_2, k_3, k_1) + term(k_3, k_1, k_2) # mfw the permutation is cyclic

  def mbs_tree(self, k_1, k_2, k_3, z):
    '''
    Parameters:
    k_1, k_2, k_3 : floats
      wavenumbers
    z : float
      redshift

    Returns:
    float
    '''
    term = lambda p_1, p_2, p_3 : 2 * self.F_2_tree(p_1, p_2, p_3, z) * self.mps(p_1, z) * self.mps(p_2, z)
    return term(k_1, k_2, k_3) + term(k_2, k_3, k_1) + term(k_3, k_1, k_2) # mfw the permutation is cyclic
  
  # approximate expression of wigner_3j symbol based on appendix in paper "Cosmological parameters from lensing power spectrum and bispectrum tomography"
  def wigner_3j_approx_nocheck(self, l1, l2, l3):

    l1 = float(l1)
    l2 = float(l2)
    l3 = float(l3)

    L = l1 + l2 + l3

    # does not check if L is even or triangle inequalities, does so in lbs_f and lbs_der directly instead to save on computing bispec if result should be zero anyway

    L_half = L / 2.0
    
    # Common factors
    factor = (-1)**L_half * (2 * 3.141592653589793)**(-0.5) * np.exp(3.0 / 2) * (L + 2)**(-0.25)
    
    # Power term for each fraction
    term1 = (L_half - l1 + 0.5) / (L_half - l1 + 1)
    term2 = (L_half - l2 + 0.5) / (L_half - l2 + 1)
    term3 = (L_half - l3 + 0.5) / (L_half - l3 + 1)

    # Raising to required powers
    term1_pow = term1 ** (L_half - l1 + 1.0 / 4.0)
    term2_pow = term2 ** (L_half - l2 + 1.0 / 4.0)
    term3_pow = term3 ** (L_half - l3 + 1.0 / 4.0)

    # The denominator terms
    denominator = ((L_half - l1 + 1)**0.25 * (L_half - l2 + 1)**0.25 * (L_half - l3 + 1)**0.25)
    
    # The final result
    return factor * term1_pow * term2_pow * term3_pow / denominator
  
  def lbs(self, l1, l2, l3, types):
    '''
    lensing bispectrum
    
    Parameters:
    k_1, k_2, k_3: floats
        wavenumbers
    types: iterable of length 3 containing 'shear' and 'convergence'
    
    Returns:
    float
    '''
    if (l1 + l2 + l3) % 2 == 0 and l3 <= l1 + l2 and l1 - l2 <= l3 and l2 - l1 <= l3:

      # convert to float to prevent integer overflow problems
      l1 = float(l1)
      l2 = float(l2)
      l3 = float(l3)

      def integrand(chi):
          return chi**2 * self.scale_factor(chi)**(-3) * self.window_func(chi, types[0]) * self.window_func(chi, types[1]) * self.window_func(chi, types[2]) \
          * self.mbs(l1/chi, l2/chi, l3/chi, self.results.redshift_at_comoving_radial_distance(chi))
      
      integration_result = scipy.integrate.quad(integrand, 0, self.chi_last_scatter, limit = 500, epsabs = 1e-02, epsrel = 1e-02)[0]

      wigner_factor = np.abs(self.wigner_3j_approx_nocheck(l1, l2, l3)) # TODO: check if this is okay
      sqrt_factor = np.sqrt((2.0 * l1 + 1) * (2.0 * l2 + 1) * (2.0 * l3 + 1) / (4 * np.pi))
      fraction_factor = 1 / (l1**2 * l2**2 * l3**2)
      const_factor = self.rho_bar ** 3 * 8

      # print(f"exact: Tuple (l1, l2, l3): ({l1}, {l2}, {l3})")
      # print(f"exact: Integration Result: {integration_result}")
      # print(f"exact: Wigner Factor: {wigner_factor}")
      # print(f"exact: Sqrt Factor: {sqrt_factor}")
      # print(f"exact: Fraction Factor: {fraction_factor}")
      # print(f"exact: Const Factor: {const_factor}")
      
      return wigner_factor * sqrt_factor * fraction_factor * const_factor * integration_result
    else:
      #raise ValueError(f"l's did not satisfy triangle inequalities: {(l1, l2, l3)}")
      return 0.
    
  def lbs_tree(self, l1, l2, l3, types):
    '''
    lensing bispectrum
    
    Parameters:
    k_1, k_2, k_3: floats
        wavenumbers
    types: iterable of length 3 containing 'shear' and 'convergence'
    
    Returns:
    float
    '''
    if (l1 + l2 + l3) % 2 == 0 and l3 <= l1 + l2 and l1 - l2 <= l3 and l2 - l1 <= l3:

      # convert to float to prevent integer overflow problems
      l1 = float(l1)
      l2 = float(l2)
      l3 = float(l3)

      def integrand(chi):
          return chi**2 * self.scale_factor(chi)**(-3) * self.window_func(chi, types[0]) * self.window_func(chi, types[1]) * self.window_func(chi, types[2]) \
          * self.mbs_tree(l1/chi, l2/chi, l3/chi, self.results.redshift_at_comoving_radial_distance(chi))
      
      integration_result = scipy.integrate.quad(integrand, 0, self.chi_last_scatter, limit = 500, epsabs = 1e-02, epsrel = 1e-02)[0]

      wigner_factor = np.abs(self.wigner_3j_approx_nocheck(l1, l2, l3)) # TODO: check if this is okay
      sqrt_factor = np.sqrt((2.0 * l1 + 1) * (2.0 * l2 + 1) * (2.0 * l3 + 1) / (4 * np.pi))
      fraction_factor = 1 / (l1**2 * l2**2 * l3**2)
      const_factor = self.rho_bar ** 3 * 8

      # print(f"exact: Tuple (l1, l2, l3): ({l1}, {l2}, {l3})")
      # print(f"exact: Integration Result: {integration_result}")
      # print(f"exact: Wigner Factor: {wigner_factor}")
      # print(f"exact: Sqrt Factor: {sqrt_factor}")
      # print(f"exact: Fraction Factor: {fraction_factor}")
      # print(f"exact: Const Factor: {const_factor}")
      
      return wigner_factor * sqrt_factor * fraction_factor * const_factor * integration_result
    else:
      #raise ValueError(f"l's did not satisfy triangle inequalities: {(l1, l2, l3)}")
      return 0.
    
  
  
  def window_func(self, chi, type):
    '''
    window functions following general def W(chi) := int_{chi}^{chi_{max}} d chi' p(chi') (chi' - chi)(chi'chi)

    Parameters:
    chi : float
      comoving radial distance
    type : string starting with 's' for 'shear' or 'c' for 'convergence'
      which window function to use

    Returns:
    float    
    '''
    if type[0] == 'c':
      return (self.chi_last_scatter - chi)/(self.chi_last_scatter * chi)
    elif type[0] == 's':
      return scipy.integrate.quad(
          lambda chi_prime : self.galaxy_density(self.results.redshift_at_comoving_radial_distance(chi_prime)) * (chi_prime - chi) / (chi_prime * chi),
          chi, self.chi_galaxy_source, epsabs = 1e-3, epsrel = 1e-3)[0]
    else:
      raise Exception('Unknown window function type, string needs to start with s or c')

if __name__ == '__main__':
  # for testing purposes
  results = lensing_spectra()
  print(results.window_func(200, 's'))