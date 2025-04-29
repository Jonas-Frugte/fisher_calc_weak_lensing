print(f'Now running setup.py as: {__name__}.')

import camb
import numpy as np
import scipy

#======================================

matter_bispec_fit_params = (0.484, 3.740, -0.849, 0.392, 1.013, -0.575, 0.128, -0.722, -0.926)

class lensing_spectra:
  '''
  Calculates relevant features of a cosmology given input parameters and matter bispectrum fitting constants
  '''
  def __init__(self, H0=67.4, ombh2=0.0223, omch2=0.119, ns=0.965, As=2.13e-9, mnu=0.06, w0=-1, redshifts=np.linspace(1100, 0, 200), matter_bispec_fit_params = matter_bispec_fit_params):
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
    self.matter_bispec_fit_params = matter_bispec_fit_params # params of matter bispectrum fitting function as (a_1, ..., a_9)
    self.k_nl_max = 2000

    self.obtain_results()



    self.z_at_chi = self.results.redshift_at_comoving_radial_distance

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

    self.k_nl = self.find_k_nl()
    print(f"k_nl = {self.k_nl}.")


  def obtain_results(self):
    '''
    Sets cosmological results
    '''
    
    print('Obtaining cosmological results.')

    pars = camb.CAMBparams()
    lmax = 2000
    pars.set_cosmology(H0=self.H0, ombh2=self.ombh2, omch2=self.omch2, mnu=self.mnu, omk=0, tau=0.063, neutrino_hierarchy='normal')
    pars.InitPower.set_params(As=self.As, ns=self.ns, r=0)
    pars.InitPower.pivot_scalar = 0.05
    pars.set_for_lmax(lmax, lens_potential_accuracy=0)
    pars.set_matter_power(redshifts=self.redshifts, kmax=self.k_nl_max, nonlinear=True)
    pars.set_dark_energy(w=self.w0)
    self.pars = pars

    kmax_interp = 2000

    # power spectrum in units of Mpc^3, k in units of Mpc^-1, with input (z, k) on RHS 
    print('Creating linear mps')
    linear_mps_swapped = camb.get_matter_power_interpolator(
      self.pars,
      nonlinear=False, 
      hubble_units=False, 
      k_hunit=False,
      zmin = 0,
      zmax = 1100,
      kmax = kmax_interp,
      extrap_kmax=True).P
    
    print('Creating mps')
    mps_swapped = camb.get_matter_power_interpolator(
      self.pars,
      nonlinear=True, 
      hubble_units=False, 
      k_hunit=False,
      zmin = 0,
      zmax = 1100,
      kmax = kmax_interp,
      extrap_kmax=True).P
    
    self.mps = lambda k, z : mps_swapped(z, k)
    self.linear_mps = lambda k, z : linear_mps_swapped(z, k)

    print('Creating rest of results')
    self.results = camb.get_results(pars)

    print('Cosmological results created.')

    
    pass

  # def lps_cc_CAMB(self, lmax=10000):
  #   """
  #   Calculate the lensing potential power spectrum (C_phi_phi) using CAMB for the CMB (convergence) case.

  #   Parameters:
  #   -----------
  #   lmax : int
  #       Maximum multipole for the calculation.

  #   Returns:
  #   --------
  #   ell : ndarray
  #       Array of multipole moments.
  #   C_phi_phi : ndarray
  #       Array of lensing potential power spectrum values.
  #   """
  #   # Get the lensing potential power spectrum
  #   # if raw_cl = false then returns L^2(L+1)^2 / (2pi) times the result
  #   lens_potential_cls = self.results.get_lens_potential_cls(lmax=lmax, raw_cl=True)
  #   ell = np.arange(lmax + 1)
  #   C_phi_phi = lens_potential_cls[:, 0]

  #   return ell, C_phi_phi
  
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
    result = scipy.misc.derivative(
      lambda ln_k : np.log(self.linear_mps(np.exp(ln_k), z)),
      np.log(k),
      dx = 1e-4
    )
    # if result > 2:
    #   print('ns_eff > 2:', result, k, z)
    return result

  def sigma_8(self, z):
    '''
    Parameters:
    z : float

    Returns:
    sigma_8 : float
    '''
    # np.interp only accepts increase x values, hence why we need to invert the matrices, self.redshifts is in order of decreasing z values
    return np.interp(z, self.redshifts[::-1], self.results.get_sigma8()[::-1])
  
  def find_k_nl(self):
    error_tol = 0.1
    kmax = 1
    num_samples = int(kmax / 0.00001)  # stepsize

    ls = np.linspace(1e-2, kmax, num_samples)
    # toshiya uses factor of 4 pi here in paper but 1/(2 pi^2) agrees with his results
    vals = np.array([np.abs(1 - (2 * np.pi**2)**(-1) * l**3 * self.linear_mps(l, 0.)) for l in ls])
    
    vals_min_index = np.argmin(vals)
    error = vals[vals_min_index]
    result = ls[vals_min_index]

    if vals_min_index == len(ls) - 1:
        raise ValueError("Minimum found at maximum index; consider increasing kmax or reducing stepsize.")
    if error > error_tol:
        raise ValueError(f"Error tolerance not reached: error = {error} > {error_tol}")

    return result
      
  def q(self, k):
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
    return k / self.k_nl

  def Q_3(self, n):
    '''
    used in matter bispec fitting func

    Parameters:
    n : float

    Returns:
    float
    '''
    result = (4 - 2**n) / (1 + 2**(n + 1))
    # if result < 0:
    #   print(result, n)
    return max(0., result)
  
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
  
  def lps(self, l, types):
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
    result =  9 * (1/l)**4 * self.omega_m**2 * self.H0**4 * self.lightspeed_kms**(-4) * scipy.integrate.quad(integrand, 1e-5, self.chi_last_scatter, epsabs=1e-30, epsrel=1e-3, limit = 50)[0]

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

    numerator = 1 + self.sigma_8(z)**a_6 * np.sqrt(0.7 * self.Q_3(self.ns_eff(k, z))) * (self.q(k) * a_1)**(self.ns_eff(k, z) + a_2)
    denominator = 1 + (self.q(k) * a_1)**(self.ns_eff(k, z) + a_2)
  
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

    numerator = 1 + 0.2 * a_3 * (self.ns_eff(k, z) + 3) * (self.q(k) * a_7) ** (self.ns_eff(k, z) + 3 + a_8)
    denominator = 1 + (self.q(k) * a_7) ** (self.ns_eff(k, z) + 3.5 + a_8)

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
    
    numerator = 1 + (4.5 * a_4 / (1.5 + (self.ns_eff(k, z) + 3)**4)) * (self.q(k) * a_5)**(self.ns_eff(k, z) + 3 + a_9)
    denominator = 1 + (self.q(k) * a_5)**(self.ns_eff(k, z) + 3.5 + a_9)

    return numerator / denominator
  
  def law_cosines(self, x, y, z):
    # gives cosine of angle between vector x and y, where we know the magnitudes of x, y, z and that x + y + z = 0 vector
    return -1 * (x**2 + y**2 - z**2) / (2 * x * y)
  
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
    # set a = b = c = 1 and use linear matter power spectrum
    term_1 = (5. / 7.)
    term_2 = 0.5 * self.law_cosines(k1, k2, k3) * (k1 / k2 + k2 / k1)
    term_3 = (2. / 7.) * self.law_cosines(k1, k2, k3)**2
    return term_1 + term_2 + term_3
  
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
  
  def mbs(self, k_1, k_2, k_3, z, model = 'nonlinear'):
    '''
    Parameters:
    k_1, k_2, k_3 : floats
      wavenumbers
    z : float
      redshift

    Returns:
    float
    '''
    if model == 'nonlinear':
      term = lambda p_1, p_2, p_3 : 2 * self.F_2(p_1, p_2, p_3, z) * self.mps(p_1, z) * self.mps(p_2, z)

    if model == 'tree' or model == 'linear':
      term = lambda p_1, p_2, p_3 : 2 * self.F_2_tree(p_1, p_2, p_3, z) * self.linear_mps(p_1, z) * self.linear_mps(p_2, z)

    return term(k_1, k_2, k_3) + term(k_2, k_3, k_1) + term(k_3, k_1, k_2) # mfw the permutation is cyclic
  
  # approximate expression of wigner_3j symbol, same as used by toshiya (2016)
  def wigner_3j_approx_nocheck(self, l1, l2, l3):

    l1 = float(l1)
    l2 = float(l2)
    l3 = float(l3)

    L = ( l1 + l2 + l3 ) / 2

    # does not check if 2L is even or triangle inequalities, does so in lbs_f and lbs_der directly instead to save on computing bispec if result should be zero anyway

    factor = (-1)**L * np.sqrt(np.e**3 / (2 * np.pi)) * (L + 1)**(-0.25)
    term = lambda li : (L-li+1)**(-0.25) * ( (L-li+0.5) / (L-li+1) )**(L-li+0.25)
    return factor * term(l1) * term(l2) * term(l3)

  def lbs_flat(self, l1, l2, l3, types, model = 'nonlinear'):
    '''
    lensing bispectrum
    
    Parameters:
    k_1, k_2, k_3: floats
        wavenumbers
    types: iterable of length 3 containing 'shear' and 'convergence'
    
    Returns:
    float
    '''

    # convert to float to prevent integer overflow problems
    l1 = float(l1)
    l2 = float(l2)
    l3 = float(l3)

    def integrand(chi):
        return chi**2 * self.scale_factor(chi)**(-3) * self.window_func(chi, types[0]) * self.window_func(chi, types[1]) * self.window_func(chi, types[2]) \
        * self.mbs(l1/chi, l2/chi, l3/chi, self.results.redshift_at_comoving_radial_distance(chi), model = model)
    
    integration_result = scipy.integrate.quad(integrand, 0, self.chi_last_scatter, limit = 500, epsabs = 1e-40, epsrel = 1e-3)[0]

    fraction_factor = 8 / (l1**2 * l2**2 * l3**2)
    const_factor = self.rho_bar ** 3
    
    return fraction_factor * const_factor * integration_result

  def lbs(self, l1, l2, l3, types, model = 'nonlinear'):
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

      wigner_factor = np.abs(self.wigner_3j_approx_nocheck(l1, l2, l3)) # TODO: check if this is okay for nondiag components of Fisher matrices
      sqrt_factor = np.sqrt((2.0 * l1 + 1) * (2.0 * l2 + 1) * (2.0 * l3 + 1) / (4 * np.pi))
      
      return wigner_factor * sqrt_factor * self.lbs_flat(l1, l2, l3, types, model = model)
    else:
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