"""
functions for the cold matter halo density profile
"""

#packages needed
import numpy as np
import scipy 
from scipy import optimize

#own packages
from cosmology.basic_cosmology import *
from cosmology.overdensities import *
from cosmology.variance import *

#find z_formation given by Mead 2020 eq. 21
def func_z_formation(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma, f = 0.01):
    """
    k_sigma is in units of h/Mpc, PS_sigma in (Mpc/h)^3 and M in solar_mass/h
    returns the formation redshift of a halo which is defined as
    in HMCode2020: https://arxiv.org/abs/2009.01858 in eq. 21
    """
    z = cosmo_dic['z']
    Omega_m_0 = cosmo_dic['Omega_m_0']
    Omega_w_0 = cosmo_dic['Omega_w_0']
    D_norm =  func_D_z_norm(z, Omega_m_0, Omega_w_0)
    def func_find_root(x, Mass):
        nu = func_nu(f*Mass, k_sigma, PS_sigma, Omega_0_sigma)
        return func_D_z_norm(x, Omega_m_0, Omega_w_0) - D_norm * nu
    if isinstance(M, (int, float)) == True:
        #test if we find a root, if not by definition formation redshift is set to given z.
        if func_find_root(z, M)*func_find_root(100., M) > 0.:
            return z
        else:
            z_f = optimize.brentq(func_find_root, z, 100., args = (M))
            return z_f
    else:
        #test also if we can find a root
        return np.array([optimize.brentq(func_find_root, z, 100., args=(m)) if func_find_root(z, m)*func_find_root(100., m) < 0 else z for m in M])



def func_conc_param(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma, axion_dic=None):
    """
    k_sigma is in units of h/Mpc, PS_sigma in (Mpc/h)^3 and M in solar_mass/h
    NOTE: Omega_0 must match with chosen PS_sigma
    returns the concentration parameter as defined in
    https://arxiv.org/abs/2009.01858 in eq. 20
    """
    B = 4.#5.196
    if 'gamma_1' in cosmo_dic and 'gamma_2' in cosmo_dic:
        # Implementing 2111.01199 Eq. (33) for mixed dark matter
        #Only apply correction for cold DM halos > M_cut
        if axion_dic == 'ignore':
            #print('Assuming axions affect arbitrarily low-mass cold halos')
            correction_array = np.ones_like(M)
        else:
            #cold_halo_cut = axion_dic['M_cut'] * cosmo_dic['omega_d_0'] / cosmo_dic['omega_ax_0']
            correction_array = (M > axion_dic['M_cut'])

        gamma_1 = cosmo_dic['gamma_1']
        gamma_2 = cosmo_dic['gamma_2']
        f_ax = cosmo_dic['omega_ax_0'] / (cosmo_dic['omega_ax_0']+cosmo_dic['omega_d_0'])
        M0 = 1.6e10 * (cosmo_dic['m_ax']/1e-22)**(-4/3) * cosmo_dic['h'] # to convert to Msun/h
        factor = 1 + (f_ax * (((1+gamma_1*M0/M)**(-gamma_2)) - 1.) * correction_array)
        #print('Gamma correction =', factor)
        B *= factor

    conc_param = B * (1 + func_z_formation(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma) )/(1+cosmo_dic['z']) 
    #print('Concentration - mass relation =', conc_param, M)
    return conc_param


#function for the normaliation factor in NFW profile
def func_for_norm_factor(x):
    """
    normalisation function for the NFW profile
    """
    return (- x/(1+x)) + np.log(1+x)


#density profile in k space (fourietrafo)
def func_dens_profile_kspace(M, k, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = False, axion_dic=None):
    """
    k, k_sigma units of h/Mpc, M in solar_mass/h and PS, PS_sigma in (Mpc/h)^3 
    NOTE: be carefull, we have two k's: k is the k, where the function is evaluated 
    and k_sigma is needed for for sigma(M, z) (the same is true for the PSs)
    NOTE: Omega_0 must match with chosen PS_sigma
    returns Fourier trafo of NFW profile (dimensionless) at k as given in my masterthesis eq. 4.20
    """
    #eta is a halo shape parameter introduced my Mead in https://arxiv.org/abs/2009.01858 in Tab2
    if eta_given == True:
        eta = hmcode_dic['eta']
        nu = func_nu(M, k_sigma, PS_sigma, Omega_0_sigma)
    else:
        eta = np.array([0.])
        nu = 1.
    
    R_vir = func_r_vir(cosmo_dic['z'], M, Omega_0, cosmo_dic['Omega_m_0'], cosmo_dic['Omega_w_0'])
    concentration = func_conc_param(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma, axion_dic=axion_dic) / (nu**eta)
    k_R_vir = np.outer(R_vir, k)
    a = np.outer(R_vir/concentration, k)
    
    def sin_integral(x):
        return scipy.special.sici(x)[0]
    def cos_integral(x):
        return scipy.special.sici(x)[1]
    
    summand1 = np.cos(a) * (cos_integral(a+k_R_vir) - cos_integral(a))
    summand2 = np.sin(a) * (sin_integral(a+k_R_vir) - sin_integral(a))
    summand3 = - np.sin(k_R_vir) / (a+k_R_vir)
    
    dens_profile_kspace = 1. / func_for_norm_factor(concentration)[:, None] * (summand1 + summand2 + summand3)
    
    #np.savetxt("dens_profile_kspace.txt",dens_profile_kspace)
    #np.savetxt("k.txt",k)
    #print("dens_profile_kspace",dens_profile_kspace)
    #print("shape dens_profile_kspace",np.shape(dens_profile_kspace))
    #print("shape of dens_profile_kspace", np.shape(dens_profile_kspace))
    return dens_profile_kspace


#delta_char for the NFW profile
def func_delta_char(M, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = False, axion_dic=None): 
    """
    k_sigma units of h/Mpc, M in solar_mass/h and PS_sigma in (Mpc/h)^3 
    returns NFW profile in h^2 * M_sun/Mpc^3 at k in h/Mpc
    as given in my masterthesis eq. 4.15 (left)
    """
    #eta is a halo shape parameter introduced my Mead in https://arxiv.org/abs/2009.01858 in Tab2
    if eta_given == True:
        eta = hmcode_dic['eta']
        nu = func_nu(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma)
    else:
        eta = np.array([0.]) 
        nu = 1.   
    Delta_vir = func_Delta_vir(cosmo_dic['z'], Omega_0, cosmo_dic['Omega_m_0'], cosmo_dic['Omega_w_0'])
    
    concentration = func_conc_param(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma, axion_dic=axion_dic) / (nu**eta)
    delta_char = func_rho_comp_0(Omega_0) * Delta_vir * concentration **3 / (3. * func_for_norm_factor(concentration))
    return delta_char


#density profile in real space
def NFW_profile(M, r, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = False, axion_dic=None):
    """
    r in units of Mpc/h, k_sigma in h/Mpc, M in solar_mass/h PS_sigma in (Mpc/h)^3
    returns NFW denisty profile in units of solar_mass/Mpc^3*h^2 at radius r
    """
    #eta is a halo shape parameter introduced my Mead in https://arxiv.org/abs/2009.01858 in Tab2
    if eta_given == True:
        eta = hmcode_dic['eta']
        nu = func_nu(M, k_sigma, PS_sigma, Omega_0)
    else:
        eta = np.array([0.]) 
        nu = 1.
        
    concentration = func_conc_param(M, k_sigma, PS_sigma, cosmo_dic, Omega_0_sigma, axion_dic=axion_dic) / (nu**eta)
    normalisation = func_delta_char(M, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = eta_given, axion_dic=axion_dic)
    r_s = func_r_vir(cosmo_dic['z'], M, Omega_0, cosmo_dic['Omega_m_0'], cosmo_dic['Omega_w_0']) / concentration
    
    NFW_func = 1 /((r/r_s) * (1+r/r_s)**2)
    
    return normalisation * NFW_func

