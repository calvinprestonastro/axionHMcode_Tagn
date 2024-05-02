
print("check")
'''
functions added to change to the halo model window function, following HMCODE2020 https://arxiv.org/pdf/2009.01858

this changes the window function (the fourier transform of the halo matter density profile) used in calculating the 1-halo and 2-halo term, following Eq. 25 of this paper

Note: HMCODE2020 and axionHMCODE (https://arxiv.org/pdf/2209.13445) use slightly different definitions of the halo density profile, differing by a factor of 1/M (picked up again in the definition of the 1-halo term)
'''
from cold_density_profile import *

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
    return dens_profile_kspace


def func_dens_profile_kspace_baryons(M, k, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, Om_m, Om_c, Om_b, Mb, fstar, eta_given = False, axion_dic=None):
    '''
    Normalised Fourier transform for NFW profile, including baryonic effects
    Equation (25) from Mead et al. (2021)
    Calls _win_NFW
    Changes Wk to account for halo deformation, stars and lower gas in haloes
    Changes made to make compatable with axionHMCODE 
    '''
    Wk = func_dens_profile_kspace(M, k, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = False, axion_dic=None)
    fg = (Om_b/Om_m-fstar)*(M/Mb)**2/(1.+(M/Mb)**2) # Gas content (Eq. 24 with beta=2)
    M=np.array([M])
    print(len(M))
    for i in range(0, len(M) ):
        Wk[i,:] = (Om_c/Om_m+fg[i])*Wk[i,:]
        print(i)
    Wk += fstar
    return Wk

def _get_feedback_parameters(T_AGN:float) -> dict:
    '''
    Maps one-Param baryon feedback model from HMCode2020 to 6 baryonic parameters
    Uses parameters from Table 5 in Mead et al. (2021)
    This fit was obtained using the vanilla halo model! 
    If the hmcode tweaks are used, different values are likely needed.
    All values are reduced down to to one via best-fitting to the BAHAMAS model
    '''
    theta = np.log10(T_AGN/np.power(10, 7.8))
    params = {
        'B0': 3.44-0.496*theta,
        'Bz': -0.0671-0.0371*theta,
        'Mb0': np.power(10, 13.87+1.81*theta), # [Msun/h]
        'Mbz': -0.108+0.195*theta,
        'f0': (2.01-0.3*theta)*1e-2,
        'fz': 0.409+0.0224*theta,
    }
    return params

'''
def _get_feedback_suppression(k:np.array, zs:np.array, CAMB_results:camb.CAMBdata, T_AGN:float, 
                              Mmin=1e0, Mmax=1e18, nM=256, verbose=False) -> np.ndarray:
'''
'''
    Calculates the ratio of the powerspectrum with baryonic effects to that of dark-matter-only
    Assumes the one-parameter T_AGN model from HMCode2020
    Warning: Since the fit for the baryonic effects was obtained with the vanilla halo model, 
    it is not safe to set tweaks=True below
    '''
'''
    Pk_gravity = power(k, zs, CAMB_results, T_AGN=None, Mmin=Mmin, Mmax=Mmax, nM=nM, 
                       tweaks=False, verbose=verbose)
    Pk_feedback = power(k, zs, CAMB_results, T_AGN=T_AGN, Mmin=Mmin, Mmax=Mmax, nM=nM, 
                        tweaks=False, verbose=verbose)
    return Pk_feedback/Pk_gravity

'''
