'''
functions added to change to the halo model window function, following HMCODE2020 https://arxiv.org/pdf/2009.01858

this changes the window function (the fourier transform of the halo matter density profile) used in calculating the 1-halo and 2-halo term, following Eq. 25 of this paper

Note: HMCODE2020 and axionHMCODE (https://arxiv.org/pdf/2209.13445) use slightly different definitions of the halo density profile, differing by a factor of 1/M (picked up again in the definition of the 1-halo term)
'''
from cold_density_profile import *

def func_dens_profile_kspace_baryons(M, k, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, Tagn, eta_given = False, axion_dic=None):
    '''
    Normalised Fourier transform for NFW profile, including baryonic effects
    Equation (25) from Mead et al. (2021)
    Changes Wk to account for halo deformation, stars and lower gas in haloes
    Changes made to make compatable with axionHMCODE 
    '''
    print("These changes are working!")
    Wk = func_dens_profile_kspace(M, k, k_sigma, PS_sigma, cosmo_dic, hmcode_dic, Omega_0, Omega_0_sigma, eta_given = False, axion_dic=None)

    print("M[0]",M[0])
    print("M[99]",M[99])

    print("Wk[0,:] before changes",Wk[0,:])
    print("Wk[99,:] before changes",Wk[99,:])

    ### Apply baryons below
    
    feedback_params = _get_feedback_parameters(Tagn) # Calling feedback parameters, calibrated on the BAHAMAS simulations, based off of T_agn

    #fg = ((cosmo_dic['omega_b_0']/cosmo_dic['omega_m_0']-feedback_params['f0'])*(M/feedback_params['Mb0'])**2/(1.+(M/feedback_params['Mb0'])**2))
    
    #print("Now trying to print fg!", fg)
    print("printing func_rho_comp_0(Omega_0)", func_rho_comp_0(Omega_0))
    print("printing feedback_params['f0']", feedback_params['f0'])
    print("printing  feedback_params['f0']*M/(func_rho_comp_0(Omega_0))", feedback_params['f0']*M/(func_rho_comp_0(Omega_0)))
    print("M/func_rho_comp_0(Omega_0)",M/func_rho_comp_0(Omega_0))

    Mb=feedback_params['Mb0']
    Om_b=cosmo_dic['omega_b_0']
    Om_c=cosmo_dic['omega_d_0']
    Om_m=cosmo_dic['omega_m_0']
    fstar=feedback_params['f0']
    
    fg = (Om_b/Om_m-fstar)*(M/Mb)**2/(1.+(M/Mb)**2) # Gas content (Eq. 24 with beta=2)
    WK = ((Om_c/Om_m+fg)*Wk.T).T
    WK = ( WK.T + (fstar /( M/(func_rho_comp_0(Omega_0)) ) ) ).T
    
    #WK = (( cosmo_dic['omega_d_0']/cosmo_dic['omega_m_0'] + fg )*Wk.T ).T 
    #WK = (WK.T + feedback_params['f0']/(M/(func_rho_comp_0(Omega_0)))).T
    #WK = (WK.T + feedback_params['f0']).T

    
    print("Now trying to print new Wk!",WK)
    print("Now trying to print new WK[0,:]!",WK[0,:])
    print("Now trying to print new WK[99,:]!",WK[99,:])

    print("shape of old Wk",np.shape(Wk))
    print("shape of new WK",np.shape(WK))

    Wk=WK

    return Wk

def _get_feedback_parameters(T_AGN):
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