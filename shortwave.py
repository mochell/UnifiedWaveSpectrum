import numpy as np

from .unified import c, u_star
from .longwave import L_PM

def F_m(k, k_m):
    """
    Short-wave side effect function.
    """
    return np.exp(-0.25 * ((k / k_m) - 1)**2)

def c_m(k_m=370.8, g=9.81, return_wavenumber=False):
    """
    Phase speed at the high-frequency spectral peak.

    Parameters:
    k_m (float): Secondary (gravity-capillary) peak wave number in the curvature spectrum. Default is 370.88 rad/m.
    g (float): Gravitational acceleration. Default is 9.81 m/s^2.
    """
    if return_wavenumber:
        return np.sqrt(2 * g / k_m), k_m
    else:
        return np.sqrt(2 * g / k_m)

def a_m_ustar(u_star, c_m):
    """
    Generalized Phillips-Kitaigorodskii equilibrium range parameter for short waves.

    Parameters:
    u_star (float): Friction velocity
    c_m (float): Phase speed at the high-frequency spectral peak
    """
    if u_star <= c_m:
        return 1e-2 * (1 + np.log(u_star / c_m))
    else:
        return 1e-2 * (1 + 3 * np.log(u_star / c_m))

def a_m(u_10, k_m=None, C_d = 0.0015):
    """
    Generalized Phillips-Kitaigorodskii equilibrium range parameter for short waves.

    Parameters:
    u_10 (float): Wind speed at 10 m height
    k_m (float): Secondary (gravity-capillary) peak wave number in the curvature spectrum. Default is None.
    c_m (float): Phase speed at the high-frequency spectral peak
    """
    if k_m is None:
        c_m_value = c_m()
    else:
        c_m_value = c_m(k_m)

    return a_m_ustar(u_star(u_10), c_m_value)

def S_h(k, a_m, k_p, k_m=None, return_F_m=False):
    """
    High-wavenumber elevation spectrum.

    Parameters:
    k (float): Wavenumber
    a_m (float): Generalized Phillips-Kitaigorodskii equilibrium range parameter for short waves.
    k_m (float): Secondary (gravity-capillary) peak wave number in the curvature spectrum. Default is None.
    F_m (float): Short-wave side effect function
    """
    if k_m is None:
        c_m_value, k_m = c_m(return_wavenumber=True)
    else: 
        c_m_value = c_m(k_m)

    c_value = c(k)
    F_m_value = F_m(k, k_m) * L_PM(k, k_p)

    if return_F_m:
        return 0.5 * a_m * ( c_m_value / (k**3 * c_value) ) * F_m_value, F_m_value
    else:
        #return 0.5 * a_m * ( c_m_value / (k**3 * c_value) ) * F_m_value
        return 0.5 * a_m  * ( c_m_value / c_value )  * F_m_value / k**3