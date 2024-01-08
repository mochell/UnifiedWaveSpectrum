import numpy as np
#import numba

from .unified import c, u_star
from .longwave import k_p
from numba import njit

#@njit(parallel=True)
def Delta(k, k_p, u_10, c_m = 0.23):
    """
    Unified spreading function.

    Parameters:
    k (float): wavenumber
    k_p (float): Primary peak wavenumber
    u_10 (float): Wind speed at 10 m height
    c_m (float): 2nd-peak phase speed. Default is 0.23 m/s.
    """
    return Delta_ustar(k, k_p, u_star(u_10), c_m=0.23)


#@njit(parallel=True)
def Delta_ustar(k, k_p, u_star, c_m=0.23):
    """
    Unified spreading function.

    Parameters:
    k (float): Wavenumber
    k_p (float): Primary peak wavenumber
    u_star (float): Friction velocity
    c_m (float): 2nd-peak phase speed. Default is 0.23 m/s.
    """
    rho = 1000.   # Density of water in kg/m^3
    S = 0.072     # Surface tension of water in N/m
    g = 9.81      # Gravitational acceleration in m/s^2

    #c_value = c(k)
    c_value = np.sqrt(g/k + S/rho*k)

    #c_p = c(k_p)
    c_p = np.sqrt(g/k_p + S/rho*k_p)

    k_m = np.sqrt(rho*g/S)
    c_m = np.sqrt(g/k_m + S/rho*k_m)

    a_0 = np.log(2) / 4  # 0.1732#  # minimum value from LH 1963
    a_p = 4
    a_m = 0.13 * u_star / c_m  # consistent with Cox and Munk 1954

    return np.tanh( a_0 + a_p * (c_value / c_p)**2.5 + a_m * (c_m / c_value)**2.5 )


def Phi(k, phi, k_p, u_10):
    """
    Directional distribution function.

    Parameters:
    k (float): Wavenumber
    phi (float): Directional angle
    Delta (float): Unified spreading function
    """
    D = Delta(k, k_p, u_10, c_m = 0.23)
    return 1 / (2 * np.pi) * (1 + D * np.cos(2 * phi))  
    #return np.where((phi > - np.pi* 3/8.) & (phi < np.pi * 3/8.),  1 / (2 * np.pi) * (1 + D * np.cos(2 * phi))   , 0 )
    #return np.where((phi > - np.pi/2.) & (phi < np.pi/2.),  1 / (2 * np.pi) * (1 + D * np.cos(2 * phi))   , 0 )