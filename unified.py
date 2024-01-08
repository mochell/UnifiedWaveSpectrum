
import numpy as np
# %% Unified 2D spectrum
#@numba.jit(nopython=True)

def Psi(k, varphi, S_l, S_h, Phi):
    """
    This function calculates the unified Displacement spectrum for given wave number (k) and angle (varphi) vectors.
    
    Parameters:
    k (numpy.ndarray): Wave number vector
    varphi (numpy.ndarray): Angle vector in radians

    S_l (float): Low-wave curvature spectrum
    S_h (float): Short-wave curvature spectrum
    Delta (float): Unified spreading function

    Returns:
    numpy.ndarray: The 2D matrix of the unified spectrum at the given parameters
    """
    #k, varphi = np.meshgrid(k, varphi)
    #return (1 / (2 * np.pi * k)) * (S_l + S_h) * (1 + Delta * np.cos(2 * varphi))
    return (S_l + S_h) * Phi # (1 + Delta * np.cos(2 * varphi))

# %% General Deep water wave functions

# #@numba.jit(nopython=True)
# def c(k):
#     """
#     phase speed for a given wave number (k).
    
#     Parameters:
#     k (numpy.ndarray): Wave number vector
#     g (float): Gravitational acceleration

#     Returns:
#     numpy.ndarray: The phase speed at the given wave number
#     """
#     g= 9.81
#     return np.sqrt(g / k)

# check c definitions with this:
#https://github.com/pakodekker/oceansar/blob/ffc6f13be34fef8ab9e7570a005b130e6462688e/oceansar/spec/elfouhaily.py#L43
def c(k, k_m = 370.88):
    """
    phase speed for a given wave number (k).
    
    Parameters:
    k (numpy.ndarray): Wave number vector
    k_m (float): Secondary (gravity-capillary) peak wave number in the curvature spectrum. Default is 370.88 rad/m.

    Returns:
    numpy.ndarray: The phase speed at the given wave number
    """
    g= 9.81
    return np.sqrt( (g / k) * (1  + k/k_m))




def C_d_neutral(u_10):
    """
    Neutral drag coefficient.

    Parameters:
    u_10 (float): Wind speed at 10 m height
    """
    return (0.8 + 0.065 * u_10) * 1e-3

def u_star(u_10):
    """
    Friction velocity.

    Parameters:
    u_10 (float): Wind speed at 10 m height
    """
    return np.sqrt( C_d_neutral(u_10) ) * u_10