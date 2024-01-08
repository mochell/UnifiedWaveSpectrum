import numpy as np
#import numba

from .unified import c
# Inverse wave age parameter
def Omega_1D_simple(u10, c_p):
    """
    This function calculates the inverse-wave age parameter for a given phase speed at spectral peak (c_p).
    
    Parameters:
    u10 (float): Wind speed
    c_p (float): Phase speed at spectral peak

    Returns:
    float: The inverse-wave age parameter
    """
    return u10 / c_p

def Omega_1D(X, method='Elfouhaily'):
    """
    This function calculates the inverse-wave age parameter for a given method.
    
    Parameters:
    X (float): Non-dimensional fetch
    method (str): Method to calculate Omega. Default is 'Elfouhaily'.

    Returns:
    float: The inverse-wave age parameter
    """
    if method == 'Elfouhaily':
        X_0 = 2.2 * 10**4
        return 0.84 * np.tanh((X / X_0)**0.4)**-0.75
    elif method == 'Donelan':
        C_1 = 1.81 * 10**-2 * X**0.4
        return 0.91 * np.tanh(C_1)**-0.75
    elif method == 'Hasselmann':
        return 22 * X**-0.33
    elif method == 'Donelan2':
        return 11.6 * X**-0.23
    else:
        raise ValueError("Invalid method. Expected one of: 'Elfouhaily', 'Donelan', 'Hasselmann', 'Donelan2'")


def Omega_2D(u_10, c_p, method='PiCLES', theta=None):
    """
    This function calculates the inverse-wave age parameter for a g iven method.
    
    Parameters:
    u_10 (float or numpy.ndarray): Wind speed at 10m height
    c_p (float or numpy.ndarray): Phase speed at spectral peak
    method (str): Method to calculate Omega. Default is 'PiCLES'.
    theta (float, optional): Angle between the wind and peak wave direction. Required if u_10 and c_p are scalars.

    Returns:
    float: The inverse-wave age parameter
    """

    if method == 'directional':
        if isinstance(u_10, np.ndarray) and isinstance(c_p, np.ndarray):
            # If u_10 and c_p are vectors, calculate theta
            theta = np.arccos(np.dot(u_10, c_p) / (np.linalg.norm(u_10) * np.linalg.norm(c_p)))        
            alpha_p = np.linalg.norm(u_10) / np.linalg.norm(c_p)
            return np.where((theta > - np.pi/2.) & (theta < np.pi/2.),  alpha_p * np.cos(theta) , 0 )

        elif theta is None and (not isinstance(u_10, np.ndarray) or not isinstance(c_p, np.ndarray)):
            raise ValueError("Theta is required for the 'directional' method when u_10 and c_p are scalars.")
        else:
            alpha_p = u_10 / c_p
            return np.where((theta > - np.pi/2.) & (theta < np.pi/2.),  alpha_p * np.cos(theta) , 0 )

    elif method == 'PiCLES':
        if not isinstance(u_10, np.ndarray) or not isinstance(c_p, np.ndarray):
            raise ValueError("u_10 and c_p are required to be vectors for PiCLES method.")
        omega = (c_p[0] * u_10[0] - c_p[1] * u_10[1]) / (2 * np.linalg.norm(c_p)**2 )
        return np.where( omega > 0, omega, 0 )

    else:
        raise ValueError("Invalid method. Expected one of: 'directional', 'PiCLES'")

# %% Fetch-law k_p predictors:

def k_p(Omega_1D, u_10):
    """
    This function calculates the spectral peak wave number for a given inverse-wave age parameter.
    
    Parameters:
    Omega_1D (float): Inverse-wave age parameter
    u_10 (float): Wind speed at 10m height

    Returns:
    float: The spectral peak wave number
    """
    k_0 = 9.81  / u_10**2
    return Omega_1D**2  * k_0 


# %% Low-wave curvature spectrum

def S_l(k, k_p, Omega, return_F_p=False):
    """
    This function calculates the long-wave elevation spectrum for a given wave number (k).
    
    Parameters:
    k (numpy.ndarray): Wave number vector
    k_p (float): Peak wave number
    Omega (float): Inverse-wave age parameter

    Returns:
    numpy.ndarray: The long-wave elevation spectrum at the given parameters
    """
    g = 9.81
    a_p = cal_a_p(Omega)
    c_p = c(k_p)
    F_p = L_PM(k, k_p) * J_p(k, k_p, Omega) * np.exp(-Omega / np.sqrt(10) * (np.sqrt(k / k_p) - 1))
    if return_F_p:
        return (a_p / (2 * k**3)) * (c_p / c(k)) * F_p, F_p
    else:
        return (a_p / (2 * k**3)) * (c_p / c(k)) * F_p


def cal_a_p(Omega):
    """
    This function calculates the generalized Phillips-Kitaigorodskii equilibrium range parameter for long waves.
    
    Parameters:
    Omega (float): Inverse-wave age parameter

    Returns:
    float: The generalized Phillips-Kitaigorodskii equilibrium range parameter for long waves
    """
    return 6 * 10**-3 * np.sqrt(Omega)

# Pierson-Moskowitz spectrum
def L_PM(k, k_p):
    """
    Pierson-Moskowitz spectrum

    Parameters:
    k (numpy.ndarray): Wave number vector
    k_p (float): Peak wave number

    Returns:
    numpy.ndarray: The Pierson-Moskowitz spectrum at the given parameters
    """
    return np.exp(-5 / 4 * (k_p / k)**2)

# %% JONSWAP peak functions

def J_p(k, k_p, Omega):
    """
    JONSWAP peak function

    Parameters:
    k (numpy.ndarray): Wave number vector
    k_p (float): Peak wave number
    Omega (float): Inverse-wave age parameter

    Returns:
    numpy.ndarray: The JONSWAP peak function at the given parameters
    """
    return gamma(Omega)**Gamma(k, k_p, Omega)

def gamma(Omega):
    """
    Parameter for the spectral peak base for the JONSWAP peak function.
    """

    if 0.84 < Omega < 1:
        return 1.7
    elif 1 < Omega < 5:
        return 1.7 + 6 * np.log(Omega)
    else:
        return np.nan  # Return NaN for values of Omega outside the specified range

def Gamma(k, k_p, Omega):
    """
    Parameter for the spectral peak width for the JONSWAP peak function.
    """
    sigma = 0.08 * (1 + 4 / Omega**3)
    return np.exp(-((np.sqrt(k / k_p) - 1)**2) / (2 * sigma**2))


# %% High-wave curvature spectrum



