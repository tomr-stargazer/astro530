

from __future__ import division

import numpy as np


def u_function(frequency, central_frequency, doppler_width):
    """
    u is a parameterized form of frequency that is used in the Hjerting function.

    It's differenced from the line's central frequency and then divided
    by the local Doppler width of the line.

    Make sure you provide all the inputs in the same units (e.g. Hz).

    Parameters
    ----------
    frequency : float
        The frequency $\nu$ you want to calculate u for.
    central_frequency : float
        The central frequency $\nu_0$ of the spectral line we're looking at.
    doppler_width : float
        The doppler width $\Delta \nu_D$ at the local level, in frequency units.

    Returns
    -------
    u : float
        The `u` parameter that is used in the Hjerting function.
        It's unitless.

    """

    nu = frequency
    nu_0 = central_frequency
    delta_nu = doppler_width

    return (nu - nu_0) / delta_nu

def hjerting_gaussian_approximation(damping_parameter, u):
    """
    Calculates the H(0, u) Hjerting Gaussian approximation.

    """

    return np.exp(-u**2)


def hjerting_H1_approximation(damping_parameter, u):
    """
    Calculates the H_1 polynomial Hjerting function approximation.

    """

    a = damping_parameter

    # Two constants from the homework
    c1 = 0.56491
    c2 = 0.846

    return a * (c1/u**2 + c2/u**4)
    

def hjerting_piecewise_approximation(damping_parameter, u):
    """
    Approximates the Hjerting function correctly, everywhere.

    At small |u|, it uses the Gaussian approximation; 
    at larger |u| it uses the H_1 approximation plus the Gaussian approximation.

    We are doing a piecewise approximation to avoid the behavior of
    H_1 near u=0.

    """

    if damping_parameter > 5e-3:
        raise ValueError(
            "This piecewise approximation breaks for large values of `a`!")
    
    a = damping_parameter

    if np.abs(u) <= 2:
        return hjerting_gaussian_approximation(a, u)
    else:
        return (hjerting_gaussian_approximation(a, u) + 
                hjerting_H1_approximation(a, u))
