

from __future__ import division

import numpy as np
from scipy.integrate import quad, trapz
from scipy.constants import c
import scipy.special
import astropy.constants as const
from astropy.units.quantity import Quantity

from function_definitions import planck_function, eddington_temperature, E2


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
    c1 = 0.56419
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
            "This piecewise approximation breaks for large values of `a`!\n" 
            "Tom Rice has made the code throw an error if a > 5e-3.")
    
    a = damping_parameter

    if np.abs(u) <= 2:
        return hjerting_gaussian_approximation(a, u)
    else:
        return (hjerting_gaussian_approximation(a, u) + 
                hjerting_H1_approximation(a, u))

hjerting_function = np.vectorize(hjerting_piecewise_approximation)


def flux_with_line(continuum_optical_depth, frequency, effective_temperature,
                   relative_line_center_opacity, line_center_frequency=6e14, 
                   doppler_width=4e9, damping_parameter=1e-4):
    """
    Calculates the flux at `tau_c` noting an absorption line at `nu_0`.

    This is meant to approximate function_definitions.flux for frequencies
    far away from nu_0.

    """

    tau_c = continuum_optical_depth
    nu = frequency
    T_eff = effective_temperature
    zeta = relative_line_center_opacity
    a = damping_parameter

    u = u_function(nu, line_center_frequency, doppler_width)

    line_absorption_factor = (1 + zeta * hjerting_function(a, u))

    temp_at_t = lambda t: eddington_temperature(t, T_eff)

    flux = (
        2 * np.pi * quad(
            (lambda t: planck_function(temp_at_t(t), nu) * 
             E2(line_absorption_factor*(t - tau_c))),
            tau_c, np.inf)[0] -
        2 * np.pi * quad(
            (lambda t: planck_function(temp_at_t(t), nu) * 
             E2(line_absorption_factor*(tau_c - t))),
            0, tau_c)[0] )

    return flux
