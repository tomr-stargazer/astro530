"""
This is a script for the computational parts of Problem Set 1.

Some tips given here: http://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html

I'm trying out the Quantity / constants modules from astropy, 
with the unit handling and everything. It's going well so far.

"""

from __future__ import division

import numpy as np
from scipy.integrate import quad
import astropy.constants as const
from astropy.units.quantity import Quantity


def planck_function(effective_temperature, frequency):
    r"""
    Computes the Planck's Law intensity per frequency at T_eff.

    T_eff : `effective_temperature`
    nu : `frequency`

    Planck's Law:
    B_\nu(T_{eff}) = \frac{2 h \nu^3}{ c^2 } ( \exp(h \nu / k_B T_{eff}) - 1)^{-1} 
    (see: http://en.wikipedia.org/wiki/Planck's_law)

    Parameters
    ----------
    effective_temperature : float (in units Kelvin)
        The effective temperature $T_{eff}$ of the blackbody 
        under consideration. 
    frequency : float (in units seconds^{-1})
        The frequency $\nu$ of the light being considered.

    Returns
    -------
    intensity : astropy.units.quantity.Quantity
        The intensity $B_\nu(T_{eff})$ (also known as Spectral Radiance: 
        http://en.wikipedia.org/wiki/Radiance ) of the blackbody radiation
        at the given frequency. Given in units of Hz^3 J s^3 / (m^2), which is
        identical to units of W / (m^2 Hz sr) .

    """

    # Converting user-supplied arguments into Quantity values.
    Teff = Quantity(effective_temperature, 'K')
    nu = Quantity(frequency, 'Hz')
    
    h = const.h.to('J s') # Planck's constant h, in Joules * seconds.
    c = const.c.to('m/s') # Speed of light c, in meters / seconds.
    k_B = const.k_B.to('J/K') # Boltzmann constant k_B, in Joules / kelvin.

    intensity = ((2 * h * (nu)**3) / c**2 / 
                 (np.exp( (h * nu) / (k_B * Teff)) - 1)
                 )

    return intensity


def E2(x):
    """ Second exponential integral. See Lee's notes. """
    
    e2 = x * quad( lambda t: t**(-2) * np.exp(-t), x, np.inf)[0]
    return e2


def eddington_temperature(optical_depth, effective_temperature):
    """ 
    Uses the Eddington approximation to give a local temperature
    as a function of observer's optical depth and the star's T_eff.
    
    \tau : `optical_depth`
    T_{eff} : `effective_temperature`

    T^4 = (3/4) T_{eff}^4 (\tau + 2/3)
    
    """
    
    tau = optical_depth
    T_eff = effective_temperature
    
    T4 = (3/4) * T_eff**4 * (tau + 2/3)
    
    return T4**(1/4)


def flux(optical_depth, frequency, effective_temperature):
    r"""
    Computes the Flux per frequency at an optical depth for T_eff.

    T_eff : `effective_temperature`
    nu : `frequency`
    tau : `optical_depth`

    Flux: 
    .. math::
    F_\nu(\tau_\nu) = 2 \pi \int_{\tau_\nu}^{\infty} S_\nu(t_\nu) E_2(t_\nu - \tau_\nu) d t_\nu - 2 \pi \int_0^{\tau_\nu} S_\nu(t_\nu) E_2(\tau_\nu - t_\nu)dt_\nu

    in this case, S_\nu = B_\nu (Planck's function, as shown below):
    
    Planck's Law:
    .. math::
    B_\nu(T_{eff}) = \frac{2 h \nu^3}{ c^2 } ( \exp(h \nu / k_B T_{eff}) - 1)^{-1} 
    (see: http://en.wikipedia.org/wiki/Planck's_law)

    and E_2 is given by:
    .. math::
    E_2 = x \int_x^\infty y^{-2} e^{-y} dy
    

    Parameters
    ----------
    optical_depth : float (unitless)
        The optical depth $\tau$ through which an observer views the flux.
    frequency : float (in units seconds^{-1})
        The frequency $\nu$ of the light being considered.
    effective_temperature : float (in units Kelvin)
        The effective temperature $T_{eff}$ of the blackbody 
        under consideration. 

    Returns
    -------
    flux : float
        The flux $F_\nu(\tau)$ observed at the given optical depth inside
        the star with T_eff at the given frequency.

    """

    tau = optical_depth
    nu = frequency
    T_eff = effective_temperature

    temp_at_t = lambda t: eddington_temperature(t, T_eff)

    flux = (
        2 * np.pi * quad(
            (lambda t: planck_function(temp_at_t(t), nu) * E2(t - tau)),
            tau, np.inf)[0] -
        2 * np.pi * quad(
            (lambda t: planck_function(temp_at_t(t), nu) * E2(tau - t)),
            0, tau)[0] )

    return flux


def test_planck_function():

    # See: http://wolfr.am/1ghoEBo for Wolfram Alpha example result
    assert round(planck_function(6000, 200e12), 8+3) == 2.985e-8
    

def test_eddington_temperature():

    # Practically by definition, T(tau = 2/3) == T_eff
    assert eddington_temperature(2/3, 6000) == 6000
    assert eddington_temperature(2/3, 10) == 10

    assert eddington_temperature(1, 6000) > 6000
    assert eddington_temperature(1/2, 6000) < 6000
