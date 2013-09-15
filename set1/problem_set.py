"""
This is a script for the computational parts of Problem Set 1.

"""

from __future__ import division

import numpy as np
import astropy.constants as const
from astropy.units.quantity import Quantity


def planck_function(effective_temperature, frequency):
    """
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
    intensity : float
        The intensity $B_\nu(T_{eff})$ (also known as Spectral Radiance: 
        http://en.wikipedia.org/wiki/Radiance ) of the blackbody radiation
        at the given frequency.

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


def test_planck_function():

    # See: http://wolfr.am/1ghoEBo for Wolfram Alpha example result
    assert round(planck_function(6000, 200e12), 8+3) == 2.985e-8
    
    
