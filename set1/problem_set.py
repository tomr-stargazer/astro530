"""
This is a script for the computational parts of Problem Set 1.

"""

from __future__ import division

import numpy as np
#import scipy.constants
import astropy.constants as const


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

    h = const.h.to('J s')#.value # Planck's constant h, in Joules * seconds.
    c = const.c.to('m/s')#.value # Speed of light c, in meters / seconds.
    k_B = const.k_B.to('J/K')#.value # Boltzmann constant k_B, in Joules / kelvin.

    intensity = ((2 * h * (frequency)**3) / c**2 / 
                 (np.exp( (h * frequency) / (k_B * effective_temperature)) - 1)
                 )

    return intensity


def test_planck_function():

    assert 1 == 1 #placeholder

    #    assert planck_function(6000, 2e12) == 
    
    
