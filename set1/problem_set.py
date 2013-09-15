"""
This is a script for the computational parts of Problem Set 1.

"""

from __future__ import division

import scipy.constants
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
    effective_temperature : float
        The effective temperature $T_{eff}$ of the blackbody 
        under consideration.
    frequency : float
        The frequency $\nu$ of the light being considered.

    Returns
    -------
    intensity : float
        The intensity $B_\nu(T_{eff})$ (also known as Spectral Radiance: 
        http://en.wikipedia.org/wiki/Radiance ) of the blackbody radiation
        at the given frequency.

    """

    intensity = 2 * const.h * (frequency)**3 

    return intensity


def test_planck_function():

    assert 1 == 1 #placeholder
    
