"""
This is a script for the computational parts of Problem Set 1.

"""

from __future__ import division


def planck_function(effective_temperature, frequency):
    """
    Computes the Planck's Law intensity per frequency at T_eff.

    T_eff : `effective_temperature`
    nu : `frequency`

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

    
