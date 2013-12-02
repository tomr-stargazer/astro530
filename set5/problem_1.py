"""
This is code to do the computations & plotting for 1c.

Problem 1 is about a rotating, collapsing protostellar cloud of 
one solar mass.

"""

from __future__ import division

import numpy as np

from astropy.constants import k_B, m_p
from astropy.units.quantity import Quantity

def sound_speed(temperature, mean_molecular_weight=1.4):
    """
    Computes the sound speed $c_s$ given an input temperature.

    c_s = \sqrt( kT / \mu m_H )

    Parameters
    ----------
    temperature : float
        Temperature in units kelvin of the gas.
    mean_molecular_weight : float, optional
        Mean molecular weight (\mu) of the gas under consideration.
        Default is 1.4, corresponding to some mix of H and He.

    Returns
    -------
    c_s : float
        The sound speed corresponding to the given temperature.

    """

    T = Quantity(temperature, unit='K')
    mu = mean_molecular_weight

    c_s = (k_B * T / (mu * m_p) )**(1/2)
    
    return c_s.to('m/s')
    
