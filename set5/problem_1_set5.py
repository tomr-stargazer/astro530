"""
This is code to do the computations & plotting for 1c.

Problem 1 is about a rotating, collapsing protostellar cloud of 
one solar mass.

"""

from __future__ import division

import numpy as np

from astropy.constants import k_B, m_p, G, M_sun
from astropy.units.quantity import Quantity

def sound_speed(temperature, mean_molecular_weight=2.3):
    """
    Computes the sound speed $c_s$ given an input temperature.

    c_s = \sqrt( kT / \mu m_H )

    Parameters
    ----------
    temperature : float
        Temperature of the gas in units kelvin.
    mean_molecular_weight : float, optional
        Mean molecular weight (\mu) of the gas.
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
    
def centrifugal_radius(time, temperature, angular_velocity_factor):
    """
    Computes the centrifugal radius $r_c$ of the infalling envelope.

    r_c(t) = c_s t^3 / 16 * \Omega_0^2

    Parameters
    ----------
    time : float
        Time coordinate (t) in seconds.
    temperature : float
        Temperature (T) of the envelope in kelvin.
    angular_velocity_factor : float
        Multiplier on the breakup angular frequency (\Omega_b) 
        of the rotation.

    Returns
    -------
    r_c : float
        Centrifugal radius in meters.

    """

    t = Quantity(time, unit='s')
    T = temperature

    c_s = sound_speed(T)

    r_o = 0.5 * G * M_sun / c_s**2
    angular_velocity = angular_velocity_factor * (G * M_sun / r_o**3)**(1/2)
    
    omega = Quantity(angular_velocity, unit='s-1')

    r_c = c_s * t**3 / 16 * omega**2

    return r_c.to('m')
