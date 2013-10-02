"""
Code relating to problem 1 in Astro 530 HW #2.

Problem 1 is about the Saha equation and ions of various atoms
in three different stellar atmospheres.

The Saha equation:
https://en.wikipedia.org/wiki/Saha_ionization_equation

but I have followed Carroll and Ostlie's choice to place 
the electron number density $n_e$ on the denominator of the right 
side of the equation, rather than the numerator of the left side.

Thus:

..math::
\frac{n_{i+1}}{n_i} = \frac{2}{n_e\Lambda^{3}}\frac{g_{i+1}}{g_i}\exp\left[-\frac{(\epsilon_{i+1}-\epsilon_i)}{k_BT}\right]

"""

from __future__ import division

import numpy as np

import astropy.constants as const
from astropy.units.quantity import Quantity


Model_1 = {'T_eff': 4500,
           'log g': 4.0,
           'T': 4469.7,
           'n_e': 4.75e12,
           'n_H': 6.6e16}

Model_2 = {'T_eff': 5750,
           'log g': 4.0,
           'T': 5609.3,
           'n_e': 9.55e12,
           'n_H': 6.0e16}

Model_3 = {'T_eff': 15000,
           'log g': 4.0,
           'T': 14630,
           'n_e': 3.17e14,
           'n_H': 3.67e14}


def saha_equation(degeneracy_of_lower_state,
                  degeneracy_of_upper_state,
                  ionization_energy_of_lower_state,
                  ionization_energy_of_upper_state,
                  temperature,
                  electron_number_density,
                  units_of_energy='eV'):
    """
    Calculates the Saha equation for a given set of parameters.

    Returns
    -------
    saha_ratio : float
        The ratio n_{i+1} / n_i, the output of the Saha equation.

    """

    g_i = degeneracy_of_lower_state
    g_ip1 = degeneracy_of_upper_state

    epsilon_i = Quantity(ionization_energy_of_lower_state, units_of_energy)
    epsilon_ip1 = Quantity(ionization_energy_of_upper_state, units_of_energy)

    n_e = Quantity(electron_number_density, '1 / (m3)')
    T = Quantity(temperature, 'K')

    h = const.h.to(units_of_energy+' s')
    k_B = const.k_B.to(units_of_energy+' / (K)')

    # "Lambda" is just a name for this bundle of parameters and constants.
    lambda_squared = (h**2 / (2 * np.pi * const.m_e * k_B * T)).to('m2')

    saha_ratio = ( 2 / lambda_squared**(3/2) *
                   g_ip1 / g_i *
                   np.exp(-(epsilon_ip1 - epsilon_i) /(k_B * T))
                   ) / n_e

    return saha_ratio
