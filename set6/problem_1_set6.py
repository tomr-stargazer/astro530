"""
Computes and plots the observed SED for a disk model from last homework.

"""

from __future__ import division

from problem_3_set5 import (temperature_at_R, surface_density_at_R_t, 
                            M_d_at_t_in_R)


def dust_opacity_at_nu(nu):
    """ Calculates dust opacity kappa_nu at wavelength nu. """

    kappa_nu = 0.1 * (nu / 1e12)

    return kappa_nu

