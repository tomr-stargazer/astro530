"""
Computes and plots the observed SED for a disk model from last homework.

"""

from __future__ import division

import numpy as np
from scipy.integrate import quad
from astropy.constants import c, R_sun
from astropy.units.quantity import Quantity

from astro530.set5.problem_3_set5 import (temperature_at_R, 
                                          surface_density_at_R_t, 
                                          M_d_at_t_in_R)
from astro530.set1.function_definitions import planck_function


def dust_opacity_at_nu(nu):
    """ Calculates dust opacity kappa_nu at wavelength nu. """

    kappa_nu = 0.1 * (nu / 1e12)

    return kappa_nu

R_outer = Quantity(10**2.25, 'AU')
R_inner = Quantity(0.04, 'AU')

t = Quantity(1e6, 'yr')

C = 0.1 / M_d_at_t_in_R(t, R_outer)

def disk_flux_at_R_nu(radius, frequency):
    """ Computes the disk spectral flux at radius R for frequency nu. """

    # Looks like B_nu (T(R)) * 1-exp(-tau)
    # for tau = sigma * kappa

    R = radius
    nu = frequency

    tau = surface_density_at_R_t(R, t, C) * dust_opacity_at_nu(nu)

    blackbody_flux = planck_function(temperature_at_R(R), nu).to('erg s-1 cm-2 Hz-1')

    disk_flux = blackbody_flux * (1 - np.exp(-tau))

    return disk_flux


def total_disk_flux_at_nu(nu):
    """ Integrates over the whole disk. """

    flux_at_R_times_R = lambda R: disk_flux_at_R_nu(R, nu)*R

    distance = Quantity(100, 'pc')

    pre_distanced_flux = Quantity(
        2*np.pi*quad(flux_at_R_times_R, R_inner, R_outer)[0],
        'erg s-1 cm-2 Hz-1 AU2')

    flux =  (pre_distanced_flux/distance**2)

    return flux.to('erg cm-2 Hz-1 s-1')

def flux_from_star(nu):

    blackbody_flux = planck_function(5000, nu).to('erg s-1 cm-2 Hz-1')

    radius = 2 * R_sun
    distance = Quantity(100, 'pc')

    flux = 2*np.pi*(radius**2/distance**2) * blackbody_flux

    return flux.to('erg cm-2 Hz-1 s-1')

def plot_SED_disk():

    wavelength_array = Quantity(np.logspace(-1, 2, 20)*3e-6, 'm')
    nu_array = c / wavelength_array

    SED_array = np.zeros_like(nu_array)
    
    for nu, i in zip(nu_array, range(len(nu_array))):

        SED_array[i] = nu*(total_disk_flux_at_nu(nu))

    return (wavelength_array, SED_array)

def plot_SED_star():

    wavelength_array = Quantity(np.logspace(-1, 2, 20)*3e-6, 'm')
    nu_array = c / wavelength_array

    SED_array = np.zeros_like(nu_array)
    
    for nu, i in zip(nu_array, range(len(nu_array))):

        SED_array[i] = nu*(flux_from_star(nu))

    return (wavelength_array, SED_array)
        
