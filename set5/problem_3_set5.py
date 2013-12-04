"""
Let's compute and plot things for Problem 3.

Problem 3 is about a "purely viscous" accretion disk in which
the similarity solution holds; nu is proportional to R at all times.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from astropy.units.quantity import Quantity
from astropy.constants import M_sun, G

from problem_1_set5 import sound_speed

def temperature_at_R(radius):
    """
    Returns the temperature in kelvin at radius R in AU.

    """

    R = radius

    return 300 * (R)**(-1/2)

def angular_velocity_at_R(radius, central_mass=1):
    """
    Returns the angular velocity at radius R in AU.

    """

    R = Quantity(radius, unit='AU')
    M = central_mass*M_sun

    omega = (G * M / R**3)**(1/2)

    return omega
    

def viscosity_at_R(radius, alpha=0.01):
    """
    Returns the viscosity nu at radius R.
    
    """
    
    R = radius
    T = temperature_at_R(R)
    c_s = sound_speed(T)
    omega = angular_velocity_at_R(R)

    nu = alpha*c_s**2 / omega

    return nu

def t_V_at_t(time, scaling_radius=10):
    """
    Returns the quantity t_V at time t.

    """

    R_1 = Quantity(scaling_radius, unit='AU')

    t = time

    t_V = t * 3*viscosity_at_R(R_1) / R_1**2 + 1

    return t_V
    
def surface_density_at_R_t(radius, time, normalization_constant=1,
                           scaling_radius=10):
    """
    Returns the surface density at radius R in AU.

    """

    R = radius
    t = time
    C = normalization_constant
    R_1 = Quantity(scaling_radius, unit='AU')
    nu_R1 = viscosity_at_R(scaling_radius)

    t_V = t_V_at_t(time)

    return C * (t_V)**(-3/2)*R_1 / (3*np.pi*R*nu_R1) * np.exp(-R*R_1/t_V)

def Mdot_at_t_R(time, radius, normalization_constant=1,
                scaling_radius=10):
    """
    Returns the accretion rate Mdot at time t and radius R.

    """

    C = normalization_constant
    t = time
    R = radius
    R_1 = scaling_radius

    t_V = t_V_at_t(time)
    
    Mdot = C * t_V * np.exp(-R/(R_1 * t_V))*(1 - 2*R/(R_1 * t_V))

def problem_3b():
    """
    Consider a 1 solar mass central star with R_1 = 10 AU and
    a disk mass of 0.1 M_sun at the initial time t=0.

    Assume the viscosity can be parameterized by nu = alpha c_s^2/Omega
    where alpha = 0.01 and c_s is the sound speed for a gas of mean
    molecular weight 2.3.
    
    Furthermore, assume that the disk is vertically isothermal and 
    its temperature is dominated by irradiation from the central star,
    such that T(R) = 300 (R/1AU)^(-1/2) K.

    Compute and plot the surface density Sigma(R), the mass 
    accretion rate Mdot(R), and the mass interior to R, M_d(R), 
    at t=0.5, 1, and 3e6 years. Make the plots of Sigma and M_d in log-log
    plots for R measured in AU. You will have to plot M_d linearly, so
    rescale at each time.
    
    """


    
    
