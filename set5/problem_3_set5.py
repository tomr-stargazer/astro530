"""
Let's compute and plot things for Problem 3.

Problem 3 is about a "purely viscous" accretion disk in which
the similarity solution holds; nu is proportional to R at all times.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from problem_1_set5 import sound_speed

def temperature_at_R(radius):
    """
    Returns the temperature in kelvin at radius R in AU.

    """

    R = radius

    return 300 * (R)^(-1/2)

def surface_density_at_R(radius):
    """
    Returns the surface density at radius R in AU.

    """

    pass

def t_V_at_t(time):
    """
    Returns the quantity t_V at time t.

    """

    pass

def Mdot_at_t_R(time, radius):
    """
    Returns the accretion rate Mdot at time t and radius R.

    """

    pass

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

    
