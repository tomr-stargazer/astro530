"""
Let's compute and plot things for Problem 3.

Problem 3 is about a "purely viscous" accretion disk in which
the similarity solution holds; nu is proportional to R at all times.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad

from astropy.units.quantity import Quantity
from astropy.constants import M_sun, G

from problem_1_set5 import sound_speed

def temperature_at_R(radius):
    """
    Returns the temperature in kelvin at radius R in AU.

    """

    R = Quantity(radius, unit='AU')

    T = 300 * (R.value)**(-1/2)

    return Quantity(T, unit='K')

def angular_velocity_at_R(radius, central_mass=1):
    """
    Returns the angular velocity Omega at radius R in AU.

    """

    R = Quantity(radius, unit='AU')
    M = central_mass*M_sun

    omega = (G * M / R**3)**(1/2)

    return omega.to('s-1')
    

def viscosity_at_R(radius, alpha=0.01):
    """
    Returns the viscosity nu at radius R.
    
    """
    
    R = Quantity(radius, unit='AU')
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

    t = Quantity(time, unit='yr')

    t_V = (t * (3*viscosity_at_R(R_1) / R_1**2)) + Quantity(1, '')

    return t_V.decompose()
    
def surface_density_at_R_t(radius, time, normalization_constant=1,
                           scaling_radius=10):
    """
    Returns the surface density at radius R in AU.

    """

    R = Quantity(radius, 'AU')
    t = Quantity(time, 'yr')
    C = normalization_constant
    R_1 = Quantity(scaling_radius, unit='AU')
    nu_R1 = viscosity_at_R(R_1)

    t_V = t_V_at_t(time)

    return C * (t_V)**(-3/2)*R_1 / (3*np.pi*R*nu_R1) * np.exp(-R*R_1/t_V)

def Mdot_at_t_R(time, radius, normalization_constant=1,
                scaling_radius=10):
    """
    Returns the accretion rate Mdot at time t and radius R.

    """

    C = normalization_constant
    t = Quantity(time, 'yr')
    R = Quantity(radius, 'AU')
    R_1 = Quantity(scaling_radius, 'AU')

    t_V = t_V_at_t(t)

    Mdot = -C * (t_V)**(-3/2) * (np.exp(-R/(R_1 * t_V)) *
                                 (Quantity(1,'') - 2*R/(R_1 * t_V)))

    return Mdot

def M_d_at_t_in_R(time, radius, normalization_constant=1):
    """
    Returns the disk mass M_d contained in radius R at time t.
    
    """

    t = time
    R = radius
    C = normalization_constant

    # integrate 2 pi R Sigma from zero to R
    r_times_Sigma_at_t = lambda r: r*surface_density_at_R_t(r, t, C)

    return 2*np.pi* quad( r_times_Sigma_at_t, 0, R)[0]
    

def problem_3b(disk_mass=0.1):
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

    t_list = [0.5e6, 1e6, 3e6]
    r_array = np.logspace(-1, 2.25, 20)

    # calculate the normalization dudes
    C_list = []
    for t in t_list:

        C = disk_mass / M_d_at_t_in_R(t, r_array.max())

        C_list.append( C )

    fig_1 = plt.figure()
    ax1 = plt.gca()

    fig_2 = plt.figure()
    ax2 = plt.gca()

    fig_3 = plt.figure()
    ax3 = plt.gca()

    for t, C in zip(t_list, C_list):
        # Sigma(R)
        ax1.plot(r_array, surface_density_at_R_t(r_array, t, 
                                                 normalization_constant=C),
                 label=r'$%s \times 10^6$ yr' % (t/1e6))

        # M_d(R)
        ax2.plot(r_array, np.vectorize(M_d_at_t_in_R)(t, r_array, C),
                 label=r'$%s \times 10^6$ yr' % (t/1e6))
        
        # Mdot(R)
        ax3.plot(r_array, Mdot_at_t_R(t, r_array, normalization_constant=C),
                 label=r'$%s \times 10^6$ yr' % (t/1e6))

    ax1.legend(loc='lower left')
    ax2.legend(loc='lower right')
    ax3.legend(loc='upper left')
    
    ax1.set_ylabel(r"$\Sigma(R)$")
    ax2.set_ylabel(r"$M_d(R)$")
    ax3.set_ylabel("$\dot{M}(R)$")

    ax1.set_title("Surface density of disk")
    ax2.set_title("Mass of disk within radius R")
    ax3.set_title("Mass accretion rate of disk")

    for ax in [ax1, ax2, ax3]:
        ax.set_xlim(r_array.min(), r_array.max())
        ax.set_xlabel("Radius (AU)")        

    ax1.loglog()
    ax2.loglog()
    ax3.semilogx()

    ax1.set_ylim(1e-15, 1)
    
    plt.show()

def problem_3c(disk_mass=0.1):
    """
    Compute and plot the behavior of M_d(t) and Mdot(t, R->0).
    Give masses and mass accretion rates in units of M_sun 
    and M_sun yr^-1, respectively.

    """

    # M_d(t) is probably M_d(R_max) as a function of time.

    r_outer = 10**2.25

    t_array = np.logspace(0, 7, 100)
    C_list = []
    for t in t_array:

        C = disk_mass / M_d_at_t_in_R(t, r_outer)
        C_list.append( C )
        

    plt.figure()

    # M_d(t)
    plt.plot(t_array,
             np.vectorize(M_d_at_t_in_R)(t_array, r_outer, C_list))

    plt.figure()

    # Mdot(t, R->0)
    plt.plot(t_array, Mdot_at_t_R(t_array, 0))

    plt.show()
