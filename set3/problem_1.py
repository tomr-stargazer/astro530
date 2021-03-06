"""
Does computational things for problem 1 in HW #3 in astro530.

Borrows from the code I wrote in HW #1, so I'm really glad I spent 
the time to make it high quality and tested, etc back then.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
from astropy.units.quantity import Quantity
from scipy.integrate import trapz

from function_definitions import flux
from function_definitions_hw3 import *

line_center_wavelength = Quantity(5000, 'Angstrom')
line_center_frequency = const.c/line_center_wavelength

def problem_1a():
    """
    Calculate the emergent flux at 5000 Angstroms (0.5 microns) for
    a gray atmosphere, using the Eddington approximation, 
    and an effective temperature of 5780 K, with the source function
    as the Planck function.

    """    
    result = flux(0, line_center_frequency, 5780)

    # these units are not DIRECTLY tested anywhere. I currently trust them. 
    print result, "in units of W / (m^2 Hz sr)"

    return result

def problem_1b():
    """
    Plot H(a, u) from u=0 to u=10, for a=1e-3 and a=1e-4, in loglog space.

    """

    fig = plt.figure()

    u_array = np.arange(0,10, 0.01)

    hv = np.vectorize(hjerting_piecewise_approximation)

    plt.plot(u_array, hv(1e-4, u_array), label=r"$a = 10^{-4}$" )
    plt.plot(u_array, hv(1e-3, u_array), '--', label=r"$a = 10^{-3}$" )

    plt.legend()

    plt.xlabel("u")
    plt.ylabel(r"$\phi_\nu$", rotation='horizontal')
    plt.title("HW #3, Problem 1b: Normalized Line Profile Function. Tom Rice.")
    
    plt.loglog()

    plt.show()

    return fig

def problem_1c(sampling=100):
    """
    Plot the line profile for a diversity of different zetas.

    Do it once for zero optical depth and once for high optical depth.

    """

    fig1 = plt.figure()
    fig2 = plt.figure()

    zeta_list = [0.1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
    a = 1e-4
    T_eff = 5780

    delta_nu = 4e9 # corresponding to v_thermal = 2 km/s    

    vflux_with_line = np.vectorize(flux_with_line)

    nu_array = (6+5*np.logspace(-10, -3, sampling))*1e14
    ratio_array = np.zeros_like(nu_array)

    equivalent_width_array = np.zeros_like(zeta_list)
    equivalent_width_array_lambda = np.zeros_like(zeta_list)    

    tau_c_list = [0, 10]

    for tau_c, fig in zip(tau_c_list, (fig1, fig2)):

        ax = fig.add_subplot(1,1,1)

        for j, zeta in zip(range(len(zeta_list)), zeta_list):

            for i in range(len(ratio_array)):
                ratio_array[i] = (flux_with_line(
                    tau_c, nu_array[i], T_eff, zeta) /
                    flux(tau_c, np.median(nu_array), T_eff) )

            ax.plot(nu_array, ratio_array, label=r"$\zeta = 10^{%d}$" % 
                    np.log10(zeta))

            equivalent_width_array[j] = trapz(1-ratio_array, nu_array)*2
            equivalent_width_array_lambda[j] = trapz(
                1-ratio_array, -(const.c / nu_array))*2

            print "done with zeta = %f" % zeta

        ax.set_title(r"HW #3, Problem 1c: Line Profiles at $\tau_c = %d$ star of $T_{eff} = 5780 K$. Tom Rice." % tau_c)
        ax.set_xlabel(r"Frequency $\nu$ (Hz)")
        ax.set_ylabel(r"$F_\nu / F_c$")

        ax.legend()

        ax.set_xlim(nu_array.min(), nu_array.max())
        ax.set_ylim(0,1)

        break

    plt.show()
        
    return (fig1, fig2, zeta_list, 
            equivalent_width_array, equivalent_width_array_lambda)

