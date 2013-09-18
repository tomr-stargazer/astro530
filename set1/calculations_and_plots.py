"""
Module that uses function_definitions to make plots and do calculations.

"""

from __future__ import division

import matplotlib.pyplot as plt
from scipy.constants import c

from function_definitions import *

# Some constants

T_eff = 5780 # kelvin; similar to the Sun
tau_array = [3, 1, 0.3, 0.1]

def problem_3b(spacing=20):
    """
    Calculates and plots the figure for 3b.

    """
    
    frequency_array = np.logspace(12, 17, num=spacing) # In units of Hz. 
    wavelength_array = c / frequency_array #Spans a wavelength range <10nm to >100um

    # Make a figure
    fig = plt.figure()

    for tau in tau_array:

        flux_array = np.zeros_like(frequency_array)
        for i, nu in zip( range(len(flux_array)), frequency_array):

            flux_array[i] = flux(tau, nu, T_eff)

        plt.plot( wavelength_array*1e6, 
                  flux_array*frequency_array, 
                  label=r"$\tau$ = %s" % tau,
                  lw=2)

    plt.xlabel(r"Wavelength $\lambda$ (microns $\mu m$)")
    plt.ylabel(r"$\nu \cdot \mathcal{F}_\nu$ (Flux times frequency)")

    plt.title("Problem Set 1, problem #3b")

    plt.legend(loc = 'lower right')
    plt.loglog()

    plt.xlim(5e-2, 1e1)
    plt.ylim(1, 1e8)

    plt.text(0.3, 1e4, r"$T_{eff} = 5780$ K", fontsize=18)

    plt.show()


def problem_3c():
    """
    Calculates and plots stuff for problem 3c.

    """
    tau_array = np.logspace(-5, 1, 6)

    for tau in tau_array:

        flux_at_tau = lambda nu: flux(tau, nu, T_eff)

        integrated_flux_at_tau = quad(flux_at_tau, 1e10, 1e18)[0]

        print "Integrated flux at tau=%s: \n %s" % (tau, integrated_flux_at_tau)

    
