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

from function_definitions import flux
from function_definitions_hw3 import *

def problem_1a():
    """
    Calculate the emergent flux at 5000 Angstroms (0.5 microns) for
    a gray atmosphere, using the Eddington approximation, 
    and an effective temperature of 5780 K, with the source function
    as the Planck function.

    """

    wavelength = Quantity(5000, 'Angstrom')
    frequency = const.c/wavelength

    print frequency.to('Hz')
    
    result = flux(0, frequency, 5780)

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
