"""
Does computational things for problem 1 in HW #3 in astro530.

Borrows from the code I wrote in HW #1, so I'm really glad I spent 
the time to make it high quality and tested, etc back then.

"""

from __future__ import division

import numpy as np
from scipy.integrate import quad, trapz
from scipy.constants import c
import scipy.special
import astropy.constants as const
from astropy.units.quantity import Quantity

from function_definitions import flux

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
