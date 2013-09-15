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

frequency_array = np.logspace(12, 17) # In units of Hz. 
wavelength_array = c / frequency_array #Spans a wavelength range <10nm to >100um

for tau in tau_array:

    flux_array = np.zeros_like(frequency_array)
    for i, nu in zip( range(len(flux_array)), frequency_array):

        flux_array[i] = flux(tau, nu, T_eff)

    plt.plot(wavelength_array, flux_array*frequency_array, label=r"$\tau$ = %s" % tau)


plt.xlabel("Wavelength (meters)")
plt.ylabel("Flux * frequency")

plt.title("Problem Set 1, problem #3b")

plt.legend(loc = 'lower right')
plt.loglog()
plt.show()
