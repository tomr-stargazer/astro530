from __future__ import division

from problem_1 import *

k_B = const.k_B

def test_saha_equation():

    # This test uses Caroll & Ostlie p. 218, example 8.1.5
    lower_state_degeneracy = 2
    upper_state_degeneracy = 1

    chi_01 = 13.6

    P_e = Quantity(20, 'N / m2')

    # What's n_e for a given P_e and T?
    n_e = lambda T: P_e / (k_B*Quantity(T, 'K'))

    # Test a bunch of ionization fractions for various temperatures!
    for T in np.arange(5,15,0.5)*1000.:
    
        saha = saha_equation( lower_state_degeneracy,
                              upper_state_degeneracy,
                              chi_01, T, n_e(T))

        print T, saha, (saha / (1.0+saha.value))

        # Caroll & Ostlie say that when P_e = 20 N m^-2,
        # 50% ionization occurs at T~9600. Let's make this the test.
        if T < 9600:
            assert saha.value/(1+saha.value) < .5
        elif T > 9600:
            assert saha.value/(1+saha.value) > .5
