""" Tests problem1 set6. """

from astro530.set6.problem_1_set6 import *

def test_dust_opacity_at_nu():

    assert dust_opacity_at_nu(1e12) == 0.1
    

    assert dust_opacity_at_nu(1e15) == 0.1 * 1e3
