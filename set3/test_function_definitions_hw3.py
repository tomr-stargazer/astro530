
from function_definitions_hw3 import *

def test_u_function():

    nu = 1e15
    nu_0 = 2e15
    delta_nu = 1e14

    assert u_function(nu, nu_0, delta_nu) == -10

    
def test_hjerting():

    assert hjerting_gaussian_approximation(100, 5) == np.exp(-5**2)

    assert (hjerting_H1_approximation(1e-3, 5) == 
            (1e-3 * (0.56419 / 5**2 + 0.846 / 5**4)) )

    assert (hjerting_piecewise_approximation(1e-3, 1) == 
            hjerting_gaussian_approximation(100, 1))

    assert (hjerting_piecewise_approximation(1e-3, 5) ==
            hjerting_H1_approximation(1e-3, 5) + 
            hjerting_gaussian_approximation(1e-3, 5))
    
