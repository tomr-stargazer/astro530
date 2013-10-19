
from function_definitions_hw3 import *

def test_u_function():

    nu = 1e15
    nu_0 = 2e15
    delta_nu = 1e14

    assert u_function(nu, nu_0, delta_nu) == -10

    
