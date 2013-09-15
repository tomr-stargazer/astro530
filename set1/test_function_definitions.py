from function_definitions import *

def test_planck_function():

    # See: http://wolfr.am/1ghoEBo for Wolfram Alpha example result
    assert round(planck_function(6000, 200e12), 8+3) == 2.985e-8
    

def test_eddington_temperature():

    # Practically by definition, T(tau = 2/3) == T_eff
    assert eddington_temperature(2/3, 6000) == 6000
    assert eddington_temperature(2/3, 10) == 10

    assert eddington_temperature(1, 6000) > 6000
    assert eddington_temperature(1/2, 6000) < 6000

def test_E2():
    # See Wolfram Alpha
    assert round(E2(1), 4) == 0.1485
