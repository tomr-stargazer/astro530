"""
Some basic unit tests for Problem 1.

"""

from problem_1_set5 import *

def test_sound_speed():

    by_hand_calculation_of_cs_for_T100_mu1 = 908.02
    
    assert (round(sound_speed(100,1), -1) ==
            round(by_hand_calculation_of_cs_for_T100_mu1, -1) )

    cs_test2 = np.sqrt(k_B * 50 / (2*m_p))

    assert sound_speed(50,2).value == cs_test2

    
