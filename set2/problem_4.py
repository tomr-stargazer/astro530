"""
Problem 4.

Plotting a thing.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

def J_nu(tau, a=1, b=0.5, lamda=0.01):

    first_term = ((b / np.sqrt(3) - a)/(1 + np.sqrt(lamda)) * 
                  np.exp( - np.sqrt(lamda * 3) * tau))
    second_term = a + b*tau

    return first_term + second_term

def plot_J_nu():

    pass

    
