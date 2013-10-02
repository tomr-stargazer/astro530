"""
Code relating to problem 1 in Astro 530 HW #2.

Problem 1 is about the Saha equation and ions of various atoms
in three different stellar atmospheres.

The Saha equation:
https://en.wikipedia.org/wiki/Saha_ionization_equation

but I have followed Carroll and Ostlie's choice to place 
the electron number density $n_e$ on the denominator of the right 
side of the equation, rather than the numerator of the left side.

Thus:

..math::
\frac{n_{i+1}}{n_i} = \frac{2}{n_e\Lambda^{3}}\frac{g_{i+1}}{g_i}\exp\left[-\frac{(\epsilon_{i+1}-\epsilon_i)}{k_BT}\right]

"""

from __future__ import division

import numpy as np

import astropy.constants as const
from astropy.units.quantity import Quantity
from astropy.table import Table, Column


Model_1 = {'T_eff': 4500,
           'log g': 4.0,
           'T': 4469.7,
           'n_e': 4.75e12,
           'n_H': 6.6e16}

Model_2 = {'T_eff': 5750,
           'log g': 4.0,
           'T': 5609.3,
           'n_e': 9.55e12,
           'n_H': 6.0e16}

Model_3 = {'T_eff': 15000,
           'log g': 4.0,
           'T': 14630,
           'n_e': 3.17e14,
           'n_H': 3.67e14}

model_list = [Model_1, Model_2, Model_3]


def saha_equation(degeneracy_of_lower_state,
                  degeneracy_of_upper_state,
                  ionization_potential_between_states,
                  temperature,
                  electron_number_density,
                  units_of_energy='eV'):
    """
    Calculates the Saha equation for a given set of parameters.

    Returns
    -------
    saha_ratio : float
        The ratio n_{i+1} / n_i, the output of the Saha equation.

    """

    g_i = degeneracy_of_lower_state
    g_ip1 = degeneracy_of_upper_state

    chi_i = Quantity(ionization_potential_between_states, units_of_energy)

    n_e = Quantity(electron_number_density, '1 / (cm3)')
    T = Quantity(temperature, 'K')

    h = const.h.to(units_of_energy+' s')
    k_B = const.k_B.to(units_of_energy+' / (K)')

    # "Lambda" is just a name for this bundle of parameters and constants.
    lambda_squared = (h**2 / (2 * np.pi * const.m_e * k_B * T)).to('cm2')

    saha_ratio = (( 2 / lambda_squared**(3/2) *
                    g_ip1 / g_i *
                    np.exp(-(chi_i) /(k_B * T))
                    ) / n_e).value

    return saha_ratio

def number_density_of_each_state_from_saha(saha_list, n_H, abundance):
    """
    Gives number density of each ionization state given ALL of the
    outputs of the Saha equation for one atom.

    Thanks enormously to Danny Zhu for recommending I not turn this
    into a total recursive clusterfuck. Multiplication's easier.

    "(if N=5 and the ratios are r_1,r_2,r_3,r_4) 
    {1, r_1, r_1 r_2, r_1 r_2 r_3, r_1 r_2 r_3 r_4} / (1 + r_1 + r_1 r_2 + r_1 r_2 r_3 + r_1 r_2 r_3 r_4) ?" --dzhu

    """
    
    numerator_array = np.cumprod([1]+saha_list)
    denominator = reduce(np.add, numerator_array)

    return numerator_array/denominator * abundance * n_H

def compute_grid_of_number_densities():
    """
    Computes a grid of number densities for each model. 

    Returns three tables to you.

    """

    # let's load in that table
    atomic_table = np.loadtxt('Atomic_data_table.txt', skiprows=2)

    atomic_weights = atomic_table[:,0]
    abundances = atomic_table[:,1]

    degeneracy_of_states = [atomic_table[:,x] for x in range(2,8)]
    ionization_energy_per_level = [atomic_table[:,x] for x in range(8,13)]

    output_table_list = []

    for model in model_list:

        model_table = Table()

        grid_of_densities = np.zeros((len(degeneracy_of_states),
                                      len(atomic_weights))) * np.nan

        print "Model parameters: "+str(model.items())

        for i in range(atomic_weights.size):

            # For this atom, save all of the saha ratios.
            saha_list = []

            for j in range(len(ionization_energy_per_level)):

                if ionization_energy_per_level[j][i] == 999:
                    continue
                
                saha_ratio = saha_equation(
                    degeneracy_of_states[j][i],
                    degeneracy_of_states[j+1][i],
                    ionization_energy_per_level[j][i],
                    model['T'], model['n_e'])

                saha_list.append(saha_ratio)

            number_densities_for_atom =number_density_of_each_state_from_saha(
                saha_list, model['n_H'], abundances[i])

            for j in range(len(number_densities_for_atom)):
                grid_of_densities[j][i] = number_densities_for_atom[j]


        model_table.add_column(Column(name="AW", data=atomic_weights))
        for i in range(grid_of_densities.shape[0]):
            model_table.add_column(Column(name="N_%d" % i,
                                          data=grid_of_densities[i]))

        output_table_list.append(model_table)
        
    return output_table_list

def write_outputs_to_disk():

    table1, table2, table3 = compute_grid_of_number_densities()

    table1.write("model1_grid.txt", format='ascii')
    table2.write("model2_grid.txt", format='ascii')
    table3.write("model3_grid.txt", format='ascii')


if __name__ == "__main__":
    write_outputs_to_disk()
