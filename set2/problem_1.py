"""
Code relating to problem 1 in Astro 530 HW #2.

Problem 1 is about the Saha equation and ions of various atoms
in three different stellar atmospheres.

The Saha equation:
https://en.wikipedia.org/wiki/Saha_ionization_equation

..math::
\frac{n_{i+1}n_e}{n_i} = \frac{2}{\Lambda^{3}}\frac{g_{i+1}}{g_i}\exp\left[-\frac{(\epsilon_{i+1}-\epsilon_i)}{k_BT}\right]

"""

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
