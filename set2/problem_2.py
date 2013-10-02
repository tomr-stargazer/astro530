"""
Code relating to problem 1 in Astro 530 HW #2.

Problem 2 is about the Saha equation and the number of electrons.

"""

from __future__ import division

import numpy as np

from problem_1 import compute_grid_of_number_densities

tables = compute_grid_of_number_densities()

for table, n in zip(tables, range(1,4)):

    # sums up the number of electrons freed by each ion species.
    n_e = sum([np.nansum(table['N_%d' % i] * i) for i in range(6)])

    print n_e, "Model %d" % n
