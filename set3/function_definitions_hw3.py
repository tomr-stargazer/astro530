

from __future__ import division

import numpy as np


def u_function(frequency, central_frequency, doppler_width):
    """
    u is a parameterized form of frequency that is used in the Hjerting function.

    It's differenced from the line's central frequency and then divided
    by the local Doppler width of the line.

    Make sure you provide all the inputs in the same units (e.g. Hz).

    Parameters
    ----------
    frequency : float
        The frequency $\nu$ you want to calculate u for.
    central_frequency : float
        The central frequency $\nu_0$ of the spectral line we're looking at.
    doppler_width : float
        The doppler width $\Delta \nu_D$ at the local level, in frequency units.

    Returns
    -------
    u : float
        The `u` parameter that is used in the Hjerting function.
        It's unitless.

    """

    nu = frequency
    nu_0 = central_frequency
    delta_nu = doppler_width

    return (nu - nu_0) / delta_nu
