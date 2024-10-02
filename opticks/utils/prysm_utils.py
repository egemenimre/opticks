# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Package for prysm integration utilities and helpers.

"""

import copy

from pint import Unit
from prysm._richdata import RichData

from opticks import u


def richdata_with_units(rich_data: RichData, dx_units: Unit = Unit("mm")) -> RichData:
    """Generates a deepcopy of the input `RichData` object with units.

    Adds units to `dx` and `wavelength` (if available).
    The `data` structure is already without units by definition.

    Parameters
    ----------
    rich_data : RichData
        input object without units
    dx_units: Unit
        units to be used for `dx` spacing parameter

    Returns
    -------
    RichData
        output object with units
    """

    # data is a simple ndarray without units
    data: RichData = copy.deepcopy(rich_data.data)

    # inter-sample spacing, mm for PSF and cy/mm for MTF
    dx = rich_data.dx * dx_units

    # wavelength of light, um
    if rich_data.wavelength:
        wavelength = rich_data.wavelength * u.um
    else:
        wavelength = None

    return RichData(data, dx, wavelength)
