# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
__version__ = "0.1.0"

import numpy as np
import portion as Portion
from astropy import units as u
from astropy.units import Quantity


class UnitInterval(Portion.Interval):
    """Subclass the Interval class to handle units.

    This fixes the problem between the infinity definitions
    of Portion and the unit handler (astropy, pint etc.)."""

    @classmethod
    def from_atomic(cls, left, lower, upper, right):

        # change the inf to np.inf, assign units if possible
        if isinstance(lower, Portion.const._PInf):
            lower = np.inf
            if isinstance(upper, u.Quantity):
                lower = lower * upper.unit
        elif isinstance(lower, Portion.const._NInf):
            lower = -1 * np.inf
            if isinstance(upper, u.Quantity):
                lower = lower * upper.unit

        if isinstance(upper, Portion.const._PInf):
            upper = np.inf
            if isinstance(lower, u.Quantity):
                upper = upper * lower.unit
        elif isinstance(upper, Portion.const._NInf):
            upper = -1 * np.inf
            if isinstance(lower, u.Quantity):
                upper = upper * lower.unit

        return super().from_atomic(left, lower, upper, right)


# Create the Portion Interval API with the new shorthand
P = Portion.create_api(UnitInterval)

# shorthand for Quantity
Q_ = Quantity
