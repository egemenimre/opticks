# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
__version__ = "0.1.0"

from pathlib import Path

from pint import UnitRegistry

# Init units
u = UnitRegistry()
Q_ = u.Quantity

# define cycles in the MTF context
u.define("cycle = 1 * turn = cy")
