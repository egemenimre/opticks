# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
__version__ = "0.1.0"

from typing import TypeAlias

import pint
import pint.registry_helpers
import portion as P
from pint.facets.plain.quantity import check_implemented, ireduce_dimensions


class PortionCompQuantity(pint.UnitRegistry.Quantity):
    """Quantity modified for Portion Compatibility."""

    @check_implemented
    def compare(self, other, op):
        # sort out the portion P.inf incompatibility
        # check for portion Inf and replace with Quantity object
        if isinstance(other, P.const._PInf):
            other = self._REGISTRY.Quantity("inf") * self.u
        elif isinstance(other, P.const._NInf):
            other = self._REGISTRY.Quantity("-inf") * self.u

        return super().compare(other, op)

    @check_implemented
    @ireduce_dimensions
    def _mul_div(self, other, magnitude_op, units_op=None):
        # check for portion Inf and replace with Quantity object
        other = self.__convert_portion_inf__(other)

        return super()._mul_div(other, magnitude_op, units_op)

    def __convert_portion_inf__(self, other):
        if isinstance(other, P.const._PInf):
            return self._REGISTRY.Quantity("inf")
        elif isinstance(other, P.const._NInf):
            return self._REGISTRY.Quantity("-inf")
        else:
            return other


class PortionCompRegistry(
    pint.registry.GenericUnitRegistry[PortionCompQuantity, pint.UnitRegistry.Unit]
):
    """Registry modified for Portion Compatibility."""

    Quantity: TypeAlias = PortionCompQuantity
    Unit: TypeAlias = pint.UnitRegistry.Unit

    # Explicit super() here, otherwise `setup_matplotlib` not found
    def setup_matplotlib(self, enable: bool = True) -> None:
        pint.UnitRegistry.setup_matplotlib(self, enable)

    # Explicit super() here, otherwise `check` and `wraps` not found
    check = pint.registry_helpers.check
    # wraps = registry_helpers.wraps
    wraps = pint.registry_helpers.wraps


# Init units with Portion compatibility
# u = pint.UnitRegistry()
u = PortionCompRegistry()
Q_ = u.Quantity

# define cycles in the MTF context
u.define("cycle = 1 * turn = cy")
