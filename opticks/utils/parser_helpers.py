# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Package for Pydantic and YAML helpers.

"""

from typing import Annotated

from astropy.units import Quantity
from pydantic import AfterValidator, BeforeValidator, PlainSerializer

from opticks import u

# add bpp as a unit
u.add_enabled_units([u.def_unit("bpp", u.bit / u.pix)])


def is_quantity(quantity_text: str) -> bool:
    """
    Checks whether the `quantity_text` can be parsed as a valid `Quantity` object.

    Parameters
    ----------
    quantity_text : str
        Text that will be parsed as a `Quantity` object

    Returns
    -------
    is_quantity : bool
        `True` if text can be parsed, `False` otherwise
    """
    try:
        Quantity(quantity_text)
    except ValueError:
        return False
    else:
        return True


def _parse_quantity(v):
    """Pydantic BeforeValidator for Quantity objects."""
    if isinstance(v, Quantity):
        return v
    # Only Python 3.6+ supports underscores in numeric literals
    return Quantity(str(v).replace("_", ""))


def _serialize_quantity(v: Quantity) -> str:
    """Pydantic PlainSerializer for Quantity objects."""
    return str(v)


PydanticQty = Annotated[
    Quantity,
    BeforeValidator(_parse_quantity),
    PlainSerializer(_serialize_quantity, return_type=str),
]
"""Annotated Quantity type for use in Pydantic models.

Parses strings like '12905 mm' into astropy Quantity objects,
and serializes them back to strings.
"""


def _validate_positive_quantity(v: Quantity) -> Quantity:
    """Pydantic AfterValidator that requires a strictly positive Quantity."""
    if v.value <= 0:
        raise ValueError(f"Value must be positive, got {v}")
    return v


PositivePydanticQty = Annotated[
    Quantity,
    BeforeValidator(_parse_quantity),
    AfterValidator(_validate_positive_quantity),
    PlainSerializer(_serialize_quantity, return_type=str),
]
"""Annotated Quantity type for use in Pydantic models, restricted to positive values."""
