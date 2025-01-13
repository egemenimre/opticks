# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


import numpy as np
from astropy.units import FunctionUnitBase, Quantity, UnitBase
from numpy.typing import ArrayLike


def quantity_from_list(
    data: list[Quantity], unit: UnitBase | FunctionUnitBase = None
) -> Quantity:
    """Converts a list of values to a Quantity object.

    The list should contain only Quantity objects that can be
    mutually converted (for example all of them should be length).
    If there are different units in the list, the first unit is used.
    If `unit` is given, then the data is converted to that unit.

    If the data is not a list of Quantity objects, then returns
    the data with no changes.

    Parameters
    ----------
    data : list[Quantity]
        The data to be converted
    unit : UnitBase | FunctionUnitBase, optional
        Target unit

    Returns
    -------
    Quantity
        The data as a Quantity object
    """

    if isinstance(data, list) and any(isinstance(x, Quantity) for x in data):

        # find the first unit
        for x in data:
            if isinstance(x, Quantity):
                unit = x.unit
                # found one unit, break
                break

        # all items should be Quantity, otherwise this will raise an Exception
        return Quantity(data, unit)
    else:
        return data


def split_value_and_unit(
    data: Quantity | float | ArrayLike,
) -> tuple[float | np.ndarray, UnitBase | FunctionUnitBase]:
    """Splits the value and the unit (if available).

    The input object structure is preserved. If the input is a list or
    ndarray then the output will also be a list or ndarray.

    If the `data` is not a Quantity object, then the unit will be set to
    `None`. Initialising an object with `Quantity(5, None)` will just
    result in a dimensionless Quantity object.

    Parameters
    ----------
    data : Quantity | float | np.ArrayLike
        The data to be split

    Returns
    -------
    tuple[float | np.ndarray, UnitBase | FunctionUnitBase]
        The value (without units) and the unit (if available)
    """
    if not isinstance(data, Quantity):
        unit = None
        data_val = data
    else:
        data_val = data.value
        unit = data.unit

    return (data_val, unit)


def split_value_and_force_unit(
    data: Quantity | float | ArrayLike,
    tgt_unit: UnitBase | FunctionUnitBase,
    equivalencies=[],
) -> tuple[float | np.ndarray, UnitBase | FunctionUnitBase]:
    """Splits the value and the unit, converting the data to the target unit.

    If the data has no units, then the target unit is assigned to the
    data.

    Parameters
    ----------
    data : Quantity | float | np.ArrayLike
        The data to be split
    tgt_unit : _type_, optional
        target unit
    equivalencies : list, optional
        Equivalencies to be used during conversion

    Returns
    -------
    tuple[float | np.ndarray, UnitBase | FunctionUnitBase]
        The value (without units) and the unit
    """

    if not isinstance(data, Quantity):
        # no units given
        unit = tgt_unit
        data_val = data
    else:
        # force into target units
        data_val = data.to_value(tgt_unit, equivalencies=equivalencies)
        unit = tgt_unit

    return (data_val, unit)
