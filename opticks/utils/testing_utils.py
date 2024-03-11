# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Package for testing utilities.

"""
from typing import Optional

import numpy as np
from pint import Quantity
from pint.testing import assert_allclose as pint_assert_allclose
from pint.testing import assert_equal as pint_assert_equal


def assert_equal(first: Quantity, second: Quantity, msg: Optional[str] = None) -> None:
    """
    Compares `first` and `second` items and checks whether they are equal.

    `first` and `second` are array-like objects.

    Essentially a wrapper for the pint `assert_equal`.

    Parameters
    ----------
    first : Quantity
        item with numpy array like Quantity object
    second : Quantity
        item with numpy array like Quantity object
    msg : Optional[str]
        human-readable message
    """
    # just a wrapper for the pint assert_equal
    pint_assert_equal(first, second, msg=msg)


def assert_allclose(
    first: Quantity,
    second: Quantity,
    rtol: float = 1e-07,
    atol: Quantity = None,
    msg: Optional[str] = None,
) -> None:
    """
    Compares `first` and `second` items and checks whether they are equal within a certain tolerance.

    `first` and `second` are array-like objects. Tolerance can be relative (dimensionless) or absolute (`Quantity`).

    Essentially a wrapper for the pint `assert_allclose`.

    Parameters
    ----------
    first : Quantity
        item with numpy array like Quantity object
    second : Quantity
        item with numpy array like Quantity object
    rtol : float
        relative tolerance
    atol : `Quantity`
        absolute tolerance with the same dimensionality as `first`and `second`
    msg : Optional[str]
        human-readable message
    """
    if atol is None:
        # Use rtol or atol for errors given in percent or absolute.
        # atol does not accept units, but assumed by pint as the unit of the first item.
        pint_assert_allclose(first, second, rtol=rtol, msg=msg)

    else:

        # atol should be a scalar
        if np.isscalar(atol.m):

            # If atol is requested, convert atol to the units of the first item
            # and get the magnitude. If atol and first are incompatible, will raise
            # a DimensionalityError.
            atol = atol.to(first.u).m

            # Use rtol or atol for errors given in percent or absolute.
            # atol does not accept units, but assumed by pint as the unit of the first item.
            pint_assert_allclose(first, second, rtol=rtol, atol=atol, msg=msg)

        else:
            raise ValueError("atol is not a scalar")
