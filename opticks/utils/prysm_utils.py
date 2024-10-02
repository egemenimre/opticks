# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Package for prysm integration utilities and helpers.

"""

import copy

from numpy import ndarray
from pint import Quantity, Unit
from prysm._richdata import RichData
from prysm.coordinates import cart_to_polar, make_xy_grid

from opticks import u


class Grid:

    x: ndarray = None
    """Cartesian x coordinate."""

    y: ndarray = None
    """Cartesian y coordinate."""

    _r: ndarray = None
    """Radial coordinate."""

    _t: ndarray = None
    """Azimuthal coordinate"""

    dx: Quantity | float
    """Inter-sample spacing."""

    shape = None
    """Number of samples per dimension.

    If a scalar value, broadcast to both dimensions.
    Order is numpy axis convention, (row, col).
    Type is int or tuple of int."""

    def __init__(self, shape, dx: Quantity | float) -> None:
        """Create an x, y grid from -1, 1 with n number of samples.

        Parameters
        ----------
        shape : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions.  Order is numpy axis convention, (row, col)
        dx : Quantity | float
            inter-sample spacing

        """
        # if True, return meshgrid of x,y; else return 1D vectors (x, y)
        grid = True

        self.shape = shape
        self.dx = dx

        # Cartesian grid
        self.x, self.y = make_xy_grid(shape, dx=dx, grid=grid)

    @classmethod
    def from_size(cls, shape, size: Quantity | float) -> "Grid":
        """Generates a `Grid` object from the size of one side.

        If sample sizes are not equal in the two dimensions,
        the larger sample size is taken:
        dx = size / max(shape)

        Parameters
        ----------
        shape : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions. Order is numpy axis convention, (row, col)
        size : Quantity | float
            size of the grid on one side

        Returns
        -------
        Grid
            new Grid object
        """
        # convert shape to tuple if single int
        if not isinstance(shape, tuple):
            shape = (shape, shape)

        dx = size / max(shape)

        return Grid(shape, dx)

    def polar(self) -> tuple[ndarray, ndarray]:
        """Gets the polar grid with the given internal sample points.

        Returns
        -------
        rho, phi: tuple[ndarray, ndarray]
            tuple of radial coordinate and azimuthal coordinate, respectively
        """

        # lazy init polar coords
        if self._r is None:
            self._r, self._t = cart_to_polar(self.x, self.y)

        return self._r, self._t

    def cartesian(self) -> tuple[ndarray, ndarray]:
        """Gets the cartesian grid with the given internal sample points.

        Returns
        -------
        x, y: tuple[ndarray, ndarray]
            tuple of x and y cartesian coordinates, respectively
        """
        return self.x, self.y

    @property
    def r(self) -> ndarray:
        """Gets the radial coordinates."""

        # lazy init polar coords
        if self._r is None:
            self.polar()

        return self._r

    @property
    def t(self) -> ndarray:
        """Gets the azimuthal coordinates."""

        # lazy init polar coords
        if self._t is None:
            self.polar()

        return self._t

    def strip_units(self, units: Unit = Unit("mm")) -> "Grid":
        """Converts the Grid object to the default units.

        Converts the internal parameters to float ndarrays
        and returns a deepcopy of the `Grid` object.

        Azimuthal coordinates are in radians.

        Parameters
        ----------
        units : Unit, optional
            requested unit, by default "mm"

        Returns
        -------
        Grid
            Grid object with float ndarrays
        """

        dx = self.dx

        if self.has_units:
            # we already have units, convert them to the requested ones
            dx = dx.m_as(units)

        return Grid(self.shape, dx)

    @property
    def has_units(self) -> bool:
        """Checks whether the Grid internal data have units or
        are plain float ndarrays."""
        if isinstance(self.x, Quantity):
            return True
        else:
            return False


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
