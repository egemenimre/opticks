# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Package for prysm integration utilities and helpers.

"""

import copy

from astropy.units import Quantity, Unit
from numpy import ndarray
from prysm._richdata import RichData
from prysm.coordinates import cart_to_polar, make_xy_grid
from prysm.polynomials import ansi_j_to_nm, sum_of_2d_modes, zernike_nm_sequence

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

    dx: Quantity
    """Inter-sample spacing."""

    shape = None
    """Number of samples per dimension.

    If a scalar value, broadcast to both dimensions.
    Order is numpy axis convention, (row, col).
    Type is int or tuple of int."""

    def __init__(self, shape, dx: Quantity) -> None:
        """Create an x, y grid from -1, 1 with n number of samples.

        Parameters
        ----------
        shape : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions.  Order is numpy axis convention, (row, col)
        dx : Quantity
            inter-sample spacing

        """
        # if True, return meshgrid of x,y; else return 1D vectors (x, y)
        grid = True

        self.shape = shape
        self.dx = dx

        # Cartesian grid
        self.x, self.y = make_xy_grid(shape, dx=dx, grid=grid)

    @classmethod
    def from_size(cls, shape, size: Quantity) -> "Grid":
        """Generates a `Grid` object from the size of one side.

        If sample sizes are not equal in the two dimensions,
        the larger sample size is taken:
        dx = size / max(shape)

        Parameters
        ----------
        shape : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions. Order is numpy axis convention, (row, col)
        size : Quantity
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
            New Grid object with float ndarrays
        """

        if self.has_units:
            # we already have units, convert them to the requested ones
            dx = self.dx.to_value(units)
        else:
            # deep copy internal data
            dx = copy.deepcopy(self.dx)

        return Grid(self.shape, dx)

    @property
    def has_units(self) -> bool:
        """Checks whether the internal data have units or
        are plain float ndarrays."""
        if isinstance(self.x, Quantity):
            return True
        else:
            return False


class OptPathDiff:

    def __init__(self, opd: ndarray):
        """Init with Optical Path Difference.

        The input data should be in a format that can be used
        in `prysm`. Therefore, it is recommended to initialise
        this object with `from_zernike` or similar methods.

        Parameters
        ----------
        opd : ndarray
            Optical Path Difference data (in nm)
        """
        self.data = opd

    @classmethod
    def from_zernike(
        cls, wfe_rms: list[Quantity], aperture_diameter: Quantity, grid: Grid
    ) -> "OptPathDiff":
        """Computes the Optical Path Difference via Zernike Polynomials.

        Generates the Zernike Polynomials to the order corresponding
        to the number of coefficients (e.g., 9 coefficients = mode 8)
        sums them properly, adding the WFE RMS Zernike coefficients.

        The result is the monochromatic OPD for a single location
        on the PSF plane.

        Parameters
        ----------
        wfe_rms : list[Quantity]
            ANSI list of aberration coefficients (WFE RMS) (in nm)
        diameter : Quantity
            aperture diameter in mm
        grid : Grid
            aperture grid (in mm and rad)

        Returns
        -------
        OptPathDiff
            Optical Path Difference (OPD)
        """

        # mode n = (n+1) elements
        elems = len(wfe_rms)

        # Generate the (n,m) tuples in ANSI order
        nms = [ansi_j_to_nm(i) for i in range(0, elems)]

        # radial coords normalised by aperture radius
        # normalisation required by the polynomials
        ap_radius = aperture_diameter / 2.0

        r, t = grid.polar()

        # reduce to dimensionless
        rho = (r / ap_radius).decompose()

        # compute the polynomials (dimensionless)
        # t should be in radians
        mode = list(zernike_nm_sequence(nms, rho.value, t.value))

        # monochromatic OPD with multiple aberrations
        # wfe_rms and opd units are should be in nm
        opd = sum_of_2d_modes(mode, wfe_rms.to_value(u.nm))

        return OptPathDiff(opd * u.nm)

    def strip_units(self, units: Unit = Unit("nm")) -> "OptPathDiff":
        """Converts the OptPhaseDiff object to the default units.

        Converts the internal parameters to float ndarrays
        and returns a deepcopy of the `OptPhaseDiff` object.

        Parameters
        ----------
        units : Unit, optional
            requested unit, by default "nm"

        Returns
        -------
        OptPathDiff
            New OptPathDiff object with float ndarrays
        """

        if self.has_units:
            # we already have units, convert them to the requested ones
            opd = self.data.to_value(units)
        else:
            # deep copy internal data without units
            opd = copy.deepcopy(self.data)

        return OptPathDiff(opd)

    @property
    def has_units(self) -> bool:
        """Checks whether the internal data have units or
        are plain float ndarrays."""
        if isinstance(self.data, Quantity):
            return True
        else:
            return False


def richdata_has_units(rich_data: RichData) -> bool:
    """Checks whether `RichData` object has units."""
    if isinstance(rich_data.dx, Quantity):
        return True
    else:
        return False


def richdata_with_units(rich_data: RichData, dx_units: Unit = Unit("mm")) -> RichData:
    """Generates a deepcopy of the input `RichData` object with units.

    Adds units to `dx` and `wavelength` (if available).
    The `data` structure is already without units by definition.

    If the input `rich_data` has units, raises a `ValueError`
    exception.

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

    # first check for units
    if richdata_has_units(rich_data):
        raise ValueError("Input RichData object already has units.")

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
