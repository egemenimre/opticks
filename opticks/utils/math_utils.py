# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from enum import Enum

from astropy.units import FunctionUnitBase, Quantity, UnitBase
from numpy import ndarray
from scipy.interpolate import (
    Akima1DInterpolator,
    CubicHermiteSpline,
    CubicSpline,
    PchipInterpolator,
    PPoly,
)

from opticks.utils.unit_utils import (
    quantity_from_list,
    split_value_and_force_unit,
    split_value_and_unit,
)


class InterpolatorWithUnitTypes(Enum):
    """Enumerator for the interpolator types available."""

    AKIMA = Akima1DInterpolator
    """Akima1D Interpolator with default method (monotonic)."""
    MAKIMA = Akima1DInterpolator
    """Akima1D Interpolator with Modified Akima method (monotonic)."""
    PCHIP = PchipInterpolator
    """PCHIP Interpolator (monotonic)."""
    CUBICSPL = CubicSpline
    """Cubic Spline Interpolator."""
    CUBICHERMITESPL = CubicHermiteSpline
    """Cubic Hermite Spline Interpolator."""


class PPolyWithUnits(PPoly):
    """Subclass of `PPoly`, adding units functionality.

    The cubic and monotone interpolators in
    scipy all derive from `PPoly`, therefore this is
    equivalent to an interpolator with units.
    The interpolator is initialised elsewhere and
    its coefficients are used to initialise this object.

    Piecewise polynomial in terms of coefficients and breakpoints

    The polynomial between ``x[i]`` and ``x[i + 1]`` is written in the
    local power basis::

        S = sum(c[m, i] * (xp - x[i])**(k-m) for m in range(k+1))

    where ``k`` is the degree of the polynomial.

    Parameters
    ----------
    c : ndarray, shape (k, m, ...)
        Polynomial coefficients, order `k` and `m` intervals.
    x : ndarray, shape (m+1,)
        Polynomial breakpoints. Must be sorted in either increasing or
        decreasing order.
    extrapolate : bool or 'periodic', optional
        If bool, determines whether to extrapolate to out-of-bounds points
        based on first and last intervals, or to return NaNs. If 'periodic',
        periodic extrapolation is used. Default is True.
    axis : int, optional
        Interpolation axis. Default is zero.
    x_unit : UnitBase or FunctionUnitBase, optional
        unit of the x axis, by default None
    y_unit : UnitBase or FunctionUnitBase, optional
        unit of the y axis, by default None
    """

    def __init__(
        self,
        c,
        x,
        extrapolate=None,
        axis=0,
        x_unit: UnitBase | FunctionUnitBase = None,
        y_unit: UnitBase | FunctionUnitBase = None,
    ):

        super().__init__(c, x, extrapolate, axis)

        self.x_unit = x_unit
        self.y_unit = y_unit

    @classmethod
    def from_ppoly(
        cls,
        ppoly: type[PPoly],
        x_unit: UnitBase | FunctionUnitBase = None,
        y_unit: UnitBase | FunctionUnitBase = None,
    ) -> "PPolyWithUnits":
        """Upgrade an existing `PPoly` object (or from a subclass)
        with units.

        Reinitialises a new object by shallow copying the
        coefficients, breakpoints and other properties.

        Parameters
        ----------
        ppoly : type[PPoly]
            `PPoly` object (or from a subclass)
        x_unit : UnitBase or FunctionUnitBase, optional
            unit of the x axis, by default None
        y_unit : UnitBase or FunctionUnitBase, optional
            unit of the y axis, by default None

        Returns
        -------
        PPolyWithUnits
            PPoly object with units
        """
        return PPolyWithUnits(
            ppoly.c, ppoly.x, ppoly.extrapolate, ppoly.axis, x_unit, y_unit
        )

    def __call__(self, x: float | Quantity, nu=0, extrapolate=None, equivalencies=[]):

        # check x unit and generate the unitless version of x
        x, _ = split_value_and_force_unit(x, self.x_unit, equivalencies=equivalencies)

        # run the method with the unitless version of x
        # and output the result in y_unit
        return Quantity(super().__call__(x, nu, extrapolate), self.y_unit, copy=False)

    def derivative(self, nu=1) -> "PPolyWithUnits":

        y_unit = self.y_unit / self.x_unit**nu

        # compute and add the units
        return PPolyWithUnits.from_ppoly(super().derivative(nu), self.x_unit, y_unit)

    def antiderivative(self, nu=1) -> "PPolyWithUnits":

        y_unit = self.y_unit * self.x_unit**nu

        # compute and add the units
        return PPolyWithUnits.from_ppoly(
            super().antiderivative(nu), self.x_unit, y_unit
        )

    def integrate(self, a, b, extrapolate=None, equivalencies=[]):

        a, _ = split_value_and_force_unit(a, self.x_unit, equivalencies=equivalencies)
        b, _ = split_value_and_force_unit(b, self.x_unit, equivalencies=equivalencies)

        y_unit = self.y_unit * self.x_unit

        return Quantity(super().integrate(a, b, extrapolate), y_unit, copy=False)

    def roots(self, discontinuity=True, extrapolate=None) -> ndarray:

        return self.solve(0 * self.y_unit, discontinuity, extrapolate)

    def solve(
        self,
        y: float | Quantity,
        discontinuity=True,
        extrapolate=None,
        equivalencies=[],
    ) -> ndarray:

        y, _ = split_value_and_force_unit(y, self.y_unit, equivalencies=equivalencies)

        return Quantity(
            super().solve(y, discontinuity, extrapolate), self.x_unit, copy=False
        )

    def __str__(self) -> str:

        out = f"Interpolator in range: [{self.x[0] * self.x_unit:P~}, {self.x[-1] * self.x_unit:P~}]"

        out += f" (extrapolation: {self.extrapolate})."

        return out


class InterpolatorWithUnits(PPolyWithUnits):
    """Interpolator with units.

    For most usecases, the `from_ipol_method` factory
    is more convenient than this constructor.

    Parameters
    ----------
    ppoly : type[PPoly]
        Interpolator subclassing `PPoly`
    x_unit : UnitBase or FunctionUnitBase, optional
        unit of the x axis, by default None
    y_unit : UnitBase or FunctionUnitBase, optional
        unit of the y axis, by default None
    """

    def __init__(
        self,
        ipol: type[PPoly],
        x_unit: UnitBase | FunctionUnitBase = None,
        y_unit: UnitBase | FunctionUnitBase = None,
    ):

        super().__init__(ipol.c, ipol.x, ipol.extrapolate, ipol.axis, x_unit, y_unit)

    @classmethod
    def from_ipol_method(
        cls,
        ipol_type: InterpolatorWithUnitTypes,
        x: ndarray[float | Quantity],
        y: ndarray[float | Quantity],
        *args,
        axis=0,
        extrapolate=None,
        **kwargs,
    ) -> "InterpolatorWithUnits":
        """Generates and interpolator with units.

        The interpolator is chosen with the enum `InterpolatorWithUnitTypes`,
        among scipy cubic and monotonic interpolators. Interpolator
        specific parameters can be passed via `*args, **kwargs`,
        with the exception of the `method` parameter in `Akima1DInterpolator`.
        Use the enum `MAKIMA` for the Modified Akima instead.

        Parameters
        ----------
        ipol_type : InterpolatorWithUnitTypes
            Interpolator type
        x : array_like, shape (n,)
            1-D array containing values of the independent variable.
            Values must be real, finite and in strictly increasing order.
        y : array_like
            Array containing values of the dependent variable. It can have
            arbitrary number of dimensions, but the length along ``axis``
            (see below) must match the length of ``x``. Values must be finite.
        axis : int, optional
            Axis along which `y` is assumed to be varying. Meaning that for
            ``x[i]`` the corresponding values are ``np.take(y, i, axis=axis)``.
            Default is 0.
        extrapolate : {bool, 'periodic', None}, optional
            If bool, determines whether to extrapolate to out-of-bounds points
            based on first and last intervals, or to return NaNs. If 'periodic',
            periodic extrapolation is used. If None (default), it is set to True.

        Returns
        -------
        InterpolatorWithUnits
            Interpolator with units
        """

        klass = ipol_type.value

        # if Quantity input is used, it should be in the form
        # "[0.0 1.0] meter", not a list of individual Quantity objects.
        # This checks and corrects it.
        x = quantity_from_list(x)
        y = quantity_from_list(y)

        # split the value and the units
        x, x_unit = split_value_and_unit(x)
        y, y_unit = split_value_and_unit(y)

        # Init interpolator with the unitless data
        if ipol_type == InterpolatorWithUnitTypes.AKIMA:
            ipol = klass(x, y, axis=axis, extrapolate=extrapolate, method="akima")
        elif ipol_type == InterpolatorWithUnitTypes.MAKIMA:
            ipol = klass(x, y, axis=axis, extrapolate=extrapolate, method="makima")
        else:
            ipol = klass(x, y, *args, **kwargs)

        # assign units to x and y
        return InterpolatorWithUnits(ipol, x_unit, y_unit)
