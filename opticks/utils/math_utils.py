# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from enum import Enum

from numpy import ndarray
from pint import Quantity
from pint.facets.numpy.numpy_func import unwrap_and_wrap_consistent_units
from scipy.interpolate import (
    Akima1DInterpolator,
    CubicHermiteSpline,
    CubicSpline,
    PchipInterpolator,
    PPoly,
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
    x_unit : Quantity, optional
        unit of the x axis, by default None
    y_unit : Quantity, optional
        unit of the y axis, by default None
    """

    def __init__(
        self,
        c,
        x,
        extrapolate=None,
        axis=0,
        x_unit: Quantity = None,
        y_unit: Quantity = None,
    ):

        super().__init__(c, x, extrapolate, axis)

        # Extract output units
        _, out_wrap_x = unwrap_and_wrap_consistent_units(x_unit)
        _, out_wrap_y = unwrap_and_wrap_consistent_units(y_unit)

        self.x_unit = out_wrap_x(1)
        self.y_unit = out_wrap_y(1)

    @classmethod
    def from_ppoly(
        cls,
        ppoly: type[PPoly],
        x_unit: Quantity = None,
        y_unit: Quantity = None,
    ) -> "PPolyWithUnits":
        """Upgrade an existing `PPoly` object (or from a subclass)
        with units.

        Reinitialises a new object by shallow copying the
        coefficients, breakpoints and other properties.

        Parameters
        ----------
        ppoly : type[PPoly]
            `PPoly` object (or from a subclass)
        x_unit : Quantity, optional
            unit of the x axis, by default None
        y_unit : Quantity, optional
            unit of the y axis, by default None

        Returns
        -------
        PPolyWithUnits
            PPoly object with units
        """
        return PPolyWithUnits(
            ppoly.c, ppoly.x, ppoly.extrapolate, ppoly.axis, x_unit, y_unit
        )

    def __call__(self, x: float | Quantity, nu=0, extrapolate=None):

        # check x unit and generate the unitless version of x
        (_, x), _ = unwrap_and_wrap_consistent_units(self.x_unit, x)

        # generate output unit wrapper
        _, output_wrap = unwrap_and_wrap_consistent_units(self.y_unit)

        # run the method with the unitless version of x
        # and output the result in y_unit
        return output_wrap(super().__call__(x, nu, extrapolate))

    def derivative(self, nu=1) -> "PPolyWithUnits":

        # compute and add the units
        return PPolyWithUnits.from_ppoly(
            super().derivative(nu), self.x_unit, self.y_unit / self.x_unit
        )

    def antiderivative(self, nu=1) -> "PPolyWithUnits":

        # compute and add the units
        return PPolyWithUnits.from_ppoly(
            super().antiderivative(nu), self.x_unit, self.x_unit * self.y_unit
        )

    def integrate(self, a, b, extrapolate=None):

        (_, a, b), _ = unwrap_and_wrap_consistent_units(self.x_unit, a, b)
        _, output_wrap = unwrap_and_wrap_consistent_units(self.x_unit * self.y_unit)
        return output_wrap(super().integrate(a, b, extrapolate))

    def roots(self, discontinuity=True, extrapolate=None) -> ndarray:

        return self.solve(0 * self.y_unit, discontinuity, extrapolate)

    def solve(
        self, y: float | Quantity, discontinuity=True, extrapolate=None
    ) -> ndarray:

        (_, y), _ = unwrap_and_wrap_consistent_units(self.y_unit, y)
        _, output_wrap = unwrap_and_wrap_consistent_units(self.x_unit)
        return output_wrap(super().solve(y, discontinuity, extrapolate))


# TODO delete this after tests

# def __init__(
#     self,
#     c,
#     x,
#     extrapolate=None,
#     axis=0,
#     x_unit: Quantity = None,
#     y_unit: Quantity = None,
# ):

#     super().__init__(c, x, extrapolate, axis)

#     # force units even if non-dim
#     # x_unit and y_unit converted to magnitude of 1
#     if x_unit is None or not isinstance(x_unit, Quantity):
#         self.x_unit = Q_("dimensionless")
#     else:
#         self.x_unit = 1 * x_unit.u

#     if y_unit is None or not isinstance(y_unit, Quantity):
#         self.y_unit = Q_("dimensionless")
#     else:
#         self.y_unit = 1 * y_unit.u

# def __call__(self, x: float | Quantity, nu=0, extrapolate=None):

#     # check x and update if needed
#     x = self._check_input(x, self.x_unit)

#     # run the interpolation
#     result = super().__call__(x, nu, extrapolate)

#     # add the y unit to the result
#     return result * self.y_unit

# def integrate(self, a, b, extrapolate=None):

#     # check x and update if needed
#     a = self._check_input(a, self.x_unit)
#     b = self._check_input(b, self.x_unit)

#     # compute the definite integral
#     def_integ = super().integrate(a, b, extrapolate)

#     # add the units
#     return def_integ * self.x_unit * self.y_unit

# def solve(self, y=0.0, discontinuity=True, extrapolate=None):

#     # check x and update if needed
#     y = self._check_input(y, self.y_unit)

#     roots = super().solve(y, discontinuity, extrapolate)

#     # add the units
#     return roots * self.y_unit

# def _check_input(self, other, int_unit: Quantity):

#     # check compatibility of other and internal unit
#     if not int_unit.is_compatible_with(other):
#         raise ValueError(
#             f"x is not compatible with the units of the interpolator ({int_unit.u})."
#         )

#     # convert other unit to interpolator base units (if exists)
#     # if float, then int_unit is non-dim (see is_compatible_with above)
#     if isinstance(other, Quantity):
#         other = other.m_as(int_unit)

#     return other


class InterpolatorWithUnits(PPolyWithUnits):
    """Interpolator with units.

    For most usecases, the `from_ipol_method` factory
    is more convenient than this constructor.

    Parameters
    ----------
    ppoly : type[PPoly]
        Interpolator subclassing `PPoly`
    x_unit : Quantity, optional
        unit of the x axis, by default None
    y_unit : Quantity, optional
        unit of the y axis, by default None
    """

    def __init__(
        self,
        ipol: type[PPoly],
        x_unit: Quantity,
        y_unit: Quantity,
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

        The interpolator is chosen with the enum
        `InterpolatorWithUnitTypes`, among scipy cubic and
        monotonic interpolators. Interpolator specific parameters
        can be passed via `*args, **kwargs`, with the exception
        of the `method` parameter in `Akima1DInterpolator`.
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

        # strip units from x and y while extracting output units
        (x,), out_wrap_x = unwrap_and_wrap_consistent_units(x)
        (y,), out_wrap_y = unwrap_and_wrap_consistent_units(y)

        if ipol_type == InterpolatorWithUnitTypes.AKIMA:
            ipol = klass(x, y, axis=axis, extrapolate=extrapolate, method="akima")
        elif ipol_type == InterpolatorWithUnitTypes.MAKIMA:
            ipol = klass(x, y, axis=axis, extrapolate=extrapolate, method="makima")
        else:
            ipol = klass(x, y, *args, **kwargs)

        # out_wrap assigns units properly to x and y
        return InterpolatorWithUnits(ipol, out_wrap_x(1), out_wrap_y(1))
