# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import numbers

import numpy as np
import portion as P

from opticks.utils.interval_data import (
    FunctCombinationMethod,
    IntervalData,
    IntervalDataPlot,
)
from opticks.utils.math_utils import PPolyWithUnits


class OpticalMaterial:
    """Generic Optical Material Properties class.

    Optical parameters are Reflectivity, Transmissivity and
    Emissivity, defined within the same interval of wavelengths.
    They are all defined between 0 and 1.

    The constructor runs internal checks to ensure the physical properties
    below are respected. If the sanity check fails, a `ValueError` will be
    thrown, explaining which interval has an issue. The `abs_tolerance`
    value is used to offset 1 and 0 to sidestep numerical issues. Sanity
    checks can be bypassed via the `skip_checks` flag, though it is not
    recommended.

    By definition, incident light can only  be transmitted,
    absorbed or reflected. Therefore, for any given wavelength,
    the sum of transmissivity, emissivity (or absorptivity) and
    reflectivity is equal to 1.

    Therefore: T+R+A=1

    Furthermore, according to Kirchhoff's Law of thermal
    radiation, for any given wavelength, emissivity is equal to
    absorptivity.

    Therefore: A = E

    Parameters
    ----------
    reflectivity : IntervalData
        Reflectivity defined in an interval
    transmissivity : IntervalData
        Transmissivity defined in an interval
    emissivity : IntervalData
        Emissivity (or absorptivity) defined in an interval
    name : str
        Name of the material
    abs_tolerance : float
        Absolute tolerance value for sanity checks
    skip_checks : bool
        SKip any sanity checks
    """

    name: str = None

    reflectivity: IntervalData = None

    transmissivity: IntervalData = None

    emissivity: IntervalData = None

    def __init__(
        self,
        reflectivity: IntervalData,
        transmissivity: IntervalData,
        emissivity: IntervalData,
        name="Generic Optical Material",
        abs_tolerance=1e-9,
        skip_checks: bool = False,
    ):

        self.reflectivity = reflectivity
        self.transmissivity = transmissivity
        self.emissivity = emissivity
        self.name = name

        # sanity check of the properties
        if not skip_checks:
            self._sanity_check(abs_tolerance)

    def _sanity_check(self, abs_tolerance):
        """Checks the properties for adherence to the laws of physics.

        Two tests are conducted:
        1. Is each property staying within 0 and 1 (inclusive)?
        2. Is the sum of all properties staying within 0 and 1 (inclusive)?

        To this end, the properties are resampled (to get to the polynomials)
        and their values at the edges as well as values at maxima and minima
        are checked.

        The method returns nothing but will throw a `ValueError` if any of
        the checks fail, explaining which interval has an issue.
        """

        t = self.transmissivity.resample()
        r = self.reflectivity.resample()
        e = self.emissivity.resample()

        # for every wavelength verify that 1 >= T or R or E >= 0
        # upon failure these will raise a ValueError
        _property_sanity_check(t, 1.0, 0, abs_tolerance, resampled=True)
        _property_sanity_check(r, 1.0, 0, abs_tolerance, resampled=True)
        _property_sanity_check(e, 1.0, 0, abs_tolerance, resampled=True)

        # for every wavelength verify that T+R+E=1
        t.combination_method = FunctCombinationMethod.SUM
        r.combination_method = FunctCombinationMethod.SUM
        e.combination_method = FunctCombinationMethod.SUM

        summed = t.combine(r).combine(e).resample()

        # self.summed = summed

        # upon failure this will raise a ValueError
        _property_sanity_check(summed, 1.0, 1.0, abs_tolerance, resampled=True)

    @classmethod
    def init_opaque_from_refl(
        cls,
        reflectivity: IntervalData,
        name="Generic Opaque",
    ) -> "OpticalMaterial":
        """Generates an "opaque" material from Reflectance only.

        Absorptivity or Emissivity is derived from the Reflectance as
        $A = 1 - R$. Transmissivity is by definition zero.

        Parameters
        ----------
        reflectivity : IntervalData
            Reflectivity within an interval
        name : str
            Name of the matierial.

        Returns
        -------
        OpticalMaterial
            generated object
        """

        domain = reflectivity.domain()

        # transmissivity is zero
        transmissivity = IntervalData({domain: 0.0})

        # emissivity = 1-reflectivity
        # resampling is needed for the summation
        neg_ref = reflectivity.scale(-1.0).resample()

        unity = IntervalData({domain: 1.0})
        unity.combination_method = FunctCombinationMethod.SUM

        # generate emissivity
        emissivity = unity.combine(neg_ref)
        # copy properties like interpolator or sample size
        # but not the combination method
        emissivity = reflectivity.copy_properties_to(emissivity)
        emissivity.combination_method = FunctCombinationMethod.SUM

        # resampling is needed to reset the combination method
        # and to initialise the interpolators
        emissivity = emissivity.resample()

        return OpticalMaterial(
            reflectivity=reflectivity,
            transmissivity=transmissivity,
            emissivity=emissivity,
            name=name,
        )

    @classmethod
    def init_blackbody(
        cls, domain: P.Interval, name="Generic Blackbody"
    ) -> "OpticalMaterial":

        emissivity = IntervalData({domain: 1.0})
        transmissivity = IntervalData({domain: 0.0})
        reflectivity = IntervalData({domain: 0.0})

        # sanity checks can be theoretically skipped
        # but sticking to them is good practice

        return OpticalMaterial(
            reflectivity=reflectivity,
            transmissivity=transmissivity,
            emissivity=emissivity,
            name=name,
        )

    @property
    def absorptivity(self) -> IntervalData:
        """Absorptivity, which is equivalent to the Emissivity."""
        return self.emissivity

    def plot(self) -> IntervalDataPlot:  # pragma: no cover
        """Convenience method to plot `OpticalMaterial` objects.

        Returns an `IntervalDataPlot` object. The `set_plot_style`
        method can be invoked for further styling options and also
        the usual matplotlib `plot.ax` and `plot.fig` options are
        available for advanced customisation."""

        interval_data_dict = {
            "reflectivity": self.reflectivity,
            "absorptivity": self.absorptivity,
            "transmissivity": self.transmissivity,
            # "summed": self.summed,
        }

        plot = IntervalDataPlot(interval_data_dict, apply_default_style=False)

        plot.set_plot_style(title=f"{self.name} Material Optical Properties")

        return plot


def _check_validity(value, max, min, abs_tolerance):
    """Validity check function.

    Tolerances are important to make sure we don't fail on edge cases."""
    return max + abs_tolerance >= value >= min - abs_tolerance


def _property_sanity_check(
    property: IntervalData, max, min, abs_tolerance, resampled=False
):
    """Checks the property against max and min values.

    Each interval in the property is checked against max and min values
    at the edges, as well as at the maxima/minima to ensure that the values
    are within the limits.

    The method returns nothing but will throw a `ValueError` if any of
    the checks fail, explaining which interval has an issue.

    Parameters
    ----------
    property : IntervalData
        property to be tested
    max : float
        max value for the sanity check
    min : float
        min value for the sanity check
    abs_tolerance : float
        absolute tolerance value to offset max and min
    resampled : bool
        True if the `property` is already resampled, false otherwise
    """

    # resample to flatten the data and ensure the polynomials are initialised.
    # This is important for derivatives and roots later.
    if resampled:
        # data is already resampled
        data = property
    else:
        # data not resampled, resample it
        data = property.resample()

    for interval, funct in data.as_dict().items():

        interval: P.Interval

        err_flag = False
        err_text = ""

        if isinstance(funct, numbers.Number):
            # funct is constant value, just check the validity

            if not _check_validity(funct, max, min, abs_tolerance):
                err_flag = True
                err_text += f"The value of {funct} within the interval {interval} is out of bounds."

        else:
            # Interpolated function, check the entire range

            funct: PPolyWithUnits

            # Check whether all data lies between 0 and 1
            # If linear, the data at the ends should be within
            # the valid range, and that's it.
            #
            # If second order, the data at the ends and
            # the maxima/minima should be within the valid range.
            # Extrapolating the logic, check the data at the ends
            # *and* the maxima/minima.

            # check the data at the edges
            left_edge_value = funct(interval.lower)
            right_edge_value = funct(interval.upper)

            # print("values at edges: ", left_edge_value, right_edge_value)

            if not _check_validity(left_edge_value, max, min, abs_tolerance):
                err_flag = True
                err_text += (
                    f"The left edge of the interval {interval} is out of bounds."
                )

            if not _check_validity(right_edge_value, max, min, abs_tolerance):
                err_flag = True
                err_text += (
                    f"The right edge of the interval {interval} is out of bounds."
                )

            # find the maxima/minima and evaluate the values there

            # max/min points are the roots of the derivative
            # hence one less than the order of the polynomial

            deriv = funct.derivative()
            maxmin_pts = deriv.roots(extrapolate=False)

            # # check whether
            # x = np.linspace(interval.lower, interval.upper, num=100, endpoint=True)
            # y = [deriv(x_val).m for x_val in x if x_val is not np.nan]

            # print("shapes ", np.shape(x), np.shape(y))

            # print("deriv values", y)
            # print("corr coef", np.corrcoef(x.m, y))

            # delete the NaNs in the roots
            # this happens when the derivative is equal to zero for a
            # certain part of the interval.
            # See scipy.PPoly.solve() for more info.
            maxmin_pts = [
                maxmin_pt for maxmin_pt in maxmin_pts if not np.isnan(maxmin_pt)
            ]

            # check whether the function values at all maxima / minima
            # are within the bounds.

            # maxmin_values = [funct(maxmin_pt) for maxmin_pt in maxmin_pts]
            # print("maxmin pts ", maxmin_pts)
            # print("maxmin val ", maxmin_values)

            maxmin_values_are_ok = all(
                [
                    _check_validity(funct(maxmin_pt), max, min, abs_tolerance)
                    for maxmin_pt in maxmin_pts
                ]
            )

            if not maxmin_values_are_ok:
                err_flag = True
                err_text += (
                    f"At least one of the maxima / minima within the interval {interval} "
                    "is out of bounds."
                )

        if err_flag:
            # error flag is raised
            raise ValueError(
                f"{err_text} The values within the interval should stay between 0 and 1 (inclusive)."
            )
