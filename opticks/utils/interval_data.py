# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


import math
import numbers
from enum import Enum
from typing import Any, Self

import numpy as np
import portion as P
from matplotlib import pyplot as plt
from pint import Quantity

from opticks import Q_, u
from opticks.utils.math_utils import (
    InterpolatorWithUnits,
    InterpolatorWithUnitTypes,
    PPolyWithUnits,
)

FunctCombinationMethod = Enum(
    "FunctCombinationMethod",
    [("MULTIPLY", "Multiplication"), ("SUM", "Summation"), ("UNDEFINED", "Undefined")],
)
"""Interval Data Combination Method Enum."""


class IntervalData(P.IntervalDict):

    _MIN_SAMPLE_SIZE = 20
    """Minimum sample size for interpolation and resampling."""

    def __init__(self, mapping_or_iterable=None):
        super().__init__(mapping_or_iterable)

        # copy the combination method if exists
        if mapping_or_iterable is not None and isinstance(
            mapping_or_iterable, IntervalData
        ):
            self._combination_method = mapping_or_iterable.combination_method
        else:
            self._combination_method = FunctCombinationMethod.UNDEFINED

    @classmethod
    def from_math_funct(cls, math_funct, valid_interval: P.Interval) -> "IntervalData":
        """Initialises the `IntervalData` object with a math function.

        The math function can be any callable, with a signature
        `y=math_funct(x)`.

        Parameters
        ----------
        funct
            Math function with a signature `y=math_funct(x)`.
        valid_interval : Interval
            Interval of validity or domain

        Returns
        -------
        IntervalData
            The new `IntervalData` object
        """

        # init IntervalDict and assign the interval and the interpolator
        data = P.IntervalDict()

        data[valid_interval] = math_funct

        return IntervalData(data)

    @classmethod
    def from_interpolator(cls, ipol: type[PPolyWithUnits]) -> "IntervalData":
        """Initialises the `IntervalData` object with an interpolator.

        The interval of validity is the data range of the interpolator.

        Parameters
        ----------
        ipol : type[PPolyWithUnits]
            Interpolator

        Returns
        -------
        IntervalData
            The new `IntervalData` object
        """

        valid_interval = P.closed(ipol.x[0] * ipol.x_unit, ipol.x[-1] * ipol.x_unit)

        return IntervalData.from_math_funct(ipol, valid_interval)

    def as_atomic(self) -> list[tuple[P.Interval, Any]]:
        """Returns a list of tuples containing atomic intervals
        and corresponding functions."""

        return [(ival, functs) for ivals, functs in self.items() for ival in ivals]

    @property
    def combination_method(self) -> FunctCombinationMethod:
        """Combination method for the stack of functions."""
        return self._combination_method

    @combination_method.setter
    def combination_method(self, combination_method: FunctCombinationMethod) -> Self:

        if combination_method in FunctCombinationMethod:
            self._combination_method = combination_method
        else:
            raise ValueError(
                f"Requested combination method ({combination_method}) is invalid."
            )

        return self

    def get(
        self, x: Quantity | P.Interval, default=None
    ) -> type[P.IntervalDict] | Quantity | float:
        """Gets the functions at the arbitrary data point `x`.

        This does not return the numerical result at `x`, rather it returns the
        list of functions (interpolators, mathematical functions or fixed values)
        at `x`. If the mathematical combination (summation or multiplication)
        of all the functions at `x` is desired, use the `get_value` method.

        If `x` is a single value, it returns a single value (if it exists) or
        `None` (if it doesn't).

        If `x` is an `Interval`, it returns a new `IntervalDict`
        restricted to the requested interval. In that case, the `default` value
        is used to "fill the gaps" (if any) w.r.t. given `x`.

        If `x` is not covered by the bounds of the internal `Interval`,
        then returns `None`.


        Parameters
        ----------
        x : Quantity | Interval
            Requested data point or Interval
        default
            default value to be used, by default None

        Returns
        -------
        IntervalDict | Quantity | float
            an IntervalDict/IntervalData, or a single value if key is not an Interval.
        """

        return super().get(x, default=default)

    def combine(self, other: "IntervalData", *, missing=None) -> Self:

        # check whether the two IntervalDicts can be combined.
        # They can be combined when either the combination methods
        # are the same or at least one is "UNDEFINED".
        if (
            self.combination_method
            == other.combination_method
            == FunctCombinationMethod.UNDEFINED
        ):
            raise ValueError(
                "Cannot combine UNDEFINED IntervalDicts. Set one to MULTIPLY or SUM."
            )
        if self.combination_method == other.combination_method:
            combination_method = self.combination_method
        elif self.combination_method == FunctCombinationMethod.UNDEFINED:
            combination_method = other.combination_method
        elif other.combination_method == FunctCombinationMethod.UNDEFINED:
            combination_method = self.combination_method
        else:
            raise ValueError(
                f"Cannot combine conflicting IntervalDicts ({self.combination_method}"
                f"and {other.combination_method}). Resample at least one of them."
            )

        # overrides combine but the how function is fixed to append functions
        combined = super().combine(other, how=_append_math_functions, missing=missing)

        # set the combination method
        combined.combination_method = combination_method

        return combined

    def get_value(
        self,
        x: Quantity | float,
        default=None,
    ) -> Quantity | float:
        """Gets the value at the arbitrary data point `x`.

        Computes and combines all the values of the mathematical functions
        at `x`. The combination method can be a summation or a multiplication,
        depending on the `combination_method` parameter. If the
        `combination_method` is set to neither and there are multiple values
        to be combined, a `ValueError` is thrown, as the method does not
        know how to handle the combination.

        In case `x` does not exist, the value in the `default` parameter is
        returned, which is by default `None`.

        Parameters
        ----------
        x : Quantity | float
            Requested data point
        default
            default value to be used, by default None

        Returns
        -------
        Quantity | float
            Value at the requested data point
        """
        # Retrieve all math functions at x
        functs = super().get(x, default=default)

        if isinstance(functs, list):
            # there are multiple values

            # combine them as requested
            if self.combination_method == FunctCombinationMethod.MULTIPLY:
                result = _eval_functs_multiply(x, functs)
            elif self.combination_method == FunctCombinationMethod.SUM:
                result = _eval_functs_sum(x, functs)
            else:
                raise ValueError(
                    f"Invalid combination method: {self.combination_method}. "
                    f"Set the combination method to 'SUM' or 'MULTIPLY'."
                )
        else:
            # there is a single value, just return it
            result = _eval_functs_sum(x, functs)

        return result

    def resample(
        self,
        approx_stepsize: float | Quantity,
        ipol_type=InterpolatorWithUnitTypes.AKIMA,
        **kwargs,
    ) -> "IntervalData":
        """Combines the stacked values and interpolators, via
        resampling and setting up a new interpolator.

        The process decomposes the intervals into atomic
        intervals and therefore the interval structure is
        changed.

        This analyses each interval enclosure within the dictionary
        and combines the functions / values. A `None` or zero
        takes precedence, and all other functions (continuous or
        interpolator) are resampled to create a new interpolator.

        The combination method can be a summation or a multiplication,
        depending on the `combination_method` parameter.

        The stepsize for this resampling is approximate, in the
        sense that, each atomic interval is divided into an integer
        number of steps between its respective bounds.

        The interpolator is defined by the user, with the `extrapolate`
        already set to `True`. The `kwargs` are passed on to the
        interpolator definition.

        Parameters
        ----------
        approx_stepsize : float | Quantity
            approximate stepsize
        ipol_type : InterpolatorWithUnitTypes
            Interpolator type, defaults to Akima
        combination_method : FunctCombinationMethod
            Function combination method (multiplication or summation)

        Returns
        -------
        IntervalData
            The new `IntervalData` object
        """

        # decompose to atomic intervals and corresponding functions
        atomic_tuples = self.as_atomic()

        new_intdict = IntervalData()

        for interval, functs in atomic_tuples:

            # force functs into a list
            if not isinstance(functs, list):
                functs = [functs]

            if any(funct is None for funct in functs):
                # at least one funct is None
                result = None

            elif self.combination_method == FunctCombinationMethod.MULTIPLY and any(
                funct == 0 for funct in functs
            ):
                # at least one funct is zero
                result = 0

            elif all(
                (isinstance(funct, numbers.Number) or isinstance(funct, Quantity))
                for funct in functs
            ):

                # all functs are numbers
                if self.combination_method == FunctCombinationMethod.MULTIPLY:
                    result = math.prod(functs)
                elif self.combination_method == FunctCombinationMethod.SUM:
                    result = sum(functs)
                else:
                    raise ValueError(
                        f"Invalid combination method: {self.combination_method}. "
                        f"Set the combination method to 'SUM' or 'MULTIPLY'."
                    )
            else:
                # create a new resampling
                x_values, y_values = _generate_samples(
                    self.combination_method,
                    interval,
                    functs,
                    approx_stepsize,
                    self._MIN_SAMPLE_SIZE,
                )

                # init interpolator
                result = InterpolatorWithUnits.from_ipol_method(
                    ipol_type, x_values, y_values, extrapolate=True, **kwargs
                )

            # write result to the new IntervalData
            new_intdict[interval] = result

        return new_intdict


def _generate_samples(
    combination_method,
    interval,
    functs,
    approx_stepsize,
    min_sample_size,
) -> InterpolatorWithUnits:
    """Creates a resampling from the interval"""

    # generate sample size
    samples = int((interval.upper - interval.lower) / approx_stepsize)

    if samples < min_sample_size:
        samples = min_sample_size

    # generate range samples (x axis)
    x_values = np.linspace(
        interval.lower,
        interval.upper,
        num=samples,
        endpoint=True,
    )

    # generate values
    if combination_method == FunctCombinationMethod.MULTIPLY:
        y_values = [_eval_functs_multiply(x, functs) for x in x_values]
    elif combination_method == FunctCombinationMethod.SUM:
        y_values = [_eval_functs_sum(x, functs) for x in x_values]
    else:
        raise ValueError(
            f"Invalid combination method: {combination_method}. "
            f"Set the combination method to 'SUM' or 'MULTIPLY'."
        )

    return x_values, y_values


def _eval_functs_multiply(x, functs) -> Quantity | float:
    """Evaluates `x` within a set of functions via multiplication.

    Computes and multiplies all the values of the mathematical functions
    at `x`.

    Note that `x` may be outside the domain of an interpolator,
    therefore this method should be used cautiously.

    Parameters
    ----------
    x : Quantity | float
        Requested data point
    functs : Any
        functions to evaluate

    Returns
    -------
    Quantity | float
        Value at the requested data point
    """

    # upgrade functs to list if only single item is present
    if not isinstance(functs, list):
        functs = [functs]

    result = 1.0

    # loop through the functions to multiply the values
    for funct in functs:

        if funct is None:
            return None
        elif isinstance(funct, numbers.Number) or isinstance(funct, Quantity):
            # int or float
            result *= funct
        else:
            # interpolator (or other callable)
            result *= funct(x)

    return result


def _eval_functs_sum(x, functs) -> Quantity | float:
    """Evaluates `x` within a set of functions via summation.

    Computes and sums all the values of the mathematical functions
    at `x`.

    Note that `x` may be outside the domain of an interpolator,
    therefore this method should be used cautiously.

    Parameters
    ----------
    x : Quantity | float
        Requested data point
    functs : Any
        functions to evaluate

    Returns
    -------
    Quantity | float
        Value at the requested data point
    """

    # upgrade functs to list if only single item is present
    if not isinstance(functs, list):
        functs = [functs]

    result = 0.0

    # loop through the functions to multiply the values
    for funct in functs:

        if funct is None:
            return None
        elif isinstance(funct, numbers.Number) or isinstance(funct, Quantity):
            # int or float
            result += funct
        else:
            # interpolator (or other callable)
            result += funct(x)

    return result


def _append_math_functions(fx, fy) -> list | None:
    """Appends the math functions into a list for
    `IntervalDict.combine`."""

    if fx or fy:
        # neither should be None

        # make sure the result is a list
        if isinstance(fx, list):
            result = fx
        else:
            result = [fx]

        # append fy to the result
        if not isinstance(fy, list):
            result.append(fy)
        else:
            result.extend(fy)

        return result
    else:
        # at least one is None
        return None


class IntervalDataPlot:  # pragma: no cover
    """Convenience class to plot `IntervalData` objects

    Each `IntervalData` is used to generate the plot y values,
    using the dict key as the label. The optional `plot_interval`
    constrains the plot interval.

    Uses a default combination method of multiplication for
    combined interval data, though a summation can also be used.

    The plot samples each `IntervalData` object separately,
    therefore the sample points may not exactly match.
    If the points are required to match, then a `plot_interval`
    should be specified.

    The stepsize for this resampling is approximate,
    in the sense that, each atomic interval is divided into
    an integer number of steps between its respective bounds.

    Parameters
    ----------
    interval_data_dict : dict[str, IntervalData]
        dict of label and IntervalData objects
        plot_interval : Interval, optional
    plot interval, `None` means the entire domain is plotted
        for each `IntervalData` object
    approx_stepsize : float | Quantity
        approximate stepsize, if None the minimum value is used

    """

    def __init__(
        self,
        interval_data_dict: dict[str, IntervalData],
        plot_interval: P.Interval = None,
        approx_stepsize: float | Quantity = None,
    ) -> None:

        self.fig, self.ax = plt.subplots()

        self._populate_plot(interval_data_dict, plot_interval, approx_stepsize)

        # set a sensible default plot style
        self.set_plot_style()

    def _populate_plot(
        self,
        interval_data_dict: dict[str, IntervalData],
        plot_interval: P.Interval = None,
        approx_stepsize: float | Quantity = None,
    ) -> None:
        """
        Populates the plot lines using the `IntervalData` objects.

        Each `IntervalData` is used to generate the plot y values,
        using the dict key as the label. The optional `plot_interval`
        constrains the plot interval.

        The combination method of for the combined interval data
        can be summation or multiplication, depending on whatever
        is specified in the `IntervalData` object.

        The plot samples each `IntervalData` object separately,
        therefore the sample points may not exactly match.
        If the points are required to match, then a `plot_interval`
        should be specified.

        The stepsize for this resampling is approximate,
        in the sense that, each atomic interval is divided into
        an integer number of steps between its respective bounds.

        Parameters
        ----------
        interval_data_dict : dict[str, IntervalData]
            list of IntervalData and labels
        plot_interval : Interval, optional
            plot interval, `None` means the entire domain is plotted
            for each `IntervalData` object
        approx_stepsize : float | Quantity
            approximate stepsize, if None the minimum value is used
        """

        # generate IntervalData lines
        for label, interval_data in interval_data_dict.items():

            # If plot interval is defined, intersect here
            # and get the new IntervalData
            if plot_interval is not None:
                interval_data = interval_data.get(plot_interval, default=None)

            x_values = []
            y_values = []

            atomic_tuples = interval_data.as_atomic()
            # this is to order the intervals
            atomic_tuples.sort(key=lambda tup: tup[0].upper)

            if interval_data.combination_method == FunctCombinationMethod.MULTIPLY:
                # combination method: multiplication
                for interval, functs in atomic_tuples:

                    x_list, y_list = self._prep_mult_data(
                        functs, interval, approx_stepsize
                    )

                    x_values.extend(x_list)
                    y_values.extend(y_list)

            elif interval_data.combination_method == FunctCombinationMethod.SUM:
                # combination method: summation
                for interval, functs in atomic_tuples:

                    x_list, y_list = self._prep_summed_data(
                        functs, interval, approx_stepsize
                    )

                    x_values.extend(x_list)
                    y_values.extend(y_list)
            else:
                # combination method undefined
                for interval, functs in atomic_tuples:

                    if isinstance(functs, list):
                        raise ValueError(
                            f"Invalid combination method: {interval_data.combination_method}. "
                            f"Set the combination method to 'SUM' or 'MULTIPLY'."
                        )

                    x_list, y_list = self._prep_summed_data(
                        functs, interval, approx_stepsize
                    )

                    x_values.extend(x_list)
                    y_values.extend(y_list)

            # print(label)
            # for x, y in zip(x_values, y_values):
            #     print(x, y)

            # convert from list of units to a vectorised form
            if isinstance(x_values[0], Quantity):
                x_values = Q_.from_list(x_values)

            if isinstance(y_values[0], Quantity):
                y_values = Q_.from_list(y_values)

            # generate plot line
            self.ax.plot(x_values, y_values, label=label)

    def set_plot_style(
        self,
        title: str = None,
        xlabel: str = None,
        ylabel: str = None,
        height: Quantity | float = 4 * u.inch,
        width: Quantity | float = 6 * u.inch,
    ) -> None:
        """
        Sets some basic style parameters for the plot.

        Parameters
        ----------
        title : str, optional
            title of the plot, by default None
        xlabel : str, optional
            x-axis label, by default None
        ylabel : str, optional
            y-axis label, by default None
        height : Quantity | float, optional
            height of the figure (in inches), by default 4 in
        width : Quantity | float, optional
            width of the figure (in inches), by default 6 in

        """

        # set decorators
        if xlabel:
            self.ax.set_xlabel(xlabel)
        if ylabel:
            self.ax.set_ylabel(ylabel)
        if title:
            self.ax.set_title(title)

        self.fig.legend()

        # set plot formatting
        self.ax.xaxis.grid(True)
        self.ax.yaxis.grid(True)

        self.fig.tight_layout()

        if isinstance(height, Quantity):
            height = height.m_as(u.inch)
        if isinstance(width, Quantity):
            width = width.m_as(u.inch)

        self.fig.set_figheight(height)
        self.fig.set_figwidth(width)

    # def show_plot(self):
    #     plt.show()

    def _prep_summed_data(self, functs, interval, approx_stepsize):
        """Prepares the data for summed combination"""

        if not isinstance(functs, list):
            functs = [functs]

        if any(funct is None for funct in functs):
            # at least one funct is None
            x_list = [interval.lower, interval.upper]
            y_list = [None, None]
        elif all(
            (isinstance(funct, numbers.Number) or isinstance(funct, Quantity))
            for funct in functs
        ):
            # all functs are numbers
            result = sum(functs)
            x_list = [interval.lower, interval.upper]
            y_list = [result, result]
        else:
            # create a new resampling

            # create sample size
            if approx_stepsize is None:
                samples = IntervalData._MIN_SAMPLE_SIZE
            else:
                # generate sample size
                samples = int((interval.upper - interval.lower) / approx_stepsize)

                if samples < IntervalData._MIN_SAMPLE_SIZE:
                    samples = IntervalData._MIN_SAMPLE_SIZE

            # generate range samples (x axis)
            x_list = np.linspace(
                interval.lower,
                interval.upper,
                num=samples,
                endpoint=True,
            )

            # generate values
            y_list = [_eval_functs_sum(x, functs) for x in x_list]

        return x_list, y_list

    def _prep_mult_data(self, functs, interval, approx_stepsize):
        """Prepares the data for multiplied combination"""

        if not isinstance(functs, list):
            functs = [functs]

        if any(funct is None for funct in functs):
            # at least one funct is None
            x_list = [interval.lower, interval.upper]
            y_list = [None, None]
        elif any(funct == 0 for funct in functs):
            # at least one funct is zero
            x_list = [interval.lower, interval.upper]
            y_list = [0, 0]
        elif all(
            (isinstance(funct, numbers.Number) or isinstance(funct, Quantity))
            for funct in functs
        ):
            # all functs are numbers
            result = math.prod(functs)
            x_list = [interval.lower, interval.upper]
            y_list = [result, result]
        else:
            # create a new resampling

            # create sample size
            if approx_stepsize is None:
                samples = IntervalData._MIN_SAMPLE_SIZE
            else:
                # generate sample size
                samples = int((interval.upper - interval.lower) / approx_stepsize)

                if samples < IntervalData._MIN_SAMPLE_SIZE:
                    samples = IntervalData._MIN_SAMPLE_SIZE

            # generate range samples (x axis)
            x_list = np.linspace(
                interval.lower,
                interval.upper,
                num=samples,
                endpoint=True,
            )

            # generate values
            y_list = [_eval_functs_multiply(x, functs) for x in x_list]

        return x_list, y_list
