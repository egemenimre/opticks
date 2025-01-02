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
    """Minimum sample size for each atomic interval for interpolation
    and resampling."""

    _DEFAULT_SAMPLE_SIZE = 100
    """Default sample size for each atomic interval for interpolation
    and resampling."""

    ipol_type: InterpolatorWithUnitTypes = InterpolatorWithUnitTypes.AKIMA
    """Interpolator type for resampling."""

    def __init__(self, mapping_or_iterable=None, sample_size: int = None):
        """Generates a new IntervalData object.

        If no argument is given, an empty `IntervalData` is created.
        If an argument is given, and is a mapping object (e.g., another
        `IntervalData`), an new IntervalData with the same key-value pairs
        (as well as combination method and sample size) is created. If an
        iterable is provided, it has to be a list of (key, value) pairs.

        The `sample_size` parameter sets the default sample size to be used
        for each atomic interval for resampling when needed.
        If `None` is given, the one from the `mapping_or_iterable` is used.
        If it does not exist, the default value is used.

        Parameters
        ----------
        mapping_or_iterable
            optional mapping or iterable
        sample_size : int
            default sample size for resampling
        """
        super().__init__(mapping_or_iterable)

        if mapping_or_iterable is not None and isinstance(
            mapping_or_iterable, IntervalData
        ):
            # copy the combination method
            self.combination_method = mapping_or_iterable.combination_method

            # assign or copy the sample size
            self.sample_size = mapping_or_iterable.sample_size

            # copy the interpolation type
            self.ipol_type = mapping_or_iterable.ipol_type

            # copy the resampled flag
            self._is_resampled = mapping_or_iterable._is_resampled
        else:
            # set the default combination method
            self.combination_method = FunctCombinationMethod.UNDEFINED

            # set the sample size
            if sample_size is not None:
                self.sample_size = sample_size
            else:
                self.sample_size = IntervalData._DEFAULT_SAMPLE_SIZE

            # set the default resampled flag
            self._is_resampled = False

    def copy_properties_to(self, other: "IntervalData") -> "IntervalData":
        """Copies the properties (except for the `IntervalDict`) of
        `Self` to `other`.

        Returns `other` for convenience."""

        other.ipol_type = self.ipol_type
        other.combination_method = self.combination_method

        other.sample_size = self.sample_size

        # this can be risky, therefore turned off
        # the user can cll this anytime, even after a combine operation.
        # other._is_resampled = self._is_resampled

        return other

    @classmethod
    def from_math_funct(
        cls, math_funct, valid_interval: P.Interval, sample_size: int = None
    ) -> "IntervalData":
        """Initialises the `IntervalData` object with a math function.

        The math function can be any callable, with a signature
        `y=math_funct(x)`.

        The `sample_size` parameter sets the default sample size to be used
        for resampling when needed. If `None` is given, the default
        value is used.

        Parameters
        ----------
        funct
            Math function with a signature `y=math_funct(x)`.
        valid_interval : Interval
            Interval of validity or domain
        sample_size : int
            default sample size for resampling

        Returns
        -------
        IntervalData
            The new `IntervalData` object
        """

        # init IntervalDict and assign the interval and the interpolator
        data = P.IntervalDict()

        data[valid_interval] = math_funct

        return IntervalData(data, sample_size=sample_size)

    @classmethod
    def from_interpolator(
        cls, ipol: type[PPolyWithUnits], sample_size: int = None
    ) -> "IntervalData":
        """Initialises the `IntervalData` object with an interpolator.

        The interval of validity is the data range of the interpolator.

        The `sample_size` parameter sets the default sample size to be used
        for resampling when needed. If `None` is given, the default
        value is used.

        Parameters
        ----------
        ipol : type[PPolyWithUnits]
            Interpolator
        sample_size : int
            default sample size for resampling

        Returns
        -------
        IntervalData
            The new `IntervalData` object
        """

        valid_interval = P.closed(ipol.x[0] * ipol.x_unit, ipol.x[-1] * ipol.x_unit)

        # create the IntervalData object
        intervalData = IntervalData.from_math_funct(
            ipol, valid_interval, sample_size=sample_size
        )

        # set the resampled flag
        intervalData._is_resampled = True

        return intervalData

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

    @property
    def sample_size(self) -> int:
        """Sample size for each atomic interval during any resampling needed."""
        return self._sample_size

    @sample_size.setter
    def sample_size(self, sample_size: int) -> Self:

        if sample_size >= IntervalData._MIN_SAMPLE_SIZE:
            self._sample_size = sample_size
        else:
            raise ValueError(
                f"Requsted sample size of {sample_size} is less than "
                f"the allowable minimum sample size of {IntervalData._MIN_SAMPLE_SIZE}."
            )

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

        If `x` is not covered by the bounds of the `Interval`,
        the value in the `default` parameter is returned, which is
        by default `None`.

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

    def get_value(
        self,
        x: Quantity | float,
        default=None,
    ) -> Quantity | float:
        """Gets the value at the arbitrary data point `x`.

        Computes and combines all the values of the mathematical functions
        at `x`. The combination method can be a summation or a multiplication,
        depending on the internal `combination_method` parameter. If the
        `combination_method` is set to neither and there are multiple values
        to be combined, a `ValueError` is thrown, as the method does not
        know how to handle the combination.

        If `x` is not covered by the bounds of the `Interval`,
        the value in the `default` parameter is returned, which is
        by default `None`.

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

    def combine(self, other: "IntervalData", *, missing=None) -> "IntervalData":
        """Combines this `IntervalData` object with another one.

        The combination method stacks the values or functions in the
        intersecting regions, but does not evaluate them, until
        `get_value` is called.

        The combination is not possible when one `IntervalDict` is
        of Summation type and the other Multiplication type. Then
        one of them should be resampled, resetting the combination
        method, which then enables a combination between the two.

        Similarly, two Undefined types cannot be combined as the
        `get_value` method on the combined object would not know how
        to combine the various functions.

        The properties of the new object is reset, but the user can
        easily call `copy_properties_to` from one of the combined
        objects, so that its properties are exported to the new object.

        The `missing` key is used for the values that are not
        intersecting between the two domains.

        Parameters
        ----------
        other : IntervalData
            other object to be combined
        missing : _type_, optional
            if set, use this value for missing values

        Returns
        -------
        IntervalData
            the new object with properties reset

        """

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
                f"Cannot combine conflicting IntervalData ({self.combination_method} "
                f"and {other.combination_method}). Resample at least one of them."
            )

        # overrides combine but the how function is fixed to append functions
        combined = super().combine(other, how=_append_math_functions, missing=missing)

        # set the combination method
        combined.combination_method = combination_method

        # combined can no longer be of type "resampled"
        combined._is_resampled = False

        return combined

    def scale(
        self,
        scaling_value: float | Quantity,
        *,
        missing=None,
        **kwargs,
    ) -> "IntervalData":
        """Scales the `IntervalData` object with a scalar.

        The scalar may also be with units.

        If this `IntervalData` object is of a Summation type, then the
        it is resampled using the internal `sample_size` and `interpolator`
        properties.

        The properties of this object is preserved in the returned object.

        Parameters
        ----------
        scaling_value : float | Quantity
            scalar value used in scaling
        missing
            if set, use this value for missing values when calling "how", by default None
        kwargs
            Parameters to be passed on the interpolator

        Returns
        -------
        IntervalData
            the new, scaled IntervalData object
        """

        # if self is Summation, resample self
        if self.combination_method == FunctCombinationMethod.SUM:
            to_be_scaled = self.resample(**kwargs)
        else:
            to_be_scaled = self

        # generate the new IntervalData for scaling
        data = P.IntervalDict()

        validity = self.domain()
        data[validity] = scaling_value

        scale = IntervalData(data)
        scale.combination_method = FunctCombinationMethod.MULTIPLY

        # combine with self via Multiply
        scaled = scale.combine(to_be_scaled, missing=missing)

        # copy the params of self to the scaled object
        scaled = self.copy_properties_to(scaled)

        # but make sure to set the combination method to Multiplication
        scaled.combination_method = FunctCombinationMethod.MULTIPLY

        return scaled

    def resample(
        self,
        **kwargs,
    ) -> "IntervalData":
        """Combines the stacked values and interpolators, via
        resampling and setting up a new interpolator.

        The process decomposes the intervals into atomic intervals
        and therefore the interval structure is changed.

        This analyses each interval enclosure within the dictionary
        and combines the functions / values. A `None` or zero
        takes precedence, and all other functions (continuous or
        interpolator) are resampled to create a new interpolator.

        The combination method can be a summation or a multiplication,
        depending on the `combination_method` parameter.

        The properties of this object is preserved in the returned object,
        but the combination properties is set to Undefined.

        The interpolator is defined by the `ipol_type` parameter,
        with the `extrapolate` already set to `True`. The `kwargs` are
        passed on to the interpolator definition.

        For the interpolator, each atomic interval is divided into
        an integer number of steps (given by the `sample_size` parameter)
        between its respective bounds.

        Parameters
        ----------
        kwargs
            Parameters to be passed on the interpolator

        Returns
        -------
        IntervalData
            The new `IntervalData` object
        """

        # decompose to atomic intervals and corresponding functions
        atomic_tuples = self.as_atomic()

        flattened = IntervalData()

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
                if len(functs) == 1:
                    # only a single item is present, just return it
                    result = functs[0]
                elif self.combination_method == FunctCombinationMethod.MULTIPLY:
                    result = math.prod(functs)
                elif self.combination_method == FunctCombinationMethod.SUM:
                    result = sum(functs)
                else:
                    raise ValueError(
                        f"Invalid combination method: {self.combination_method}. "
                        f"Set the combination method to 'SUM' or 'MULTIPLY'."
                    )
            else:

                # create the new resampling
                x_values, y_values = _generate_samples(
                    self.combination_method, interval, functs, self.sample_size
                )

                # init interpolator
                result = InterpolatorWithUnits.from_ipol_method(
                    self.ipol_type, x_values, y_values, extrapolate=True, **kwargs
                )

            # write result to the new IntervalData
            flattened[interval] = result

        # copy the params of self to the resampled object
        flattened = self.copy_properties_to(flattened)

        # but make sure to set the combination method to Multiplication
        flattened.combination_method = FunctCombinationMethod.UNDEFINED

        # flattened is of type "resampled"
        flattened._is_resampled = True

        return flattened

    def plot(self) -> "IntervalDataPlot":  # pragma: no cover
        """Convenience method to plot `IntervalData` objects.

        Returns an `IntervalDataPlot` object. The `set_plot_style`
        method can be invoked for further styling options and also
        the usual matplotlib `plot.ax` and `plot.fig` options are
        available for advanced customisation."""

        interval_data_dict = {
            "self": self,
        }

        plot = IntervalDataPlot(interval_data_dict, apply_default_style=False)

        plot.set_plot_style(legend_off=True)

        return plot

    def integrate(self, interval: P.Interval = None) -> float | Quantity:
        """Integrates the `IntervalData` over a certain interval.

        The integration can be over a user requested interval
        (with `interval` parameter defined) or can be over the full
        interval of the `IntervalData` object.

        Parameters
        ----------
        interval : P.Interval, optional
            integration interval

        Returns
        -------
        float | Quantity
            value of the integration over the interval
        """

        if interval is None:
            # integrate over the full interval
            data = self
        else:
            data = self.get(interval)

        if not data._is_resampled:
            # resample if needed
            data = data.resample()

        sum = 0

        # go through the atomic intervals and integrate
        for interval, funct in data.as_dict().items():

            interval: P.Interval

            if isinstance(funct, numbers.Number) or isinstance(funct, Quantity):
                # funct is constant value, just multiply
                # with the interval length
                sum += funct * (interval.upper - interval.lower)

            else:
                # Interpolated function, integrate the range

                funct: PPolyWithUnits

                sum += funct.integrate(interval.lower, interval.upper)

        return sum


def _generate_samples(
    combination_method, interval, functs, samples
) -> InterpolatorWithUnits:
    """Creates a resampling from the interval"""

    # generate range samples (x axis)
    x_values = np.linspace(
        interval.lower,
        interval.upper,
        num=samples,
        endpoint=True,
    )

    # generate values
    if len(functs) == 1:
        # only a single item is present, multiplication or summation is equivalent
        y_values = [_eval_functs_multiply(x, functs) for x in x_values]
    elif combination_method == FunctCombinationMethod.MULTIPLY:
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

    if fx is None or fy is None:
        # at least one is None
        return None
    else:
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


class IntervalDataPlot:  # pragma: no cover
    """Convenience class to plot `IntervalData` objects.

    Each `IntervalData` is used to generate the plot y values,
    using the dict key as the label. The optional `plot_interval`
    constrains the plot interval.

    Uses a default combination method of multiplication for
    combined interval data, though a summation can also be used.

    The plot samples each `IntervalData` object separately,
    therefore the sample points may not exactly match.
    If the points are required to match, then a `plot_interval`
    should be specified.

    Parameters
    ----------
    interval_data_dict : dict[str, IntervalData]
        dict of label and IntervalData objects
        plot_interval : Interval, optional
    plot interval, `None` means the entire domain is plotted
        for each `IntervalData` object
    apply_default_style : bool
        applies the default style
    """

    def __init__(
        self,
        interval_data_dict: dict[str, IntervalData],
        plot_interval: P.Interval = None,
        apply_default_style: bool = True,
    ) -> None:

        self.fig, self.ax = plt.subplots()

        self._populate_plot(interval_data_dict, plot_interval)

        # set a sensible default plot style
        if apply_default_style:
            self.set_plot_style()

    def _populate_plot(
        self,
        interval_data_dict: dict[str, IntervalData],
        plot_interval: P.Interval = None,
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

        Parameters
        ----------
        interval_data_dict : dict[str, IntervalData]
            list of IntervalData and labels
        plot_interval : Interval, optional
            plot interval, `None` means the entire domain is plotted
            for each `IntervalData` object
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
                        functs, interval, interval_data.sample_size
                    )

                    x_values.extend(x_list)
                    y_values.extend(y_list)

            elif interval_data.combination_method == FunctCombinationMethod.SUM:
                # combination method: summation
                for interval, functs in atomic_tuples:

                    x_list, y_list = self._prep_summed_data(
                        functs, interval, interval_data.sample_size
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
                        functs, interval, interval_data.sample_size
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
        legend_off: bool = False,
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
        legend_off : bool, optional
            turns the legend off, by default False

        """

        # set decorators
        if xlabel:
            self.ax.set_xlabel(xlabel)
        if ylabel:
            self.ax.set_ylabel(ylabel)
        if title:
            self.ax.set_title(title)

        if legend_off is False:
            self.fig.legend(bbox_to_anchor=(1, 1), loc="upper left")
        else:
            self.fig.legend().remove()

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

    def _prep_summed_data(self, functs, interval, sample_size):
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

            # generate range samples (x axis)
            x_list = np.linspace(
                interval.lower,
                interval.upper,
                num=sample_size,
                endpoint=True,
            )

            # generate values
            y_list = [_eval_functs_sum(x, functs) for x in x_list]

        return x_list, y_list

    def _prep_mult_data(self, functs, interval, sample_size):
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

            # generate range samples (x axis)
            x_list = np.linspace(
                interval.lower,
                interval.upper,
                num=sample_size,
                endpoint=True,
            )

            # generate values
            y_list = [_eval_functs_multiply(x, functs) for x in x_list]

        return x_list, y_list
