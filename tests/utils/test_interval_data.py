# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import math
import numbers

import numpy as np
import portion as P
import pytest
from pint import Quantity

from opticks import u
from opticks.utils.interval_data import FunctCombinationMethod, IntervalData
from opticks.utils.math_utils import InterpolatorWithUnits, InterpolatorWithUnitTypes
from opticks.utils.testing_utils import assert_allclose


class TestIntervalData:

    @pytest.fixture(scope="class")
    def filter(self) -> IntervalData:
        data = P.IntervalDict()
        # interval of validity
        validity = P.closed(-10 * u.Hz, 10 * u.Hz)
        data[validity] = 0
        # data proper
        first_rng = P.closed(-4 * u.Hz, 0 * u.Hz)
        data[first_rng] = 1.0
        second_rng = P.closed(2 * u.Hz, 5 * u.Hz)
        data[second_rng] = 1.0

        filter = IntervalData(data)

        return filter

    @pytest.fixture(scope="class")
    def filter_no_units(self) -> IntervalData:
        data = P.IntervalDict()
        # interval of validity
        validity = P.closed(-10, 10)
        data[validity] = 0
        # data proper
        first_rng = P.closed(-4, 0)
        data[first_rng] = 1.0
        second_rng = P.closed(2, 5)
        data[second_rng] = 1.0

        filter = IntervalData(data)

        return filter

    @pytest.fixture(scope="class")
    def filter2(self) -> IntervalData:
        data = P.IntervalDict()
        # interval of validity
        validity = P.closed(-20 * u.Hz, 20 * u.Hz)
        data[validity] = 0.1
        # data proper
        first_rng = P.closed(-8 * u.Hz, -6 * u.Hz)
        data[first_rng] = 1.0
        second_rng = P.closed(3 * u.Hz, 7 * u.Hz)
        data[second_rng] = 1.0

        filter = IntervalData(data)

        return filter

    @pytest.fixture(scope="class")
    def main_funct(self) -> IntervalData:

        range = P.closed(-8 * u.Hz, 8 * u.Hz)

        x = np.linspace(range.lower, range.upper, num=100, endpoint=True)

        y = (0.5 * x.m) ** 2

        ipol = InterpolatorWithUnits.from_ipol_method(
            InterpolatorWithUnitTypes.AKIMA, x, y, extrapolate=True
        )

        main_funct = IntervalData.from_interpolator(ipol)

        return main_funct

    @pytest.fixture(scope="class")
    def main_funct_no_units(self) -> IntervalData:

        range = P.closed(-8, 8)

        x = np.linspace(range.lower, range.upper, num=100, endpoint=True)

        y = (0.5 * x) ** 2

        ipol = InterpolatorWithUnits.from_ipol_method(
            InterpolatorWithUnitTypes.AKIMA, x, y, extrapolate=True
        )

        main_funct = IntervalData.from_interpolator(ipol)

        return main_funct

    def test_edge(self, filter, main_funct):
        """Test the edge with multiplication."""

        # cases
        # -------
        x = [1.9999, 2.0, 2.0000001] * u.Hz

        # truth
        # -------
        truth = [0.0, 1.0, 1.00000010]

        # test set up
        # -----------------
        filter.combination_method = FunctCombinationMethod.MULTIPLY
        filtered_main = main_funct.combine(filter)

        results = [filtered_main.get_value(x_val) for x_val in x]

        # verification
        # -----------------
        np.testing.assert_allclose(results, truth, rtol=1e-10)

    def test_resampled_edge(self, filter, main_funct):
        """Test the resampled edge with multiplication."""

        # cases
        # -------
        x = [1.9999, 2.0, 2.0000001] * u.Hz

        # truth
        # -------
        truth = [0.0, 1.0, 1.00000010]

        # test set up
        # -----------------
        approx_stepsize = 0.15 * u.Hz

        filter.combination_method = FunctCombinationMethod.MULTIPLY
        filtered_main = main_funct.combine(filter)

        resampled_filt_main = filtered_main.resample(approx_stepsize)

        results = [resampled_filt_main.get_value(x_val) for x_val in x]

        # verification
        # -----------------
        np.testing.assert_allclose(results, truth, rtol=1e-8)

    def test_resampled_edge_no_units(self, filter_no_units, main_funct_no_units):
        """Test the resampled edge with multiplication - no units."""

        # cases
        # -------
        x = [1.9999, 2.0, 2.0000001]

        # truth
        # -------
        truth = [0.0, 1.0, 1.00000010]

        # test set up
        # -----------------
        approx_stepsize = 0.15

        filter_no_units.combination_method = FunctCombinationMethod.MULTIPLY
        filtered_main = main_funct_no_units.combine(filter_no_units)

        resampled_filt_main = filtered_main.resample(approx_stepsize)

        results = [resampled_filt_main.get_value(x_val) for x_val in x]

        # verification
        # -----------------
        np.testing.assert_allclose(results, truth, rtol=1e-8)

    def test_summation(self, filter, filter2):
        """Test the edge with summation."""

        # cases
        # -------
        x = [-3.0, 2.9999, 3.0, 4.5, 5.0, 5.00001] * u.Hz

        # truth
        # -------
        truth = [1.1, 1.1, 2.0, 2.0, 2.0, 1.0]

        # test set up
        # -----------------
        filter.combination_method = FunctCombinationMethod.SUM
        combined_filter = filter2.combine(filter)

        results = [combined_filter.get_value(x_val) for x_val in x]

        # verification
        # -----------------
        np.testing.assert_allclose(results, truth, rtol=1e-10)

    def test_number_check(self):
        """Checks the number & Quantity check logic with
        summation and multiplication."""

        # test summation

        functs1 = [10 * u.m, 20 * u.m, 30 * u.m]

        functs = functs1

        if not isinstance(functs, list):
            functs = [functs]

        # check whether all are numbers (including Quantity objects)
        number_check = all(
            (isinstance(funct, numbers.Number) or isinstance(funct, Quantity))
            for funct in functs
        )

        sum_functs = sum(functs)

        assert number_check

        assert_allclose(sum_functs, 60 * u.m, atol=1e-6 * u.mm)

        # test multiplication

        functs2 = [10 * u.Hz, 2, 40 * u.Hz]

        functs = functs2

        if not isinstance(functs, list):
            functs = [functs]

        # check whether all are numbers (including Quantity objects)
        number_check = all(
            (isinstance(funct, numbers.Number) or isinstance(funct, Quantity))
            for funct in functs
        )

        mult_functs = math.prod(functs)

        assert number_check

        assert_allclose(mult_functs, 800 * u.Hz**2, atol=1e-6 * u.Hz**2)

    def test_combine_err_1(self, filter, filter2):
        """Check combination type mismatch errors."""
        with pytest.raises(ValueError):
            filter.combination_method = FunctCombinationMethod.SUM
            filter2.combination_method = FunctCombinationMethod.MULTIPLY

            filter2.combine(filter)

    def test_combine_err_2(self, filter, filter2):
        """Check combination type mismatch errors."""
        with pytest.raises(ValueError):
            filter.combination_method = FunctCombinationMethod.UNDEFINED
            filter2.combination_method = FunctCombinationMethod.UNDEFINED

            filter2.combine(filter)

    def test_scaling(self, main_funct):
        """Test the scaling."""

        # cases
        # -------
        x = [-1, 0, 2] * u.Hz
        scaling = 2.0

        # truth
        # -------
        truth = [0.5, 0.0, 2.0]

        # test set up
        # -----------------
        scaled_main = main_funct.scale(scaling)

        results = [scaled_main.get_value(x_val) for x_val in x]

        # verification
        # -----------------
        np.testing.assert_allclose(results, truth, rtol=1e-10, atol=1e-17)

    def test_scaling_err(self, main_funct):
        """Check scaling type mismatch errors."""
        with pytest.raises(ValueError):
            main_funct.combination_method = FunctCombinationMethod.SUM

            main_funct.scale(2.0)
