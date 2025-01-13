# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import math
import numbers

import numpy as np
import pytest
from astropy.units import Quantity, Unit

from opticks import P, u
from opticks.utils.interval_data import FunctCombinationMethod, IntervalData
from opticks.utils.math_utils import InterpolatorWithUnits, InterpolatorWithUnitTypes


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

        y = (0.5 * x.value) ** 2

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
        sample_size = 50

        filter.combination_method = FunctCombinationMethod.MULTIPLY
        filtered_main = main_funct.combine(filter)

        filtered_main.sample_size = sample_size
        resampled_filt_main = filtered_main.resample()

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
        sample_size = 50

        filter_no_units.combination_method = FunctCombinationMethod.MULTIPLY
        filtered_main = main_funct_no_units.combine(filter_no_units)

        filtered_main.sample_size = sample_size
        resampled_filt_main = filtered_main.resample()

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

        np.testing.assert_allclose(sum_functs, 60 * u.m, atol=1e-6 * u.mm)

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

        np.testing.assert_allclose(mult_functs, 800 * u.Hz**2, atol=1e-6 * u.Hz**2)

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

    def test_integration_limited(self):
        """Integration with limited interval."""

        # cases
        # -------
        limited_rng = P.closed(4 * u.s, 6 * u.s)

        # truth
        # -------
        truth = 98.10 * u.m

        # test set up
        # -----------------
        duration = P.closed(0 * u.s, 10 * u.s)

        t = np.linspace(duration.lower, duration.upper, num=100, endpoint=True)

        v = 9.81 * Unit("m/s^2") * t

        ipol = InterpolatorWithUnits.from_ipol_method(
            InterpolatorWithUnitTypes.AKIMA, t, v, extrapolate=True
        )

        free_fall_funct = IntervalData.from_interpolator(ipol)

        results = free_fall_funct.integrate(limited_rng)

        # verification
        # -----------------
        np.testing.assert_allclose(results, truth, atol=1e-5 * u.um)

    def test_integration_complex(self):
        """Integration with complex IntervalData."""

        # test set up
        # -----------------
        data = P.IntervalDict()

        # segment 1
        duration_1 = P.closed(0 * u.s, 10 * u.s)
        vel_1 = 50.0 * u.m / u.s
        data[duration_1] = vel_1

        # segment 2
        duration_2 = P.closed(10 * u.s, 20 * u.s)
        t = np.linspace(duration_2.lower, duration_2.upper, num=100, endpoint=True)
        v = 9.81 * Unit("m/s^2") * t

        ipol = InterpolatorWithUnits.from_ipol_method(
            InterpolatorWithUnitTypes.AKIMA, t, v, extrapolate=True
        )

        data[duration_2] = ipol

        # segment 3
        duration_3 = P.closed(20 * u.s, 30 * u.s)
        vel_3 = 80.0 * u.m / u.s
        data[duration_3] = vel_3

        free_fall_funct = IntervalData(data)

        # truth
        # -------
        truth = (
            vel_1 * (duration_1.upper - duration_1.lower)
            + 0.5 * 9.81 * Unit("m/s^2") * (duration_2.upper**2 - duration_2.lower**2)
            + vel_3 * (duration_3.upper - duration_3.lower)
        )

        # verification
        # -----------------
        results = free_fall_funct.integrate()

        np.testing.assert_allclose(results, truth, atol=1e-5 * u.um)
