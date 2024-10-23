# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import numpy as np
import portion as P
import pytest

from opticks import u
from opticks.utils.interval_data import IntervalData
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

        # cases
        # -------
        x = [1.9999, 2.0, 2.0000001] * u.Hz

        # truth
        # -------
        truth = [0.0, 1.0, 1.00000010]

        # test set up
        # -----------------
        filtered_main = main_funct.combine(filter)
        results = [filtered_main.get_value(x_val) for x_val in x]

        # verification
        # -----------------
        np.testing.assert_allclose(results, truth, rtol=1e-10)

    def test_resampled_edge(self, filter, main_funct):

        # cases
        # -------
        x = [1.9999, 2.0, 2.0000001] * u.Hz

        # truth
        # -------
        truth = [0.0, 1.0, 1.00000010]

        # test set up
        # -----------------
        approx_stepsize = 0.15 * u.Hz

        filtered_main = main_funct.combine(filter)
        resampled_filt_main = filtered_main.resample(approx_stepsize)

        results = [resampled_filt_main.get_value(x_val) for x_val in x]

        # verification
        # -----------------
        np.testing.assert_allclose(results, truth, rtol=1e-8)

    def test_resampled_edge_no_units(self, filter_no_units, main_funct_no_units):

        # cases
        # -------
        x = [1.9999, 2.0, 2.0000001]

        # truth
        # -------
        truth = [0.0, 1.0, 1.00000010]

        # test set up
        # -----------------
        approx_stepsize = 0.15

        filtered_main = main_funct_no_units.combine(filter_no_units)
        resampled_filt_main = filtered_main.resample(approx_stepsize)

        results = [resampled_filt_main.get_value(x_val) for x_val in x]

        # verification
        # -----------------
        np.testing.assert_allclose(results, truth, rtol=1e-8)
