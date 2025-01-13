# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


import numpy as np
import pytest
from astropy.units import Quantity, Unit

from opticks import Q_, P, u
from opticks.utils.unit_utils import (
    quantity_from_list,
    split_value_and_force_unit,
    split_value_and_unit,
)


class TestUnitUtils:

    def test_quantity_from_list(self):
        """Test Quantity from list."""
        data = [Q_(1, "mm"), Q_(2, "m"), Q_(3, "m")]

        target = [1e-3, 2, 3] * u.m

        np.testing.assert_array_equal(quantity_from_list(data), target)

    def test_quantity_from_list_with_error(self):
        with pytest.raises(TypeError):
            """Test Quantity from list with error."""

            data = [1.0, Q_(2, "m"), Q_(3, "m")]
            quantity_from_list(data)

    def test_split_quantity_input_with_forced(self):
        """Test Quantity input with target units."""
        q = Quantity(5, "m")
        tgt_unit = u.km
        value, unit = split_value_and_force_unit(q, tgt_unit)

        assert value == 5e-3
        assert unit == tgt_unit

    def test_split_quantity_input(self):
        """Test Quantity input."""
        q = Quantity(5, "m")
        value, unit = split_value_and_unit(q)

        assert value == 5
        assert unit == q.unit

    def test_split_float_input(self):
        """Test float input."""
        value, unit = split_value_and_unit(5.0)

        assert value == 5.0
        assert unit is None

    def test_split_array_input(self):
        """Test array input."""
        arr = np.array([1, 2, 3]) * u.m
        value, unit = split_value_and_unit(arr)

        np.testing.assert_array_equal(value, np.array([1, 2, 3]))
        assert unit == arr.unit

    def test_split_array_without_units(self):
        """Test array without units."""
        arr = np.array([1, 2, 3])
        value, unit = split_value_and_unit(arr)

        np.testing.assert_array_equal(value, arr)
        assert unit is None

    def test_split_dimensionless_quantity(self):
        """Test dimensionless Quantity."""
        q = Quantity(5, None)
        value, unit = split_value_and_unit(q)

        assert value == 5
        assert unit is Unit()


class TestUnits:

    a = P.closed(0 * u.mm, 1 * u.mm)

    b = P.closed(1.2 * u.mm, 2.4 * u.mm)

    def test_basic_units(self):
        """Test basic unit support."""

        truth = "[<Quantity 0. mm>,<Quantity 1. mm>] | [<Quantity 1.2 mm>,<Quantity 2.4 mm>]"

        c = self.a | self.b

        assert truth == str(c)

    def test_inf(self):
        """Test P.inf."""

        truth = [True, True, True]

        float_interval = P.closed(0, 1)

        result = [P.inf > float_interval, P.inf > 10 * u.mm, -P.inf < self.b]

        assert truth == result

    def test_inf_def(self):
        """Test P.inf definition."""

        truth = "[<Quantity -inf Hz>,<Quantity 10. Hz>]"
        truth_oc = "(<Quantity -inf Hz>,<Quantity 10. Hz>]"

        # First way to define Inf with units
        ninf = Q_("-inf") * u.Hz
        pinf = Q_("+inf Hz")
        interval_1 = P.closed(ninf, 10 * u.Hz)

        assert truth == str(interval_1)

        # Second way to define Inf with units (may not be fully safe)
        interval_2 = P.closed(-P.inf, 10 * u.Hz)

        assert truth == str(interval_2)

        # full range Interval
        inf = P.open(ninf, pinf)

        # check intersection operation
        assert truth_oc == str(interval_1 & inf)

    def test_intervaldict(self):
        """Test IntervalDict unit support."""

        truth = "{[<Quantity -10. Hz>,<Quantity 0. Hz>) | (<Quantity 1. Hz>,<Quantity 2. Hz>) | (<Quantity 3. Hz>,<Quantity 10. Hz>]: 0, [<Quantity 0. Hz>,<Quantity 1. Hz>] | [<Quantity 2. Hz>,<Quantity 3. Hz>]: 1.0}"

        data = P.IntervalDict()
        # interval of validity
        validity = P.closed(-10 * u.Hz, 10 * u.Hz)
        data[validity] = 0
        # data proper
        first_rng = P.closed(0 * u.Hz, 1 * u.Hz)
        data[first_rng] = 1.0
        second_rng = P.closed(2 * u.Hz, 3 * u.Hz)
        data[second_rng] = 1.0

        assert truth == str(data)
