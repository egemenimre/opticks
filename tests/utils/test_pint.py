# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import portion as P

from opticks import Q_, u


class TestPint:

    a = P.closed(0 * u.mm, 1 * u.mm)

    b = P.closed(1.2 * u.mm, 2.4 * u.mm)

    def test_basic_units(self):
        """Test basic unit support."""

        truth = "[<Quantity(0, 'millimeter')>,<Quantity(1, 'millimeter')>] | [<Quantity(1.2, 'millimeter')>,<Quantity(2.4, 'millimeter')>]"

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

        truth = "[<Quantity(-inf, 'hertz')>,<Quantity(10, 'hertz')>]"
        truth_oc = "(<Quantity(-inf, 'hertz')>,<Quantity(10, 'hertz')>]"

        # First way to define Inf with units
        ninf = Q_("-inf") * u.Hz
        pinf = Q_("+inf Hz")
        interval_1 = P.closed(ninf, 10 * u.Hz)

        assert truth == str(interval_1)

        # Second way to define Inf with units (may not be fully safe)
        interval_2 = P.closed(-P.inf * u.Hz, 10 * u.Hz)

        assert truth == str(interval_2)

        # full range Interval
        inf = P.open(ninf, pinf)

        # check intersection operation
        assert truth_oc == str(interval_1 & inf)

    def test_intervaldict(self):
        """Test IntervalDict unit support."""

        truth = "{[<Quantity(-10, 'hertz')>,<Quantity(0, 'hertz')>) | (<Quantity(1, 'hertz')>,<Quantity(2, 'hertz')>) | (<Quantity(3, 'hertz')>,<Quantity(10, 'hertz')>]: 0, [<Quantity(0, 'hertz')>,<Quantity(1, 'hertz')>] | [<Quantity(2, 'hertz')>,<Quantity(3, 'hertz')>]: 1.0}"

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
