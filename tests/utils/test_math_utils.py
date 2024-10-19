# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


import numpy as np
import pytest
from scipy.interpolate import Akima1DInterpolator

from opticks import Q_, u
from opticks.utils.math_utils import InterpolatorWithUnits, InterpolatorWithUnitTypes
from opticks.utils.testing_utils import assert_allclose


class TestInterpolatorWithUnits:

    # define the polynomial with discrete samples
    t = np.linspace(-10, 10, endpoint=True) * u.s
    y = (t.m - 2) * (t.m + 3) * u.m

    @pytest.fixture(scope="class")
    def ipol(self) -> InterpolatorWithUnits:
        """Init the Akima interpolator"""

        # choose the interpolator type
        ipol_type = InterpolatorWithUnitTypes.AKIMA

        # init the interpolator
        ipol = InterpolatorWithUnits.from_ipol_method(
            ipol_type, self.t, self.y, extrapolate=True
        )

        return ipol

    def test_interpolation(self, ipol):

        # target
        tgt = Q_("14.909999999999997 meter")

        # computation
        # ------------

        # test value
        t_tgt = 4.1 * u.s

        # results
        result_1 = ipol(t_tgt)
        result_2 = ipol(t_tgt.to(u.min))

        # comparison
        assert_allclose(result_1, tgt, atol=1e-5 * u.um)
        assert_allclose(result_2, tgt, atol=1e-5 * u.um)

    def test_derivatives(self, ipol):

        # target
        tgt_deriv = Q_("9.2 m/s")
        tgt_antideriv = Q_("230.11199999999997 m·s")
        tgt_integ = Q_("0.6119999999999957 m·s")

        # computation
        # ------------

        # test value
        t_tgt = 4.1 * u.s
        a = -1 * u.s
        b = 4.1 * u.s

        # results
        deriv = ipol.derivative(1)(t_tgt)
        antideriv = ipol.antiderivative(1)(t_tgt)
        integ = ipol.integrate(a, b)

        # comparison
        assert_allclose(deriv, tgt_deriv, atol=1e-5 * u.um / u.s)
        assert_allclose(antideriv, tgt_antideriv, atol=1e-5 * u.um * u.s)
        assert_allclose(integ, tgt_integ, atol=1e-5 * u.um * u.s)

    def test_solve(self, ipol):

        # target
        tgt_solve = [-5.623475382979799, 4.6234753829798] * u.s
        tgt_roots = [-3.0, 2.0] * u.s

        # computation
        # ------------

        # test value
        solve_tgt = 20 * u.m

        # results
        solve = ipol.solve(solve_tgt)
        roots = ipol.roots()

        # comparison
        assert_allclose(solve, tgt_solve, atol=1e-5 * u.us)
        assert_allclose(roots, tgt_roots, atol=1e-5 * u.us)

    def test_interpolation_no_units(self, ipol):

        # target
        tgt = 14.910166202386884

        # computation
        # ------------

        # test value
        t_tgt = 4.1 * u.s

        # Init scipy interpolator, strip units
        scipy_ipol = Akima1DInterpolator(self.t.m, self.y.m, method="makima")

        # Init the Interpolator with units, but try deleting the units of y
        ipol2 = InterpolatorWithUnits(scipy_ipol, self.t, self.y.m)

        # results
        result_1 = ipol2(t_tgt)

        # comparison
        np.testing.assert_allclose(result_1, tgt, atol=1e-5)

    def test_interpolation_args(self, ipol):

        # target
        tgt1 = 15.034782755102038 * u.m
        tgt2 = 14.909999999999997 * u.m

        # computation
        # ------------

        # test value
        t_tgt = 4.1 * u.s

        # Array containing derivatives of the dependent variable.
        # Set with zeros, probably not realistic, but ok for this example
        dydx = np.zeros_like(self.t)

        ipol3 = InterpolatorWithUnits.from_ipol_method(
            InterpolatorWithUnitTypes.CUBICHERMITESPL, self.t, self.y, dydx=dydx
        )

        ipol4 = InterpolatorWithUnits.from_ipol_method(
            InterpolatorWithUnitTypes.CUBICSPL, self.t, self.y, bc_type="not-a-knot"
        )

        # results
        result_1 = ipol3(t_tgt)
        result_2 = ipol4(t_tgt)

        # comparison
        assert_allclose(result_1, tgt1, atol=1e-5 * u.um)
        assert_allclose(result_2, tgt2, atol=1e-5 * u.um)
