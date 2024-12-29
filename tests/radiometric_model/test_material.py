# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


import numpy as np
import portion as P
import pytest

from opticks import u
from opticks.radiometric_model.material import OpticalMaterial
from opticks.utils.interval_data import IntervalData
from opticks.utils.math_utils import InterpolatorWithUnits, InterpolatorWithUnitTypes


class TestOpticalMaterial:

    @pytest.fixture(scope="class")
    def reflectivity_correct_input(self) -> IntervalData:
        """Correct reflectivity with an interpolator in between."""

        sub_range = P.closed(1000 * u.nm, 1800 * u.nm)
        x = np.linspace(sub_range.lower, sub_range.upper, num=100, endpoint=True)
        y = x / (2500 * u.nm)

        return y

    @pytest.fixture(scope="class")
    def reflectivity_err_type_1_input(self) -> IntervalData:
        """Incorrect reflectivity as the edges of the
        interpolator in between exceeds the limits."""

        sub_range = P.closed(1000 * u.nm, 1800 * u.nm)
        x = np.linspace(sub_range.lower, sub_range.upper, num=100, endpoint=True)
        y = (x - 1200 * u.nm) * (x - 1700 * u.nm) / (400 * u.nm) ** 2 + 0.4

        return y

    @pytest.fixture(scope="class")
    def reflectivity_err_type_2_input(self) -> IntervalData:
        """Incorrect reflectivity as the minima/maxima of the
        interpolator in between exceeds the limits."""

        sub_range = P.closed(1000 * u.nm, 1800 * u.nm)
        x = np.linspace(sub_range.lower, sub_range.upper, num=100, endpoint=True)
        y = (x - 1200 * u.nm) * (x - 1500 * u.nm) * (x - 1900 * u.nm) / (
            500 * u.nm
        ) ** 3 + 1.0

        return y

    @pytest.fixture(scope="class")
    def reflectivity_err_type_3(self) -> IntervalData:
        """Incorrect reflectivity as the summation is not equal to 1."""

        # models a fixed value at 0.333
        data = P.IntervalDict()
        range = P.closed(300 * u.nm, 2500 * u.nm)
        x = np.linspace(range.lower, range.upper, num=100, endpoint=True)
        y = 0.333 + 0 * x / u.nm

        ipol = InterpolatorWithUnits.from_ipol_method(
            InterpolatorWithUnitTypes.AKIMA, x, y, extrapolate=True
        )

        data[range] = ipol

        reflectivity = IntervalData(data)

        return reflectivity

    @pytest.fixture(scope="class")
    def reflectivity_err_type_4(self) -> IntervalData:
        """Incorrect reflectivity value."""

        # models a fixed value at 0.333
        data = P.IntervalDict()
        # interval of validity
        data[P.closed(300 * u.nm, 2500 * u.nm)] = 1.2

        reflectivity = IntervalData(data)

        return reflectivity

    @pytest.fixture(scope="class")
    def reflectivity_model(self):
        """Correct reflectivity with an interpolator in between."""

        def _reflectivity_model(y):
            data = P.IntervalDict()
            # interval of validity
            data[P.closed(300 * u.nm, 2500 * u.nm)] = 0
            # data proper
            data[P.closed(300 * u.nm, 1000 * u.nm)] = 0.8

            # add an interpolator to make it more complex
            sub_range = P.closed(1000 * u.nm, 1800 * u.nm)
            x = np.linspace(sub_range.lower, sub_range.upper, num=100, endpoint=True)

            ipol = InterpolatorWithUnits.from_ipol_method(
                InterpolatorWithUnitTypes.AKIMA, x, y, extrapolate=True
            )

            data[sub_range] = ipol

            data[P.closed(1800 * u.nm, 2500 * u.nm)] = 0.3

            reflectivity = IntervalData(data)

            return reflectivity

        return _reflectivity_model

    @pytest.mark.parametrize(
        "input, expected",
        [
            (500 * u.nm, 0.2),
            (1000 * u.nm, 0.6),
            (1050 * u.nm, 0.58),
            (1500 * u.nm, 0.4),
            (1799.99 * u.nm, 0.280004),
            (1800 * u.nm, 0.7),
            (2100 * u.nm, 0.7),
        ],
    )
    def test_init_err(
        self, reflectivity_model, reflectivity_correct_input, input, expected
    ):
        """Check proper init."""

        reflectivity = reflectivity_model(reflectivity_correct_input)

        opt_mat = OpticalMaterial.init_opaque_from_refl(reflectivity)

        np.testing.assert_allclose(
            opt_mat.absorptivity.get_value(input), expected, rtol=1e-13
        )

    def test_init_err_1(self, reflectivity_model, reflectivity_err_type_1_input):
        """Check edges exceed 1 type errors."""
        with pytest.raises(ValueError):

            reflectivity = reflectivity_model(reflectivity_err_type_1_input)

            OpticalMaterial.init_opaque_from_refl(reflectivity)

    def test_init_err_2(self, reflectivity_model, reflectivity_err_type_2_input):
        """Check maxima/minima exceed 1 type errors."""
        with pytest.raises(ValueError):

            reflectivity = reflectivity_model(reflectivity_err_type_2_input)

            OpticalMaterial.init_opaque_from_refl(reflectivity)

    def test_init_err_3(self, reflectivity_err_type_3):
        """Check summation not equal to 1 errors."""
        with pytest.raises(ValueError):

            OpticalMaterial(
                reflectivity_err_type_3,
                reflectivity_err_type_3,
                reflectivity_err_type_3,
            )

    def test_init_err_4(self, reflectivity_err_type_4):
        """Check single value exceeds 1 type errors."""
        with pytest.raises(ValueError):

            OpticalMaterial.init_opaque_from_refl(reflectivity_err_type_4)
