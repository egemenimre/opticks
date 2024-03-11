# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import numpy as np
import pytest
from pint import DimensionalityError

from opticks import u
from opticks.utils.testing_utils import assert_allclose


class TestTestingUtils:
    """
    Tests for the `testing_utils` package.
    """

    def test_single_item(self):
        calc_val_with_units = 80 * u.cm
        truth_val_with_units = 0.81 * u.m
        atol = 2 * u.cm

        #  use with rtol or atol for errors given in percent or absolute.
        assert_allclose(calc_val_with_units, truth_val_with_units, atol=atol)

    def test_array(self):
        calc_val_with_units = np.asarray([80, 60]) * u.cm
        truth_val_with_units = [0.81, 0.61] * u.m
        atol = 2 * u.cm

        #  use with rtol or atol for errors given in percent or absolute.
        assert_allclose(calc_val_with_units, truth_val_with_units, atol=atol)

    def test_atol_array(self):
        with pytest.raises(ValueError):
            calc_val_with_units = np.asarray([80, 60]) * u.cm
            truth_val_with_units = [0.81, 0.61] * u.m
            atol = [1, 2] * u.cm

            #  use with rtol or atol for errors given in percent or absolute.
            assert_allclose(calc_val_with_units, truth_val_with_units, atol=atol)

    def test_incompatible_units(self):
        with pytest.raises(DimensionalityError):
            calc_val_with_units = 80 * u.cm
            truth_val_with_units = 0.81 * u.m
            atol = 2 * u.Hz

            #  use with rtol or atol for errors given in percent or absolute.
            assert_allclose(calc_val_with_units, truth_val_with_units, atol=atol)
