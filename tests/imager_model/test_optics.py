# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import pytest

from opticks import process_paths, u
from opticks.imager_model.optics import Optics
from opticks.utils.testing_utils import assert_allclose


class TestOptics:

    @pytest.fixture(scope="class")
    def optics(self) -> Optics:

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("optics.yaml")

        # different test environments work with different paths
        file_path = process_paths(file_path, file_directory, alt_file_directory)

        return Optics.from_yaml_file(file_path)

    def test_fov(self, optics):
        """Tests the full optical FoV."""

        truth = 1.7757828128191897 * u.deg
        fov = optics.full_optical_fov

        assert_allclose(fov, truth, atol=0.001 * u.mdeg)

    def test_spatial_cutoff(self, optics):
        """Tests the spatial cutoff frequency."""

        # set up
        ref_wavelength = 640 * u.nm

        truth = 78.70011623401781 * u.cy / u.mm

        # computation
        cutoff_freq = optics.spatial_cutoff_freq(ref_wavelength)

        # verification
        assert_allclose(cutoff_freq, truth, atol=0.00001 * u.cy / u.mm)
