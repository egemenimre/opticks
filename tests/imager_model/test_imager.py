# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import pytest

from opticks import process_paths, u
from opticks.imager_model.imager import Imager
from opticks.utils.testing_utils import assert_allclose


class TestImager:

    file_directory = Path("sat_pushbroom_data")
    alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
    optics_file_path = Path("optics.yaml")
    detector_file_path = Path("pan_detector.yaml")
    rw_electronics_file_path = Path("rw_electronics.yaml")

    band_id = "pan"

    @pytest.fixture(scope="class")
    def imager(self) -> Imager:
        """Loads the yaml file and inits the Imager."""

        # different test environments work with different paths
        optics_file_path = process_paths(
            self.optics_file_path, self.file_directory, self.alt_file_directory
        )
        detector_file_path = process_paths(
            self.detector_file_path, self.file_directory, self.alt_file_directory
        )
        rw_electronics_file_path = process_paths(
            self.rw_electronics_file_path, self.file_directory, self.alt_file_directory
        )

        # checks only whether the file can be opened and whether it fits the schema
        return Imager.from_yaml_file(
            optics_file_path,
            detector_file_path,
            rw_electronics_file_path,
        )

    def test_q_factor(self, imager: Imager):
        """Tests the Q factor."""

        # check values
        truth = 0.7636094674556213

        # select the PAN channel
        ref_wavelength = 500 * u.nm

        # Generate the Q Factor
        q_value = imager.q_factor(ref_wavelength, self.band_id, with_binning=True)

        # verification
        assert q_value == pytest.approx(truth, 1e-9)

    def test_fov(self, imager: Imager):
        """Tests the horiz and vert FoV."""

        # set up
        horiz_fov_truth = 1.73139509 * u.deg
        vert_fov_truth = 0.0577131695 * u.mdeg

        # computation
        horiz_fov = imager.horizontal_fov(self.band_id)
        vert_fov = imager.vertical_fov(self.band_id)

        # verification
        assert_allclose(horiz_fov, horiz_fov_truth, atol=0.001 * u.mdeg)
        assert_allclose(vert_fov, vert_fov_truth, atol=0.00001 * u.mdeg)

    def test_ifov(self, imager: Imager):
        """Tests the instantaneous FoV."""

        # set up
        truth = 0.0577131695 * u.mdeg

        # computation
        ifov = imager.ifov(self.band_id, with_binning=False)

        # verification
        assert_allclose(ifov, truth, atol=0.00000001 * u.mdeg)
