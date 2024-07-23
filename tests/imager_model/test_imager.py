# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import pytest

from opticks import process_paths, u
from opticks.imager_model.imager import Imager


class TestImager:

    file_directory = Path("sat_pushbroom_data")
    alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
    optics_file_path = Path("optics.yaml")
    detector_file_path = Path("ms_detector.yaml")
    rw_electronics_file_path = Path("rw_electronics.yaml")

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
        truth = 0.19090236686390533

        # select the B0 Blue channel
        ref_wavelength = 500 * u.nm

        # Generate the Q Factor
        q_value = imager.q_factor(ref_wavelength, "b0_blue", with_binning=True)

        # verification
        assert q_value == pytest.approx(truth, 1e-9)
