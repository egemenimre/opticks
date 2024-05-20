# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

from opticks import process_paths
from opticks.imager_model.imager import Imager


class TestImager:

    file_directory = Path("sat_pushbroom_data")
    alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
    optics_file_path = Path("optics.yaml")
    detector_file_path = Path("ms_detector.yaml")
    rw_electronics_file_path = Path("rw_electronics.yaml")

    def test_yaml_load(self):
        """Tests the yaml file load."""

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
        Imager.from_yaml_file(
            optics_file_path,
            detector_file_path,
            rw_electronics_file_path,
        )
