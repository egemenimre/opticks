# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import pytest

from opticks import process_paths, u
from opticks.imager_model.detector import Detector
from opticks.imager_model.optics import Optics
from opticks.utils.testing_utils import assert_allclose


class TestDetector:

    @pytest.fixture(scope="class")
    def detector_input(self) -> Detector:

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("pan_detector.yaml")

        # different test environments work with different paths
        file_path = process_paths(file_path, file_directory, alt_file_directory)

        return Detector.from_yaml_file(file_path)

    @pytest.fixture(scope="class")
    def optics_input(self) -> Optics:

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("optics.yaml")

        # different test environments work with different paths
        file_path = process_paths(file_path, file_directory, alt_file_directory)

        return Optics.from_yaml_file(file_path)

    def test_ifov(self, detector_input, optics_input):
        """Tests the instantaneous FoV."""
        detector = detector_input
        optics = optics_input
        # select the first channel
        channel = next(iter(detector.params.channels.all.values()))

        # set up
        truth = 0.05771756169469254 * u.mdeg

        # computation
        ifov = channel.ifov(optics, with_binning=False)

        # verification
        assert_allclose(ifov, truth, atol=0.00000001 * u.mdeg)

    def test_fov(self, detector_input, optics_input):
        """Tests the horiz and vert FoV."""
        detector = detector_input
        optics = optics_input
        # select the first channel
        channel = next(iter(detector.params.channels.all.values()))

        # set up
        horiz_fov_truth = 1.7316586464211057 * u.deg
        vert_fov_truth = 5.771756169469742e-05 * u.deg

        # computation
        horiz_fov = channel.horizontal_fov(optics)
        vert_fov = channel.vertical_fov(optics)

        # verification
        assert_allclose(horiz_fov, horiz_fov_truth, atol=0.001 * u.mdeg)
        assert_allclose(vert_fov, vert_fov_truth, atol=0.0000001 * u.mdeg)

    def test_det_read_rate(self, detector_input):
        """Tests the detector read rate."""
        detector = detector_input
        # select all channels
        channel = detector.params.channels.all.values()

        # set up
        truth = 294.348 * u["megapixel / second"]

        # computation
        pix_read_rate = detector.pix_read_rate(
            channel, with_binning=True, with_tdi=False
        )

        # verification
        assert_allclose(pix_read_rate, truth, atol=0.00000001 * u["megapixel / second"])
