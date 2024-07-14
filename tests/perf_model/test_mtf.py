# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import pytest

from opticks import process_paths, u
from opticks.imager_model.detector import Channel, Detector
from opticks.imager_model.optics import Optics
from opticks.perf_model.mtf import MTF_Model


class TestMTF:

    file_directory = Path("sat_pushbroom_data")
    alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")

    ref_wavelength = 650 * u.nm
    input_line_freq = 30 * u.cy / u.mm

    @pytest.fixture(scope="class")
    def optics(self) -> Optics:

        file_path = Path("optics.yaml")

        # different test environments work with different paths
        file_path = process_paths(
            file_path, self.file_directory, self.alt_file_directory
        )

        return Optics.from_yaml_file(file_path)

    @pytest.fixture(scope="class")
    def detector(self) -> Detector:

        file_path = Path("pan_detector.yaml")

        # different test environments work with different paths
        file_path = process_paths(
            file_path, self.file_directory, self.alt_file_directory
        )

        return Detector.from_yaml_file(file_path)

    def test_mtf_ideal_optics(self, optics):
        """Tests the ideal optics MTF."""

        # check values
        truth = 0.5196720931409163

        # Generate the MTF model and values
        mtf_model = MTF_Model.ideal_optics(self.ref_wavelength, optics)

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_aberrated_optics(self, optics):
        """Tests the aberrated optics MTF."""

        # check values
        truth = 0.4816165574868274

        # Generate the MTF model and values
        w_rms = 0.05
        mtf_model = MTF_Model.aberrated_optics(self.ref_wavelength, w_rms, optics)

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_detector_sampling(self, detector: Detector):
        """Tests the detector sampling MTF."""

        # check values
        truth = 0.7679273089188128

        # select the first channel
        channel: Channel = next(iter(detector.params.channels.all.values()))

        # Generate the MTF model and values
        mtf_model = MTF_Model.detector_sampling(channel.pixel_pitch())

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)
