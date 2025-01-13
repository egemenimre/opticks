# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import numpy as np
import pytest
from astropy.units import Unit

from opticks.imager_model.detector import Detector
from tests import process_paths


class TestDetector:

    @pytest.fixture(scope="class")
    def detector(self) -> Detector:

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("ms_detector.yaml")

        # different test environments work with different paths
        file_path = process_paths(file_path, file_directory, alt_file_directory)

        return Detector.from_yaml_file(file_path)

    def test_det_read_rate_all_bands(self, detector):
        """Tests the detector read rate."""

        # select all channels
        band_id = [id for id in detector.params.channels.all.keys()]

        # set up
        truth = 73.587 * Unit("megapixel / second")

        # computation
        pix_read_rate = detector.pix_read_rate(
            band_id, with_binning=True, with_tdi=False
        )

        # verification
        np.testing.assert_allclose(
            pix_read_rate, truth, atol=0.00000001 * Unit("megapixel / second")
        )

    def test_det_read_rate_single_band(self, detector):
        """Tests the detector read rate."""

        # select single channel
        band_id = "b0_blue"

        # set up
        truth = 18.39675 * Unit("megapixel / second")

        # computation
        pix_read_rate = detector.pix_read_rate(
            band_id, with_binning=True, with_tdi=False
        )

        # verification
        np.testing.assert_allclose(
            pix_read_rate, truth, atol=0.00000001 * Unit("megapixel / second")
        )
