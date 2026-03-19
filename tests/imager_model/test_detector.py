# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import numpy as np
import pytest
from astropy.units import Unit, isclose

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
        band_id = [id for id in detector.channels.keys()]

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

    def test_yaml_string_round_trip(self, detector):
        """Tests that a Detector survives a YAML string write/read round trip."""

        # write to YAML string, then read back
        yaml_text = detector.to_yaml_text()
        detector_reloaded = Detector.from_yaml_text(yaml_text)

        # verify scalar fields
        assert detector_reloaded.name == detector.name
        assert detector_reloaded.detector_type == detector.detector_type
        assert detector_reloaded.horizontal_pixels == detector.horizontal_pixels
        assert detector_reloaded.vertical_pixels == detector.vertical_pixels

        # verify Quantity fields
        assert isclose(detector_reloaded.pixel_pitch, detector.pixel_pitch)
        assert isclose(
            detector_reloaded.timings.frame_rate, detector.timings.frame_rate
        )
        assert isclose(
            detector_reloaded.timings.integration_duration,
            detector.timings.integration_duration,
        )

        # verify channels survive the round trip
        assert set(detector_reloaded.channels.keys()) == set(detector.channels.keys())
        for band_id, channel in detector.channels.items():
            ch_reloaded = detector_reloaded.channels[band_id]
            assert ch_reloaded.name == channel.name
            assert ch_reloaded.horizontal_pixels == channel.horizontal_pixels
            assert ch_reloaded.vertical_pixels == channel.vertical_pixels
            assert ch_reloaded.binning == channel.binning
            assert ch_reloaded.tdi_stages == channel.tdi_stages
            assert isclose(ch_reloaded.cuton_wvl, channel.cuton_wvl)
            assert isclose(ch_reloaded.cutoff_wvl, channel.cutoff_wvl)

    def test_yaml_file_round_trip(self, detector, tmp_path):
        """Tests that a Detector survives a YAML file write/read round trip."""

        # write to file, then read back
        out_file = tmp_path / "detector_round_trip.yaml"
        detector.to_yaml_file(out_file)
        detector_reloaded = Detector.from_yaml_file(out_file)

        # verify scalar fields
        assert detector_reloaded.name == detector.name
        assert detector_reloaded.detector_type == detector.detector_type
        assert detector_reloaded.horizontal_pixels == detector.horizontal_pixels
        assert detector_reloaded.vertical_pixels == detector.vertical_pixels

        # verify Quantity fields
        assert isclose(detector_reloaded.pixel_pitch, detector.pixel_pitch)
        assert isclose(
            detector_reloaded.timings.frame_rate, detector.timings.frame_rate
        )
        assert isclose(
            detector_reloaded.timings.integration_duration,
            detector.timings.integration_duration,
        )

        # verify channels survive the round trip
        assert set(detector_reloaded.channels.keys()) == set(detector.channels.keys())
        for band_id, channel in detector.channels.items():
            ch_reloaded = detector_reloaded.channels[band_id]
            assert ch_reloaded.name == channel.name
            assert ch_reloaded.horizontal_pixels == channel.horizontal_pixels
            assert ch_reloaded.vertical_pixels == channel.vertical_pixels
            assert ch_reloaded.binning == channel.binning
            assert ch_reloaded.tdi_stages == channel.tdi_stages
            assert isclose(ch_reloaded.cuton_wvl, channel.cuton_wvl)
            assert isclose(ch_reloaded.cutoff_wvl, channel.cutoff_wvl)
