# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""Tests for the ImagingChain top-level container."""

from pathlib import Path

import pytest
from numpy.testing import assert_allclose

from opticks import u
from opticks.imaging_model.imager import Imager
from opticks.imaging_model.imaging_chain import ImagingChain
from opticks.imaging_model.processing import Processing
from tests import process_paths


class TestImagingChain:
    file_directory = Path("sat_pushbroom_data")
    alt_file_directory = Path("tests", "imaging_model", "sat_pushbroom_data")

    @pytest.fixture(scope="class")
    def imager(self) -> Imager:
        optics_file_path = process_paths(
            Path("optics.yaml"), self.file_directory, self.alt_file_directory
        )
        detector_file_path = process_paths(
            Path("pan_detector.yaml"), self.file_directory, self.alt_file_directory
        )
        rw_electronics_file_path = process_paths(
            Path("rw_electronics.yaml"),
            self.file_directory,
            self.alt_file_directory,
        )
        return Imager.from_yaml_file(
            optics_file_path, detector_file_path, rw_electronics_file_path
        )

    @pytest.fixture(scope="class")
    def processing(self) -> Processing:
        processing_file_path = process_paths(
            Path("processing.yaml"), self.file_directory, self.alt_file_directory
        )
        return Processing.from_yaml_file(processing_file_path)

    # --- Construction ---

    def test_holds_imager_and_processing(self, imager: Imager, processing: Processing):
        """ImagingChain stores the exact objects passed in."""
        chain = ImagingChain(imager, processing)
        assert chain.imager is imager
        assert chain.processing is processing

    def test_processing_optional(self, imager: Imager):
        """Processing defaults to None when omitted."""
        chain = ImagingChain(imager)
        assert chain.imager is imager
        assert chain.processing is None

    # --- End-to-end smoke ---

    def test_chain_components_usable(self, imager: Imager, processing: Processing):
        """Imager and Processing remain functional through the chain."""
        chain = ImagingChain(imager, processing)
        # Imager hardware property still reachable
        assert chain.imager.optics.focal_length.to(u.mm).value > 0
        # Processing produces a working MTF model with MTF(0) = 1
        assert chain.processing is not None
        mtf_model = chain.processing.get_resampling_mtf_1d(input_pitch=0.5 * u.m)
        assert mtf_model.mtf_value(0.0 * u.cy / u.m) == pytest.approx(1.0, abs=1e-9)

    # --- Projection ---

    def test_ssd(self, imager: Imager):
        """Tests the Spatial Sample Distance projections."""
        chain = ImagingChain(imager)
        band_id = "pan"
        sat_altitude = 694.0 * u.km

        # centre
        ssd_horiz, ssd_vert = chain.spatial_sample_distance(
            sat_altitude, band_id, False, "centre"
        )
        assert_allclose(ssd_horiz, 0.699055672 * u.m, atol=0.00000001 * u.m)
        assert_allclose(ssd_vert, 0.699055672 * u.m, atol=0.00000001 * u.m)

        # centre right
        ssd_horiz, ssd_vert = chain.spatial_sample_distance(
            sat_altitude, band_id, False, "centre right"
        )
        assert_allclose(ssd_horiz, 0.699215273 * u.m, atol=0.00000001 * u.m)
        assert_allclose(ssd_vert, 0.699135473 * u.m, atol=0.00000001 * u.m)
