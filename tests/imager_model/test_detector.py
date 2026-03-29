# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import numpy as np
import pytest
from astropy.units import Unit, isclose

from opticks import u
from opticks.imager_model.detector import Detector, SensorParams
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

    @pytest.fixture(scope="class")
    def block_detector(self) -> Detector:
        """Detector with block read functionality."""

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("block_test_pan_detector.yaml")

        # different test environments work with different paths
        file_path = process_paths(file_path, file_directory, alt_file_directory)

        return Detector.from_yaml_file(file_path)

    def test_det_read_rate_all_bands(self, detector):
        """Tests the detector read rate."""

        # select all channels
        band_id = list(detector.channels.keys())

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
            assert ch_reloaded.read_blocks == channel.read_blocks
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
            assert ch_reloaded.read_blocks == channel.read_blocks
            assert isclose(ch_reloaded.cuton_wvl, channel.cuton_wvl)
            assert isclose(ch_reloaded.cutoff_wvl, channel.cutoff_wvl)

    def test_block_read_rate_with_read_blocks(self, block_detector):
        """Tests that read_blocks multiplies the pixel read rate."""

        band_id = "pan"

        # read_blocks=4, frame_rate=2452.9 Hz
        # 30000 * 4 * 2452.9 = 294.348 Mpix/s (same as original pan with read_blocks=1)
        truth = 294.348 * Unit("megapixel / second")

        pix_read_rate = block_detector.pix_read_rate(
            band_id, with_binning=True, with_tdi=False, with_read_blocks=True
        )

        np.testing.assert_allclose(
            pix_read_rate, truth, atol=0.00000001 * Unit("megapixel / second")
        )

    def test_block_read_rate_without_read_blocks(self, block_detector):
        """Tests the pixel read rate with read_blocks disabled."""

        band_id = "pan"

        # with_read_blocks=False: 30000 * 1 * 2452.9 = 73.587 Mpix/s
        truth = 73.587 * Unit("megapixel / second")

        pix_read_rate = block_detector.pix_read_rate(
            band_id, with_binning=True, with_tdi=False, with_read_blocks=False
        )

        np.testing.assert_allclose(
            pix_read_rate, truth, atol=0.00000001 * Unit("megapixel / second")
        )


class TestSensorParams:
    @pytest.fixture(scope="class")
    def detector(self) -> Detector:
        """PAN detector with sensor_params configured."""

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("pan_detector_with_sensor_params.yaml")

        file_path = process_paths(file_path, file_directory, alt_file_directory)

        det = Detector.from_yaml_file(file_path)
        assert det.sensor_params is not None
        det.sensor_params.load_absorption_data()
        return det

    def test_sensor_params_loaded(self, detector):
        """sensor_params block is parsed from YAML."""

        sp = detector.sensor_params
        assert sp is not None
        assert sp.diffusion_model == "bsi_1"
        assert isclose(sp.diffusion_length, 75 * u.um)
        assert isclose(sp.field_free_depth, 7 * u.um)
        assert isclose(sp.depletion_depth, 7 * u.um)

    def test_absorption_coeff_at_known_wavelength(self, detector):
        """get_absorption_coeff returns the tabulated value at an exact data point."""

        # 650 nm → 0.27 1/um exactly as tabulated
        alpha = detector.sensor_params.get_absorption_coeff(650 * u.nm)
        np.testing.assert_allclose(alpha.to(1 / u.um).value, 0.27, rtol=1e-4)

    def test_absorption_coeff_interpolated(self, detector):
        """get_absorption_coeff interpolates between tabulated points."""

        # 625 nm is between 600 (0.40) and 650 (0.27) — should be between them
        alpha = detector.sensor_params.get_absorption_coeff(625 * u.nm)
        alpha_val = alpha.to(1 / u.um).value
        assert 0.27 < alpha_val < 0.40

    def test_get_diffusion_mtf_1d_returns_model(self, detector):
        """get_diffusion_mtf_1d returns a usable MTF_Model_1D at a given wavelength."""

        from opticks.contrast_model.mtf import MTF_Model_1D

        mtf_model = detector.get_diffusion_mtf_1d(650 * u.nm)
        assert isinstance(mtf_model, MTF_Model_1D)

    def test_get_diffusion_mtf_1d_values_in_range(self, detector):
        """Diffusion MTF values are in [0, 1] and decrease with frequency."""

        mtf_model = detector.get_diffusion_mtf_1d(650 * u.nm)
        freq = np.linspace(0, detector.channels["pan"].nyquist_freq().value, 50)
        freq_qty = freq * u.cy / u.mm
        mtf_vals = mtf_model.mtf_value(freq_qty)

        assert np.all(mtf_vals >= 0)
        assert np.all(mtf_vals <= 1)
        # MTF should be non-increasing
        assert np.all(np.diff(mtf_vals) <= 1e-9)

    def test_yaml_round_trip_with_sensor_params(self, detector):
        """sensor_params survives a YAML string round trip."""

        yaml_text = detector.to_yaml_text()
        reloaded = Detector.from_yaml_text(yaml_text)

        orig_sp = detector.sensor_params
        assert orig_sp is not None
        sp = reloaded.sensor_params
        assert sp is not None
        assert sp.diffusion_model == orig_sp.diffusion_model
        assert isclose(sp.diffusion_length, orig_sp.diffusion_length)
        assert orig_sp.absorption_data is not None
        assert sp.absorption_data is not None
        assert sp.absorption_data.file == orig_sp.absorption_data.file

    def test_missing_sensor_params_raises(self):
        """get_diffusion_mtf_1d raises ValueError when sensor_params is absent."""

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("pan_detector.yaml")
        file_path = process_paths(file_path, file_directory, alt_file_directory)

        det = Detector.from_yaml_file(file_path)
        with pytest.raises(ValueError, match="sensor_params"):
            det.get_diffusion_mtf_1d(650 * u.nm)

    def test_missing_absorption_data_raises(self):
        """get_absorption_coeff raises ValueError when absorption data not loaded."""

        sp = SensorParams(
            diffusion_model="bsi_1",
            diffusion_length=75 * u.um,
            field_free_depth=7 * u.um,
            depletion_depth=7 * u.um,
        )
        with pytest.raises(ValueError, match="not loaded"):
            sp.get_absorption_coeff(650 * u.nm)

    def test_get_diffusion_mtf_1d_with_preset(self):
        """get_diffusion_mtf_1d works when sensor_params uses a diffusion_preset."""
        from opticks.contrast_model.mtf import MTF_Model_1D

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("pan_detector_with_preset.yaml")
        file_path = process_paths(file_path, file_directory, alt_file_directory)

        det = Detector.from_yaml_file(file_path)
        assert det.sensor_params is not None
        det.sensor_params.load_absorption_data()

        mtf_model = det.get_diffusion_mtf_1d(650 * u.nm)
        assert isinstance(mtf_model, MTF_Model_1D)

    def test_invalid_diffusion_model_label_raises(self):
        """get_diffusion_mtf_1d raises ValueError for an unknown diffusion_model label."""

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("pan_detector_with_sensor_params.yaml")
        file_path = process_paths(file_path, file_directory, alt_file_directory)

        det = Detector.from_yaml_file(file_path)
        assert det.sensor_params is not None
        det.sensor_params.load_absorption_data()
        det.sensor_params.diffusion_model = "invalid_model"

        with pytest.raises(ValueError, match="Unknown diffusion_model"):
            det.get_diffusion_mtf_1d(650 * u.nm)

    def test_mutually_exclusive_model_and_preset_raises(self):
        """SensorParams raises ValueError when both diffusion_model and diffusion_preset are set."""
        with pytest.raises(ValueError, match="mutually exclusive"):
            SensorParams(
                diffusion_model="bsi_1",
                diffusion_preset="scientific_ccd",
            )


class TestDetectorSamplingMTF:
    @pytest.fixture(scope="class")
    def ms_detector(self) -> Detector:
        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = process_paths(
            Path("ms_detector.yaml"), file_directory, alt_file_directory
        )
        return Detector.from_yaml_file(file_path)

    @pytest.fixture(scope="class")
    def binned_detector(self) -> Detector:
        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = process_paths(
            Path("binned_test_detector.yaml"), file_directory, alt_file_directory
        )
        return Detector.from_yaml_file(file_path)

    def test_get_det_sampling_mtf_1d_no_channel(self, ms_detector):
        """No channel supplied → uses native detector pixel pitch."""
        from opticks.contrast_model.mtf import MTF_Model_1D

        mtf_model = ms_detector.get_det_sampling_mtf_1d()
        assert isinstance(mtf_model, MTF_Model_1D)

    def test_get_det_sampling_mtf_1d_with_channel(self, ms_detector):
        """Channel supplied → returns MTF_Model_1D using channel effective pitch."""
        from opticks.contrast_model.mtf import MTF_Model_1D

        mtf_model = ms_detector.get_det_sampling_mtf_1d("b0_blue")
        assert isinstance(mtf_model, MTF_Model_1D)

    def test_get_det_sampling_mtf_1d_values_in_range(self, ms_detector):
        """Sampling MTF values are in [0, 1] up to Nyquist."""
        mtf_model = ms_detector.get_det_sampling_mtf_1d("b0_blue")
        nyquist = ms_detector.get_channel("b0_blue").nyquist_freq().value
        freq = np.linspace(0, nyquist, 50) * u.cy / u.mm
        mtf_vals = mtf_model.mtf_value(freq)

        assert np.all(mtf_vals >= 0)
        assert np.all(mtf_vals <= 1)

    def test_get_det_sampling_mtf_1d_binning_changes_pitch(self, binned_detector):
        """Binned channel uses effective pitch (native × binning), giving lower MTF at same freq."""
        freq = 20 * u.cy / u.mm

        mtf_native = binned_detector.get_det_sampling_mtf_1d("native").mtf_value(freq)
        mtf_binned = binned_detector.get_det_sampling_mtf_1d("binned_4x").mtf_value(
            freq
        )
        mtf_no_channel = binned_detector.get_det_sampling_mtf_1d().mtf_value(freq)

        # binned channel (40 µm pitch) drops faster than native (10 µm pitch)
        assert mtf_binned < mtf_native
        # no channel uses native pitch
        np.testing.assert_allclose(mtf_no_channel, mtf_native, rtol=1e-9)
