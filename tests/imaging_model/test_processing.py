# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""Tests for the Processing imager component."""

from pathlib import Path

import numpy as np
import pytest

from opticks import u
from opticks.imaging_model.processing import Processing, ProcessingParams
from tests import process_paths


class TestProcessing:
    file_directory = Path("sat_pushbroom_data")
    alt_file_directory = Path("tests", "imaging_model", "sat_pushbroom_data")

    @pytest.fixture(scope="class")
    def processing(self) -> Processing:
        file_path = process_paths(
            Path("processing.yaml"),
            self.file_directory,
            self.alt_file_directory,
        )
        return Processing.from_yaml_file(file_path)

    # --- YAML round-trip ---

    def test_yaml_round_trip(self, processing: Processing):
        """from_yaml_file → to_yaml_text → from_yaml_text reproduces the object."""
        yaml_text = processing.to_yaml_text()
        restored = Processing.from_yaml_text(yaml_text)
        assert restored.name == processing.name
        assert restored.processing_params is not None
        assert processing.processing_params is not None
        assert restored.processing_params.resampling_kernel == "bilinear"
        assert (
            restored.processing_params.output_pitch
            == processing.processing_params.output_pitch
        )

    def test_loads_kernel_and_pitch(self, processing: Processing):
        """YAML fixture loads kernel and output_pitch correctly."""
        assert processing.processing_params is not None
        assert processing.processing_params.output_pitch is not None
        assert processing.processing_params.resampling_kernel == "bilinear"
        assert processing.processing_params.output_pitch.to(u.m).value == pytest.approx(
            1.0
        )

    # --- get_resampling_mtf_1d ---

    def test_get_resampling_mtf_1d_dc_equal_one(self, processing: Processing):
        """MTF model returned by get_resampling_mtf_1d has MTF(0) = 1."""
        model = processing.get_resampling_mtf_1d(input_pitch=1.0 * u.m)
        assert model.mtf_value(0.0 * u.cy / u.m) == pytest.approx(1.0, abs=1e-9)

    def test_get_resampling_mtf_1d_output_pitch_override(self, processing: Processing):
        """Runtime output_pitch override changes the MTF relative to the YAML default."""
        freqs = np.array([0.1, 0.3]) * u.cy / u.m
        model_default = processing.get_resampling_mtf_1d(input_pitch=0.5 * u.m)
        model_override = processing.get_resampling_mtf_1d(
            input_pitch=0.5 * u.m, output_pitch=2.0 * u.m
        )
        # p_eff: default uses output=1m (>input=0.5m), override uses output=2m
        assert not np.allclose(
            model_default.mtf_value(freqs), model_override.mtf_value(freqs)
        )

    def test_get_resampling_mtf_1d_varying_input_pitch(self, processing: Processing):
        """Different input_pitch values produce different MTF curves (SSD variation)."""
        freq = 0.3 * u.cy / u.m
        mtf_upsample = processing.get_resampling_mtf_1d(
            input_pitch=2.0 * u.m
        ).mtf_value(freq)
        mtf_downsample = processing.get_resampling_mtf_1d(
            input_pitch=0.5 * u.m
        ).mtf_value(freq)
        # upsample: p_eff=2m → nu=0.6; downsample: p_eff=1m → nu=0.3
        assert mtf_upsample != pytest.approx(mtf_downsample)

    # --- Error cases ---

    def test_missing_processing_params_raises(self):
        """get_resampling_mtf_1d raises when processing_params is absent."""
        proc = Processing(name="bare")
        with pytest.raises(ValueError, match="processing_params"):
            proc.get_resampling_mtf_1d(1.0 * u.m)

    def test_missing_resampling_kernel_raises(self):
        """get_resampling_mtf_1d raises when resampling_kernel is not set."""
        proc = Processing(name="no_kernel", processing_params=ProcessingParams())
        with pytest.raises(ValueError, match="resampling_kernel"):
            proc.get_resampling_mtf_1d(1.0 * u.m)

    def test_missing_output_pitch_raises(self):
        """Raises when output_pitch is absent from both YAML and call arg."""
        # ProcessingParams without output_pitch would fail the model_validator,
        # so simulate by calling with no override on a config that has no pitch.
        # We must bypass the validator — patch after construction.
        proc = Processing(name="test", processing_params=ProcessingParams())
        assert proc.processing_params is not None
        proc.processing_params.resampling_kernel = "bilinear"
        with pytest.raises(ValueError, match="output_pitch"):
            proc.get_resampling_mtf_1d(1.0 * u.m)

    # --- ProcessingParams validator ---

    def test_processing_params_kernel_without_pitch_raises(self):
        """Setting resampling_kernel without output_pitch raises at construction."""
        with pytest.raises(ValueError, match="output_pitch"):
            ProcessingParams(resampling_kernel="bilinear")

    def test_processing_params_bad_kernel_params_raises(self):
        """Passing bicubic_a to a non-BICUBIC kernel raises."""
        with pytest.raises(ValueError, match="not used by kernel"):
            ProcessingParams(
                resampling_kernel="bilinear",
                output_pitch="1 m",  # pyright: ignore[reportArgumentType]
                bicubic_a=-0.5,
            )
