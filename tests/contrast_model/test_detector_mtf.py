# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import numpy as np
import pytest
from astropy.units import Quantity
from numpy.testing import assert_allclose

from opticks import u
from opticks.contrast_model.detector_mtf import (
    DetectorDiffusionModel,
    DetectorDiffusionPreset,
)
from opticks.contrast_model.mtf import MTF_Model_1D
from opticks.imaging_model.detector import Channel, Detector
from tests import process_paths


class TestDetectorMTF:
    pushbr_file_dir = Path("sat_pushbroom_data")
    pushbr_alt_file_dir = Path("tests", "imaging_model", "sat_pushbroom_data")

    input_line_freq: Quantity = 30 * u.cy / u.mm

    @pytest.fixture(scope="class")
    def detector(self) -> Detector:
        file_path = Path("pan_detector.yaml")
        file_path = process_paths(
            file_path, self.pushbr_file_dir, self.pushbr_alt_file_dir
        )
        return Detector.from_yaml_file(file_path)

    def test_mtf_detector_sampling(self, detector: Detector):
        """Tests the detector sampling MTF."""
        truth = 0.7679273089188128

        channel: Channel = detector.channels["pan"]
        mtf_model = MTF_Model_1D.detector_sampling(channel.pixel_pitch())

        assert mtf_model.mtf_value(self.input_line_freq) == pytest.approx(truth, 1e-9)

    # --- Detector diffusion MTF tests ---

    def test_mtf_diffusion_bsi_1(self):
        """Tests the BSI-1 diffusion MTF (dead back surface)."""
        truth = 0.379991599897892

        mtf_model = MTF_Model_1D.detector_diffusion(
            model=DetectorDiffusionModel.BSI_1,
            absorption_coeff=0.19 / u.um,
            diffusion_length=200 * u.um,
            field_free_depth=17 * u.um,
            depletion_depth=15 * u.um,
        )

        assert mtf_model.mtf_value(30 * u.cy / u.mm) == pytest.approx(truth, 1e-9)

    def test_mtf_diffusion_bsi_2(self):
        """Tests the BSI-2 diffusion MTF (passivated back surface)."""
        truth = 0.391074110286271

        mtf_model = MTF_Model_1D.detector_diffusion(
            model=DetectorDiffusionModel.BSI_2,
            absorption_coeff=0.19 / u.um,
            diffusion_length=20 * u.um,
            field_free_depth=5 * u.um,
            depletion_depth=1.5 * u.um,
        )

        assert mtf_model.mtf_value(30 * u.cy / u.mm) == pytest.approx(truth, 1e-9)

    def test_mtf_diffusion_bsi_3(self):
        """Tests the BSI-3 diffusion MTF (thick substrate, general S)."""
        truth = 0.833220434894560

        mtf_model = MTF_Model_1D.detector_diffusion(
            model=DetectorDiffusionModel.BSI_3,
            absorption_coeff=0.8 / u.um,
            diffusion_length=25 * u.um,
            field_free_depth=7 * u.um,
            surface_recomb_velocity=500 * u.cm / u.s,
            diffusion_coeff=15 * u.cm**2 / u.s,
        )

        assert mtf_model.mtf_value(15 * u.cy / u.mm) == pytest.approx(truth, 1e-9)

    def test_mtf_diffusion_fsi_1(self):
        """Tests the FSI-1 diffusion MTF (finite field-free bulk)."""
        truth = 0.991929258182161

        mtf_model = MTF_Model_1D.detector_diffusion(
            model=DetectorDiffusionModel.FSI_1,
            absorption_coeff=0.19 / u.um,
            diffusion_length=15 * u.um,
            field_free_depth=3 * u.um,
            depletion_depth=2 * u.um,
        )

        assert mtf_model.mtf_value(30 * u.cy / u.mm) == pytest.approx(truth, 1e-9)

    def test_mtf_diffusion_fsi_2(self):
        """Tests the FSI-2 diffusion MTF (semi-infinite bulk)."""
        truth = 0.758417884013905

        mtf_model = MTF_Model_1D.detector_diffusion(
            model=DetectorDiffusionModel.FSI_2,
            absorption_coeff=0.19 / u.um,
            diffusion_length=50 * u.um,
            depletion_depth=3 * u.um,
        )

        assert mtf_model.mtf_value(30 * u.cy / u.mm) == pytest.approx(truth, 1e-9)

    def test_mtf_diffusion_at_zero(self):
        """Tests that diffusion MTF approaches 1 at near-zero frequency."""
        mtf_model = MTF_Model_1D.detector_diffusion(
            model=DetectorDiffusionModel.BSI_1,
            absorption_coeff=0.19 / u.um,
            diffusion_length=200 * u.um,
            field_free_depth=17 * u.um,
            depletion_depth=15 * u.um,
        )

        assert mtf_model.mtf_value(0.001 * u.cy / u.mm) == pytest.approx(1.0, abs=1e-6)

    def test_mtf_diffusion_in_range(self):
        """Tests that all presets produce MTF in [0, 1] at a typical spatial frequency."""
        alpha = 0.19 / u.um  # silicon at ~700 nm

        for preset in DetectorDiffusionPreset:
            mtf_model = MTF_Model_1D.detector_diffusion_preset(
                preset=preset,
                absorption_coeff=alpha,
            )
            val = mtf_model.mtf_value(30 * u.cy / u.mm)
            assert 0 <= val <= 1, f"MTF out of range for preset {preset}: {val}"

    def test_mtf_diffusion_preset_matches_direct(self):
        """Tests that the SCIENTIFIC_CCD preset matches a direct parameter call."""
        alpha = 0.19 / u.um

        preset_model = MTF_Model_1D.detector_diffusion_preset(
            preset=DetectorDiffusionPreset.SCIENTIFIC_CCD,
            absorption_coeff=alpha,
        )
        direct_model = MTF_Model_1D.detector_diffusion(
            model=DetectorDiffusionModel.BSI_1,
            absorption_coeff=alpha,
            diffusion_length=200 * u.um,
            field_free_depth=17 * u.um,
            depletion_depth=15 * u.um,
        )

        freq = 30 * u.cy / u.mm
        assert preset_model.mtf_value(freq) == pytest.approx(
            direct_model.mtf_value(freq), 1e-9
        )

    def test_mtf_diffusion_preset_override(self):
        """Tests that overriding a preset parameter changes the result."""
        alpha = 0.19 / u.um
        freq = 30 * u.cy / u.mm

        default_model = MTF_Model_1D.detector_diffusion_preset(
            preset=DetectorDiffusionPreset.SCIENTIFIC_CCD,
            absorption_coeff=alpha,
        )
        override_model = MTF_Model_1D.detector_diffusion_preset(
            preset=DetectorDiffusionPreset.SCIENTIFIC_CCD,
            absorption_coeff=alpha,
            diffusion_length=150 * u.um,
        )

        assert default_model.mtf_value(freq) != override_model.mtf_value(freq)

    def test_mtf_diffusion_missing_params(self):
        """Tests that omitting a required parameter raises ValueError."""
        with pytest.raises(ValueError):
            MTF_Model_1D.detector_diffusion(
                model=DetectorDiffusionModel.BSI_3,
                absorption_coeff=0.8 / u.um,
                diffusion_length=25 * u.um,
                field_free_depth=7 * u.um,
                # surface_recomb_velocity omitted — required for BSI_3
                diffusion_coeff=15 * u.cm**2 / u.s,
            )

    def test_mtf_diffusion_forbidden_params(self):
        """Tests that supplying a forbidden parameter raises ValueError."""
        with pytest.raises(ValueError):
            MTF_Model_1D.detector_diffusion(
                model=DetectorDiffusionModel.BSI_1,
                absorption_coeff=0.19 / u.um,
                diffusion_length=200 * u.um,
                field_free_depth=17 * u.um,
                depletion_depth=15 * u.um,
                surface_recomb_velocity=500 * u.cm / u.s,  # forbidden for BSI_1
            )

    def test_mtf_diffusion_preset_model_mismatch(self):
        """Tests that a forbidden override for the preset's model raises ValueError."""
        with pytest.raises(ValueError):
            MTF_Model_1D.detector_diffusion_preset(
                preset=DetectorDiffusionPreset.SCIENTIFIC_CCD,  # BSI-1
                absorption_coeff=0.19 / u.um,
                surface_recomb_velocity=500 * u.cm / u.s,  # forbidden for BSI_1
            )

    def test_mtf_diffusion_range_check(self):
        """Tests that a physically unrealistic parameter set raises ValueError."""
        # BSI-2 with tiny L_b causes eta0 < 0, making MTF > 1 after normalisation
        with pytest.raises(ValueError):
            mtf_model = MTF_Model_1D.detector_diffusion(
                model=DetectorDiffusionModel.BSI_2,
                absorption_coeff=0.19 / u.um,
                diffusion_length=200 * u.um,
                field_free_depth=17 * u.um,
                depletion_depth=0.001 * u.um,  # unrealistically small → eta0 < 0
            )
            mtf_model.mtf_value(100 * u.cy / u.mm)

    # --- Detector crosstalk MTF tests ---

    def test_mtf_crosstalk_center(self):
        """Tests center-pixel crosstalk MTF with side and diagonal coupling."""
        # Analytic 1D formula: 1 - 2*(Xs + 2*Xd)*(1 - cos(2*pi*f*p))
        # evaluated at Xs=0.03, Xd=0.005, p=13 um, f=30 cy/mm.
        truth = 0.858358940577937

        mtf_model = MTF_Model_1D.detector_crosstalk(
            pixel_pitch=13 * u.um,
            crosstalk_xs=0.03,
            crosstalk_xd=0.005,
        )

        assert mtf_model.mtf_value(30 * u.cy / u.mm) == pytest.approx(truth, 1e-9)

    def test_mtf_crosstalk_xd_zero_matches_1d(self):
        """When xd=0, the model matches 1 - 2*Xs*(1 - cos(2*pi*f*p))."""
        xs = 0.03
        pitch = 13 * u.um
        p_mm = pitch.to(u.mm).value  # pyright: ignore[reportAttributeAccessIssue]
        freqs = np.array([10, 20, 30, 38.46]) * u.cy / u.mm

        mtf_model = MTF_Model_1D.detector_crosstalk(
            pixel_pitch=pitch,
            crosstalk_xs=xs,
        )

        expected = 1 - 2 * xs * (1 - np.cos(2 * np.pi * freqs.value * p_mm))
        result = mtf_model.mtf_value(freqs)
        assert_allclose(result, expected, rtol=1e-12)

    def test_mtf_crosstalk_at_zero_freq(self):
        """Crosstalk MTF is 1.0 at zero frequency (signal conserved)."""
        mtf_model = MTF_Model_1D.detector_crosstalk(
            pixel_pitch=13 * u.um,
            crosstalk_xs=0.03,
            crosstalk_xd=0.005,
        )

        assert mtf_model.mtf_value(0 * u.cy / u.mm) == pytest.approx(1.0)

    def test_mtf_crosstalk_at_nyquist(self):
        """At Nyquist (f=1/2p), MTF = 1 - 4*(Xs + 2*Xd)."""
        xs, xd = 0.03, 0.005
        pitch = 13 * u.um
        f_ny = 1 / (2 * pitch.to(u.mm).value) * u.cy / u.mm  # pyright: ignore[reportAttributeAccessIssue]

        mtf_model = MTF_Model_1D.detector_crosstalk(
            pixel_pitch=pitch,
            crosstalk_xs=xs,
            crosstalk_xd=xd,
        )

        expected = 1 - 4 * (xs + 2 * xd)
        assert mtf_model.mtf_value(f_ny) == pytest.approx(expected, 1e-12)

    def test_mtf_crosstalk_in_range(self):
        """Crosstalk MTF stays in [0, 1] for typical detector parameters."""
        freqs = np.linspace(0, 40, 100) * u.cy / u.mm

        mtf_model = MTF_Model_1D.detector_crosstalk(
            pixel_pitch=13 * u.um,
            crosstalk_xs=0.04,
            crosstalk_xd=0.01,
        )

        result = mtf_model.mtf_value(freqs)
        assert np.all(result >= 0)
        assert np.all(result <= 1)

    def test_mtf_crosstalk_invalid_negative_xs(self):
        """Negative xs raises ValueError."""
        with pytest.raises(ValueError):
            MTF_Model_1D.detector_crosstalk(
                pixel_pitch=13 * u.um,
                crosstalk_xs=-0.01,
            )

    def test_mtf_crosstalk_invalid_negative_xd(self):
        """Negative xd raises ValueError."""
        with pytest.raises(ValueError):
            MTF_Model_1D.detector_crosstalk(
                pixel_pitch=13 * u.um,
                crosstalk_xs=0.03,
                crosstalk_xd=-0.01,
            )

    def test_mtf_crosstalk_excessive_coupling(self):
        """Excessive coupling (4*xs + 4*xd >= 1) raises ValueError."""
        with pytest.raises(ValueError):
            MTF_Model_1D.detector_crosstalk(
                pixel_pitch=13 * u.um,
                crosstalk_xs=0.2,
                crosstalk_xd=0.1,
            )
