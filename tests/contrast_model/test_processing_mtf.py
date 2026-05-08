# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""Tests for the processing-stage MTF backend (resampling kernels)."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from opticks import u
from opticks.contrast_model.mtf import MTF_Model_1D
from opticks.contrast_model.processing_mtf import (
    ResamplingKernel,
    resampling_mtf_1d,
    validate_resampling_params,
)

RESAMPLING_KERNEL_PARAMS = [
    (ResamplingKernel.NEAREST_NEIGHBOR, {}),
    (ResamplingKernel.BILINEAR, {}),
    (ResamplingKernel.BICUBIC, {"bicubic_a": -0.5}),
    (ResamplingKernel.LANCZOS, {"lanczos_n": 3}),
    (ResamplingKernel.SINC, {}),
]


class TestResamplingMTF:
    # --- DC preservation: every kernel must give MTF(0) = 1 ---

    @pytest.mark.parametrize("kernel,extra", RESAMPLING_KERNEL_PARAMS)
    def test_resampling_at_zero_freq(self, kernel, extra):
        """Every kernel preserves DC: MTF(f=0) == 1."""
        val = resampling_mtf_1d(0.0 * u.cy / u.m, kernel, 1.0 * u.m, 1.0 * u.m, **extra)
        assert val == pytest.approx(1.0, abs=1e-9)

    # --- Closed-form known values ---

    def test_resampling_nearest_neighbor_known_value(self):
        """|sinc(nu)| at nu = 0.3."""
        val = resampling_mtf_1d(
            0.3 * u.cy / u.m, ResamplingKernel.NEAREST_NEIGHBOR, 1.0 * u.m, 1.0 * u.m
        )
        assert val == pytest.approx(np.sinc(0.3), rel=1e-12)

    def test_resampling_bilinear_known_value(self):
        """sinc^2(nu) at nu = 0.3."""
        val = resampling_mtf_1d(
            0.3 * u.cy / u.m, ResamplingKernel.BILINEAR, 1.0 * u.m, 1.0 * u.m
        )
        assert val == pytest.approx(np.sinc(0.3) ** 2, rel=1e-12)

    def test_resampling_sinc_passband_and_stopband(self):
        """Brick-wall sinc kernel: MTF=1 below Nyquist, 0 above."""
        freqs = np.array([0.0, 0.25, 0.49, 0.51, 1.0]) * u.cy / u.m
        val = resampling_mtf_1d(freqs, ResamplingKernel.SINC, 1.0 * u.m, 1.0 * u.m)
        assert_allclose(val, [1.0, 1.0, 1.0, 0.0, 0.0])

    def test_resampling_bilinear_equals_nn_squared(self):
        """Identity: sinc^2 == |sinc| * |sinc| over the passband."""
        freqs = np.linspace(0.05, 0.95, 11) * u.cy / u.m
        nn = resampling_mtf_1d(
            freqs, ResamplingKernel.NEAREST_NEIGHBOR, 1.0 * u.m, 1.0 * u.m
        )
        bl = resampling_mtf_1d(freqs, ResamplingKernel.BILINEAR, 1.0 * u.m, 1.0 * u.m)
        assert_allclose(nn**2, bl, rtol=1e-12)

    # --- p_eff = max(input_pitch, output_pitch) rule ---

    def test_resampling_p_eff_uses_output_when_downsample(self):
        """When p_in < p_out (downsample), kernel scales to p_out."""
        # p_in=0.5 m, p_out=1 m -> p_eff=1 m -> nu = f * 1 m
        f = 0.3 * u.cy / u.m
        val = resampling_mtf_1d(f, ResamplingKernel.BILINEAR, 0.5 * u.m, 1.0 * u.m)
        assert val == pytest.approx(np.sinc(0.3) ** 2, rel=1e-12)

    def test_resampling_p_eff_uses_input_when_upsample(self):
        """When p_in > p_out (upsample), kernel scales to p_in."""
        # p_in=2 m, p_out=1 m -> p_eff=2 m -> nu = f * 2 m
        f = 0.3 * u.cy / u.m
        val = resampling_mtf_1d(f, ResamplingKernel.BILINEAR, 2.0 * u.m, 1.0 * u.m)
        assert val == pytest.approx(np.sinc(0.6) ** 2, rel=1e-12)

    def test_resampling_p_eff_equal_pitches(self):
        """When p_in == p_out, kernel scales to that common pitch."""
        f = 0.3 * u.cy / u.m
        val = resampling_mtf_1d(f, ResamplingKernel.BILINEAR, 1.0 * u.m, 1.0 * u.m)
        assert val == pytest.approx(np.sinc(0.3) ** 2, rel=1e-12)

    # --- Sanity / range over a sweep ---

    @pytest.mark.parametrize("kernel,extra", RESAMPLING_KERNEL_PARAMS)
    def test_resampling_non_negative_and_finite(self, kernel, extra):
        """All kernels produce finite, non-negative MTF over the passband and beyond."""
        freqs = np.linspace(0.0, 2.0, 201) * u.cy / u.m
        val = resampling_mtf_1d(freqs, kernel, 1.0 * u.m, 1.0 * u.m, **extra)
        assert np.all(np.isfinite(val))
        assert np.all(val >= 0)

    def test_resampling_lanczos_passband_boost(self):
        """Lanczos exhibits mild MTF boost (>1) in mid passband - this is real."""
        freqs = np.linspace(0.05, 0.45, 41) * u.cy / u.m
        val = resampling_mtf_1d(
            freqs, ResamplingKernel.LANCZOS, 1.0 * u.m, 1.0 * u.m, lanczos_n=3
        )
        # Some values strictly exceed 1 in the Lanczos passband (edge boost).
        assert np.any(val > 1.0)

    def test_resampling_unit_agnostic(self):
        """Equivalent (freq, pitch) pairs in different units give the same MTF.

        nu = f * p_eff is dimensionless; 0.3 cy/m * 1 m = 0.3 cy/mm * 1 mm = 0.3.
        """
        val_m = resampling_mtf_1d(
            0.3 * u.cy / u.m, ResamplingKernel.BILINEAR, 1.0 * u.m, 1.0 * u.m
        )
        val_mm = resampling_mtf_1d(
            0.3 * u.cy / u.mm, ResamplingKernel.BILINEAR, 1.0 * u.mm, 1.0 * u.mm
        )
        assert val_m == pytest.approx(val_mm, rel=1e-12)

    # --- Validators ---

    def test_validate_bicubic_a_forbidden_for_other_kernels(self):
        """bicubic_a may not be supplied with NEAREST_NEIGHBOR."""
        with pytest.raises(ValueError, match="not used by kernel"):
            validate_resampling_params(
                ResamplingKernel.NEAREST_NEIGHBOR,
                bicubic_a=-0.5,
                lanczos_n=None,
            )

    def test_validate_lanczos_n_forbidden_for_other_kernels(self):
        """lanczos_n may not be supplied with BILINEAR."""
        with pytest.raises(ValueError, match="not used by kernel"):
            validate_resampling_params(
                ResamplingKernel.BILINEAR,
                bicubic_a=None,
                lanczos_n=3,
            )

    def test_validate_bicubic_a_out_of_range(self):
        """bicubic_a outside [-1.0, 0.0] raises."""
        with pytest.raises(ValueError, match="bicubic_a"):
            validate_resampling_params(
                ResamplingKernel.BICUBIC,
                bicubic_a=0.5,
                lanczos_n=None,
            )

    def test_validate_lanczos_n_too_small(self):
        """lanczos_n < 2 raises."""
        with pytest.raises(ValueError, match="lanczos_n"):
            validate_resampling_params(
                ResamplingKernel.LANCZOS,
                bicubic_a=None,
                lanczos_n=1,
            )

    def test_validate_lanczos_n_not_integer(self):
        """Non-integer lanczos_n raises."""
        with pytest.raises(ValueError, match="lanczos_n"):
            validate_resampling_params(
                ResamplingKernel.LANCZOS,
                bicubic_a=None,
                lanczos_n=2.5,  # pyright: ignore[reportArgumentType]
            )

    def test_validate_accepts_bare_kernels(self):
        """NN, Bilinear, Sinc accept no params (no error)."""
        for k in (
            ResamplingKernel.NEAREST_NEIGHBOR,
            ResamplingKernel.BILINEAR,
            ResamplingKernel.SINC,
        ):
            validate_resampling_params(k, bicubic_a=None, lanczos_n=None)


class TestResamplingMTFFactory:
    """Tests for MTF_Model_1D.resampling() factory."""

    def test_factory_returns_model_with_dc_equal_one(self):
        """Factory returns a working MTF_Model_1D with MTF(0) = 1."""
        model = MTF_Model_1D.resampling(ResamplingKernel.BILINEAR, 1.0 * u.m, 1.0 * u.m)
        assert model.mtf_value(0.0 * u.cy / u.m) == pytest.approx(1.0, abs=1e-9)

    def test_factory_string_kernel_coercion(self):
        """Passing a kernel label string is equivalent to passing the enum."""
        model_enum = MTF_Model_1D.resampling(
            ResamplingKernel.BILINEAR, 1.0 * u.m, 2.0 * u.m
        )
        model_str = MTF_Model_1D.resampling("bilinear", 1.0 * u.m, 2.0 * u.m)
        freqs = np.array([0.1, 0.3, 0.5]) * u.cy / u.m
        assert_allclose(
            model_enum.mtf_value(freqs), model_str.mtf_value(freqs), rtol=1e-12
        )

    def test_factory_id_string_contains_kernel_and_pitches(self):
        """The id string identifies the kernel and both pitches."""
        model = MTF_Model_1D.resampling(
            ResamplingKernel.LANCZOS, 2.0 * u.m, 1.0 * u.m, lanczos_n=3
        )
        assert "lanczos" in model.id
        assert "2.0 m" in model.id
        assert "1.0 m" in model.id

    def test_factory_matches_backend(self):
        """Factory MTF values match the backend dispatcher exactly."""
        p_in, p_out = 1.5 * u.m, 1.0 * u.m
        freqs = np.array([0.1, 0.25, 0.4]) * u.cy / u.m
        model = MTF_Model_1D.resampling(ResamplingKernel.BICUBIC, p_in, p_out)
        expected = resampling_mtf_1d(freqs, ResamplingKernel.BICUBIC, p_in, p_out)
        assert_allclose(model.mtf_value(freqs), expected, rtol=1e-12)

    def test_factory_rejects_bad_params(self):
        """Factory propagates validator rejection."""
        with pytest.raises(ValueError):
            MTF_Model_1D.resampling(
                ResamplingKernel.BILINEAR, 1.0 * u.m, 1.0 * u.m, bicubic_a=-0.5
            )
