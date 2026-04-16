# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import numpy as np
import pytest
from astropy.units import Quantity
from astropy.units import isclose as qty_isclose

from opticks import u
from opticks.contrast_model.mtf import MTF_Model_1D
from opticks.contrast_model.optics_mtf import FieldAberrationModel
from opticks.imager_model.optics import Optics
from tests import process_paths


class TestFieldAberrationModel:
    """Tests for the Seidel field-dependent aberration model."""

    ref_wvl: Quantity = 550 * u.nm

    # A representative coefficient set (in waves)
    w040 = 0.05
    w131 = 0.10
    w222 = 0.08
    w220 = 0.06

    @pytest.fixture(scope="class")
    def model_waves(self) -> FieldAberrationModel:
        return FieldAberrationModel.from_waves(
            self.w040, self.w131, self.w222, self.w220, self.ref_wvl
        )

    @pytest.fixture(scope="class")
    def model_quantity(self) -> FieldAberrationModel:
        ref = self.ref_wvl.to(u.nm)
        return FieldAberrationModel.from_quantity(
            self.w040 * ref,
            self.w131 * ref,
            self.w222 * ref,
            self.w220 * ref,
            self.ref_wvl,
        )

    # ---- Dual constructors produce identical state ----

    def test_dual_constructors_match(self, model_waves, model_quantity):
        """from_waves and from_quantity produce identical internal state."""
        assert qty_isclose(model_waves.w040, model_quantity.w040)
        assert qty_isclose(model_waves.w131_edge, model_quantity.w131_edge)
        assert qty_isclose(model_waves.w222_edge, model_quantity.w222_edge)
        assert qty_isclose(model_waves.w220_edge, model_quantity.w220_edge)
        assert qty_isclose(
            model_waves.reference_wavelength, model_quantity.reference_wavelength
        )

    # ---- w_rms behaviour ----

    def test_w_rms_at_zero_is_spherical_only(self, model_waves):
        """At h=0, only spherical aberration contributes."""
        from opticks.contrast_model.optics_mtf import _RMS_COEFF_SPHERICAL

        expected = abs(_RMS_COEFF_SPHERICAL * model_waves.w040)
        result = model_waves.w_rms(0.0)

        assert qty_isclose(result, expected, rtol=1e-12)

    def test_w_rms_at_edge_uses_all_terms(self, model_waves):
        """At h=1, all four terms contribute."""
        from opticks.contrast_model.optics_mtf import (
            _RMS_COEFF_ASTIGMATISM,
            _RMS_COEFF_COMA,
            _RMS_COEFF_FIELD_CURVATURE,
            _RMS_COEFF_SPHERICAL,
        )

        var = (
            (_RMS_COEFF_SPHERICAL * model_waves.w040) ** 2
            + (_RMS_COEFF_COMA * model_waves.w131_edge) ** 2
            + (_RMS_COEFF_ASTIGMATISM * model_waves.w222_edge) ** 2
            + (_RMS_COEFF_FIELD_CURVATURE * model_waves.w220_edge) ** 2
        )
        expected = np.sqrt(var).to(u.nm)
        result = model_waves.w_rms(1.0)

        assert qty_isclose(result, expected, rtol=1e-12)

    def test_w_rms_monotonically_nondecreasing(self, model_waves):
        """w_rms is monotonically non-decreasing with h."""
        h_values = np.linspace(0, 1, 50)
        rms_values = model_waves.w_rms(h_values).value

        diffs = np.diff(rms_values)
        assert np.all(diffs >= -1e-15), "w_rms is not monotonically non-decreasing"

    # ---- w_rms_waves ----

    def test_w_rms_waves_consistent(self, model_waves):
        """w_rms_waves equals w_rms / reference_wavelength."""
        h = 0.7
        expected = (
            (model_waves.w_rms(h) / model_waves.reference_wavelength).decompose().value
        )
        result = model_waves.w_rms_waves(h)

        assert result == pytest.approx(expected, rel=1e-12)

    def test_w_rms_waves_custom_wavelength(self, model_waves):
        """w_rms_waves with explicit wavelength differs from default."""
        h = 0.5
        default = model_waves.w_rms_waves(h)
        custom = model_waves.w_rms_waves(h, wavelength=700 * u.nm)

        # Different wavelengths should give different results
        assert default != pytest.approx(custom, rel=1e-6)

    # ---- to_zernikes ----

    def test_zernikes_on_axis_only_spherical(self, model_waves):
        """At (0,0), only spherical (Z11) and its defocus balance (Z4) are set."""
        z = model_waves.to_zernikes(0.0, 0.0)

        # Z11 = W_040
        assert qty_isclose(z[11], model_waves.w040, rtol=1e-12)

        # Z4 should have the defocus balance for spherical only
        # (no astigmatism or field curvature at h=0)
        assert qty_isclose(z[4], -model_waves.w040, rtol=1e-12)

        # Coma (Z7, Z8), astigmatism (Z5, Z6) should be zero
        assert qty_isclose(z[7], 0 * u.nm, atol=1e-15 * u.nm)
        assert qty_isclose(z[8], 0 * u.nm, atol=1e-15 * u.nm)
        assert qty_isclose(z[5], 0 * u.nm, atol=1e-15 * u.nm)
        assert qty_isclose(z[6], 0 * u.nm, atol=1e-15 * u.nm)

    def test_zernikes_azimuth_rotation(self, model_waves):
        """Rotating the field point by 90° rotates coma and astigmatism."""
        h = 0.6
        z_x = model_waves.to_zernikes(h, 0.0)
        z_y = model_waves.to_zernikes(0.0, h)

        # Coma: Z7(x-field) should match Z8(y-field), Z8(x-field) ≈ 0
        assert qty_isclose(z_x[7], z_y[8], rtol=1e-12)
        assert qty_isclose(z_x[8], 0 * u.nm, atol=1e-12 * u.nm)
        assert qty_isclose(z_y[7], 0 * u.nm, atol=1e-12 * u.nm)

    def test_zernikes_length(self, model_waves):
        """Returned vector has the requested number of terms."""
        z = model_waves.to_zernikes(0.5, 0.3, n_terms=20)
        assert len(z) == 20

    # ---- Integration: Optics.field_mtf_model_1d ----

    pushbr_file_dir = Path("sat_pushbroom_data")
    pushbr_alt_file_dir = Path("tests", "imager_model", "sat_pushbroom_data")

    @pytest.fixture(scope="class")
    def optics(self) -> Optics:
        file_path = Path("optics.yaml")
        file_path = process_paths(
            file_path, self.pushbr_file_dir, self.pushbr_alt_file_dir
        )
        return Optics.from_yaml_file(file_path)

    def test_field_mtf_zero_coeffs_is_ideal(self, optics):
        """With all Seidel coefficients = 0, field MTF equals ideal optics."""
        wvl = 650 * u.nm
        zero_model = FieldAberrationModel.from_waves(0, 0, 0, 0, wvl)

        field_mtf = optics.field_mtf_model_1d(zero_model, h=1.0, wavelength=wvl)
        ideal_mtf = MTF_Model_1D.ideal_optics(optics.spatial_cutoff_freq(wvl), wvl)

        freq = 30 * u.cy / u.mm
        assert field_mtf.mtf_value(freq) == pytest.approx(
            ideal_mtf.mtf_value(freq), rel=1e-9
        )

    def test_field_mtf_on_axis_is_spherical_only(self, optics):
        """At h=0, the result should match emp_model with spherical-only w_rms."""
        wvl = 650 * u.nm
        model = FieldAberrationModel.from_waves(
            self.w040, self.w131, self.w222, self.w220, wvl
        )

        field_mtf = optics.field_mtf_model_1d(model, h=0.0, wavelength=wvl)
        w_rms_waves = model.w_rms_waves(0.0, wvl)
        direct_mtf = MTF_Model_1D.emp_model_aberrated_optics(
            optics.spatial_cutoff_freq(wvl), w_rms=w_rms_waves, wavelength=wvl
        )

        freq = 30 * u.cy / u.mm
        assert field_mtf.mtf_value(freq) == pytest.approx(
            direct_mtf.mtf_value(freq), rel=1e-12
        )

    def test_field_mtf_monotonically_nonincreasing_in_h(self, optics):
        """MTF at a fixed frequency is monotonically non-increasing in h."""
        wvl = 650 * u.nm
        model = FieldAberrationModel.from_waves(
            self.w040, self.w131, self.w222, self.w220, wvl
        )

        freq = 30 * u.cy / u.mm
        h_values = [0.0, 0.5, 1.0]
        mtf_values = [
            optics.field_mtf_model_1d(model, h=h, wavelength=wvl).mtf_value(freq)
            for h in h_values
        ]

        for i in range(len(mtf_values) - 1):
            assert mtf_values[i] >= mtf_values[i + 1] - 1e-15, (
                f"MTF not non-increasing: h={h_values[i]}: {mtf_values[i]}, "
                f"h={h_values[i + 1]}: {mtf_values[i + 1]}"
            )

    # ---- Field-coordinate helpers ----

    def test_normalised_field_from_angle_axis_is_zero(self, optics):
        """A zero-angle field point maps to (0, 0)."""
        h_x, h_y = optics.normalised_field_from_angle(0 * u.deg, 0 * u.deg)
        assert h_x == pytest.approx(0.0, abs=1e-15)
        assert h_y == pytest.approx(0.0, abs=1e-15)

    def test_normalised_field_from_angle_edge_is_unity(self, optics):
        """The half-FoV angle maps to H = 1 along that axis."""
        half_fov = optics.full_optical_fov / 2.0
        h_x, h_y = optics.normalised_field_from_angle(half_fov, 0 * u.deg)
        assert h_x == pytest.approx(1.0, rel=1e-12)
        assert h_y == pytest.approx(0.0, abs=1e-15)

    def test_normalised_field_from_angle_signed(self, optics):
        """Negative angles produce negative normalised coordinates."""
        half_fov = optics.full_optical_fov / 2.0
        h_x, h_y = optics.normalised_field_from_angle(-half_fov / 2, half_fov / 2)
        assert h_x == pytest.approx(-0.5, rel=1e-12)
        assert h_y == pytest.approx(0.5, rel=1e-12)

    def test_normalised_field_from_focal_plane_xy_axis_is_zero(self, optics):
        """The optical-axis intercept maps to (0, 0)."""
        h_x, h_y = optics.normalised_field_from_focal_plane_xy(0 * u.mm, 0 * u.mm)
        assert h_x == pytest.approx(0.0, abs=1e-15)
        assert h_y == pytest.approx(0.0, abs=1e-15)

    def test_normalised_field_from_focal_plane_xy_edge_is_unity(self, optics):
        """The image radius on the focal plane maps to H = 1."""
        r_max = optics.image_diam_on_focal_plane / 2.0
        h_x, h_y = optics.normalised_field_from_focal_plane_xy(r_max, 0 * u.mm)
        assert h_x == pytest.approx(1.0, rel=1e-12)
        assert h_y == pytest.approx(0.0, abs=1e-15)

    def test_normalised_field_from_focal_plane_xy_signed(self, optics):
        """Negative coordinates produce negative normalised values."""
        r_max = optics.image_diam_on_focal_plane / 2.0
        h_x, h_y = optics.normalised_field_from_focal_plane_xy(-r_max / 2, r_max / 2)
        assert h_x == pytest.approx(-0.5, rel=1e-12)
        assert h_y == pytest.approx(0.5, rel=1e-12)
