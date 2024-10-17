# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import numpy as np
import pytest
from prysm.coordinates import cart_to_polar, make_xy_grid
from prysm.geometry import regular_polygon
from prysm.otf import mtf_from_psf
from prysm.polynomials import zernike_nm
from prysm.propagation import Wavefront

from opticks import u
from opticks.imager_model.detector import Channel, Detector
from opticks.imager_model.optics import Optics
from opticks.perf_model.mtf import MTF_Model_1D
from tests import process_paths


class TestMTF:

    pushbr_file_dir = Path("sat_pushbroom_data")
    pushbr_alt_file_dir = Path("tests", "imager_model", "sat_pushbroom_data")

    perf_model_file_dir = Path("data")
    perf_model_alt_file_dir = Path("tests", "perf_model", "data")

    ref_wavelength = 650 * u.nm
    input_line_freq = 30 * u.cy / u.mm

    @pytest.fixture(scope="class")
    def optics(self) -> Optics:

        file_path = Path("optics.yaml")

        # different test environments work with different paths
        file_path = process_paths(
            file_path, self.pushbr_file_dir, self.pushbr_alt_file_dir
        )

        return Optics.from_yaml_file(file_path)

    @pytest.fixture(scope="class")
    def detector(self) -> Detector:

        file_path = Path("pan_detector.yaml")

        # different test environments work with different paths
        file_path = process_paths(
            file_path, self.pushbr_file_dir, self.pushbr_alt_file_dir
        )

        return Detector.from_yaml_file(file_path)

    def test_mtf_ext_data(self):
        """Tests the external MTF data."""

        file_path = Path("mtf.txt")

        # different test environments work with different paths
        file_path = process_paths(
            file_path, self.perf_model_file_dir, self.perf_model_alt_file_dir
        )

        input_line_freq = 32.456 * u.cy / u.mm

        # load external data file
        freq_data, mtf_data = np.loadtxt(file_path, unpack=True)

        # check values
        truth = 0.4464741694117462

        # Generate the MTF model and values
        mtf_model = MTF_Model_1D.external_data(freq_data, mtf_data)

        mtf_value = mtf_model.mtf_value(input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_fixed_value(self):
        """Tests the fixed value MTF."""

        # check values
        truth = 0.85

        mtf_value = 0.85
        # Generate the MTF model and values
        mtf_model = MTF_Model_1D.fixed(mtf_value)

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_fixed_value_err_neg(self):
        with pytest.raises(ValueError):

            mtf_value = -0.5
            # Generate the MTF model and values
            MTF_Model_1D.fixed(mtf_value)

    def test_mtf_ideal_optics(self, optics):
        """Tests the ideal optics MTF."""

        # check values
        truth = 0.5196720931409163

        # Generate the MTF model and values
        mtf_model = MTF_Model_1D.ideal_optics(self.ref_wavelength, optics)

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_aberrated_optics(self, optics: Optics):
        """Tests the aberrated optics MTF."""

        # check values
        truth = 0.4816165574868274

        # Generate the MTF model and values
        w_rms = 0.05
        mtf_model = MTF_Model_1D.emp_model_aberrated_optics(
            self.ref_wavelength, w_rms, optics
        )

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_detector_sampling(self, detector: Detector):
        """Tests the detector sampling MTF."""

        # check values
        truth = 0.7679273089188128

        # select the pan channel
        channel: Channel = detector.params.channels.pan

        # Generate the MTF model and values
        mtf_model = MTF_Model_1D.detector_sampling(channel.pixel_pitch())

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_jitter(self, detector: Detector):
        """Tests the jitter MTF."""

        # check values
        truth = 0.9704228869250533

        # select the pan channel
        channel: Channel = detector.params.channels.pan

        # Generate the MTF model and values
        stdev_jitter = 0.1  # 10% of the pix
        mtf_model = MTF_Model_1D.jitter(channel.pixel_pitch(), stdev_jitter)

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_smear(self, detector: Detector):
        """Tests the drift/smear MTF."""

        # check values
        truth = 0.9776341205410619

        # select the pan channel
        channel: Channel = detector.params.channels.pan

        # Generate the MTF model and values
        blur_extent = 0.3  # 30% of the pix
        mtf_model = MTF_Model_1D.smear(channel.pixel_pitch(), blur_extent)

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_combined(self, optics: Optics, detector: Detector):
        """Tests the combined MTF."""

        # check values
        truth = 0.3990703920059105

        # select the pan channel
        channel: Channel = detector.params.channels.pan

        # Generate the MTF model and values
        mtf_model_1 = MTF_Model_1D.ideal_optics(self.ref_wavelength, optics)
        mtf_model_2 = MTF_Model_1D.detector_sampling(channel.pixel_pitch())

        mtf_model = MTF_Model_1D.combined(mtf_model_1, mtf_model_2)

        # mtf_value_1 = mtf_model_1.mtf_value(self.input_line_freq)
        # mtf_value_2 = mtf_model_2.mtf_value(self.input_line_freq)

        mtf_value = mtf_model.mtf_value(self.input_line_freq)

        # verification
        assert mtf_value == pytest.approx(truth, 1e-9)

    def test_mtf_from_2d(self):
        """Tests the prysm based 2D MTF."""

        # prysm mtf computation
        efl = 50
        fno = 8

        ap_samples = 256

        x, y = make_xy_grid(ap_samples, diameter=efl / fno)
        dx = x[0, 1] - x[0, 0]
        r, t = cart_to_polar(x, y)

        radius = efl / fno / 2
        rho = r / radius
        n_sides = 14

        poly_aperture = regular_polygon(n_sides, radius, x, y)

        wvl = 0.55  # mid visible band, um

        # spherical aberration
        wfe_nm_rms = wvl / 14 * 1e3  # nm, 1/14 of a wave, 1e3 = um to nm
        mode_1 = zernike_nm(4, 0, rho, t)
        opd_1 = mode_1 * wfe_nm_rms
        pup_1 = Wavefront.from_amp_and_phase(poly_aperture, opd_1, wvl, dx)
        coherent_psf = pup_1.focus(efl, Q=2)

        psf = coherent_psf.intensity
        mtf_1 = mtf_from_psf(psf)

        fx_1, mtf_2d_x_1 = mtf_1.slices(twosided=False).x

        # coma
        mode_2 = zernike_nm(3, 1, rho, t)  # only this line changed
        opd_2 = mode_2 * wfe_nm_rms
        pup_2 = Wavefront.from_amp_and_phase(poly_aperture, opd_2, wvl, dx)
        coherent_psf = pup_2.focus(efl, Q=2)
        psf = coherent_psf.intensity
        mtf_2 = mtf_from_psf(psf, psf.dx)

        fx_2, mtf_2d_x_2 = mtf_2.slices(twosided=False).x

        # opticks computation
        mtf_opt_1 = MTF_Model_1D.from_mtf_2d(mtf_1, "x")
        mtf_opt_2 = MTF_Model_1D.from_mtf_2d(mtf_2, "x")

        # verification
        np.testing.assert_allclose(mtf_opt_1.mtf_value(fx_1), mtf_2d_x_1, rtol=1e-10)
        np.testing.assert_allclose(mtf_opt_2.mtf_value(fx_1), mtf_2d_x_2, rtol=1e-10)
