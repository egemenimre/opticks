# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import numpy as np
import pytest
from prysm import _richdata, coordinates, geometry, polynomials, propagation
from prysm.coordinates import cart_to_polar, make_xy_grid
from prysm.geometry import circle
from prysm.fttools import pad2d, mdft
from prysm.propagation import focus
from prysm._richdata import RichData

from opticks import process_paths, u
from opticks.imager_model.optics import Aperture, Optics
from opticks.utils.prysm_utils import OptPathDiff
from opticks.utils.testing_utils import assert_allclose


class TestOptics:

    @pytest.fixture(scope="class")
    def optics(self) -> Optics:

        file_directory = Path("sat_pushbroom_data")
        alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
        file_path = Path("optics.yaml")

        # different test environments work with different paths
        file_path = process_paths(file_path, file_directory, alt_file_directory)

        return Optics.from_yaml_file(file_path)

    def test_fov(self, optics):
        """Tests the full optical FoV."""

        truth = 1.7757828128191897 * u.deg
        fov = optics.full_optical_fov

        assert_allclose(fov, truth, atol=0.001 * u.mdeg)

    def test_spatial_cutoff(self, optics):
        """Tests the spatial cutoff frequency."""

        # set up
        ref_wavelength = 640 * u.nm

        truth = 78.70011623401781 * u.cy / u.mm

        # computation
        cutoff_freq = optics.spatial_cutoff_freq(ref_wavelength)

        # verification
        assert_allclose(cutoff_freq, truth, atol=0.00001 * u.cy / u.mm)

    def test_circle_aperture(self, optics):
        """Tests the circle aperture."""
        diameter = optics.params.aperture_diameter.m_as(u.mm)

        samples = 256

        # prysm definition
        xi, eta = make_xy_grid(samples, diameter=diameter)
        r, t = cart_to_polar(xi, eta)

        aperture_prysm = circle(diameter / 2, r)

        # computation
        aperture = Aperture.circle_aperture(
            optics.params.aperture_diameter, samples, with_units=True
        )

        # verification
        np.testing.assert_array_equal(aperture.data, aperture_prysm, strict=True)

    def test_circle_aperture_with_obscuration(self, optics):
        """Tests the circle aperture."""
        diameter = optics.params.aperture_diameter.m_as(u.mm)

        samples = 256
        obscuration_ratio = 0.3

        # prysm definition
        xi, eta = make_xy_grid(samples, diameter=diameter)
        r, t = cart_to_polar(xi, eta)

        pm_od = circle(diameter / 2, r)
        pm_id = circle(diameter / 2 * obscuration_ratio, r)
        aperture_prysm = pm_od ^ pm_id  # or pm_od & ~pm_id

        # computation
        aperture = Aperture.circle_aperture_with_obscuration(
            optics.params.aperture_diameter, obscuration_ratio, samples, with_units=True
        )

        # verification
        np.testing.assert_array_equal(aperture.data, aperture_prysm, strict=True)

    def test_aperture_scaling(self):
        """Replicate the 'Radiometrically Correct Modeling' example
        in prysm."""

        # prysm definition
        x, y = make_xy_grid(256, diameter=2)
        r, t = cart_to_polar(x, y)
        aperture = circle(1, r)
        inc_psf = abs(focus(aperture, Q=2)) ** 2

        # unity sum
        aperture2 = aperture / np.sqrt(aperture.sum())
        inc_psf_norm_sum = abs(focus(aperture2, Q=2)) ** 2

        # unity peak
        aperture3 = pad2d(aperture, Q=2)
        aperture3 = aperture3 * (2 * np.sqrt(aperture.size) / aperture.sum())
        inc_psf_norm_peak = abs(focus(aperture3, Q=1)) ** 2
        inc_psf_norm_peak.sum(), inc_psf_norm_peak.max()

        # computation
        yaml_text = """
        name: PTest
        focal_length: 2000 mm
        aperture_diameter: 500 mm
        image_diam_on_focal_plane: 50 mm   # not used"""

        optics = Optics.from_yaml_text(yaml_text)

        pup_samples = 256

        aperture = Aperture.circle_aperture(
            optics.params.aperture_diameter, pup_samples
        )

        optics.set_aperture_model(aperture)

        ref_wvl = 500 * u.nm

        optics.add_pupil_func(ref_wvl, None)

        Q = 2
        psf_dx = ref_wvl.to(u.um) * optics.f_number / Q

        psf_samples = pup_samples * Q
        psf = optics.psf(ref_wvl, psf_dx, psf_samples=psf_samples)

        # unity sum
        optics_psf_unity_sum = Optics.from_yaml_text(yaml_text)

        aperture_unity_sum = aperture.scale_for_norm_sum_psf()
        optics_psf_unity_sum.set_aperture_model(aperture_unity_sum)

        optics_psf_unity_sum.add_pupil_func(ref_wvl, None)

        psf_unity_sum = optics_psf_unity_sum.psf(
            ref_wvl, psf_dx, psf_samples=psf_samples
        )

        # unity peak
        optics_psf_unity_peak = Optics.from_yaml_text(yaml_text)

        Q_pad = 1
        aperture_unity_peak = aperture.scale_for_norm_peak_psf(Q, Q_pad)
        optics_psf_unity_peak.set_aperture_model(aperture_unity_peak)

        optics_psf_unity_peak.add_pupil_func(ref_wvl, None)

        psf_unity_peak = optics_psf_unity_peak.psf(
            ref_wvl, psf_dx, psf_samples=psf_samples
        )

        # verification
        assert_allclose(
            [psf.data.sum(), psf.data.max()], [inc_psf.sum(), inc_psf.max()], rtol=1e-9
        )
        assert_allclose(
            [psf_unity_sum.data.sum(), psf_unity_sum.data.max()],
            [inc_psf_norm_sum.sum(), inc_psf_norm_sum.max()],
            rtol=1e-9,
        )
        assert_allclose(
            [psf_unity_peak.data.sum(), psf_unity_peak.data.max()],
            [inc_psf_norm_peak.sum(), inc_psf_norm_peak.max()],
            rtol=1e-9,
        )

    def test_psf(self):
        """Replicate the 'polychromatic propagation'
        example in prysm docs (but with fixed sampling)."""

        # prysm
        # ------
        # prysm set up
        ap_samples = 256
        res = 512  # psf_samples
        fno = 4  # F number
        efl = 150  # effective focal length
        epd = efl / fno  # aperture diameter
        r_aper = epd / 2  # aperture radius
        wvl0 = 0.550  # wavelength in microns

        # sample size in microns
        res_el = wvl0 * fno * 1.22 / 4  # 4 pixels per airy radius

        xi, eta = coordinates.make_xy_grid(ap_samples, diameter=epd)
        r, t = coordinates.cart_to_polar(xi, eta)
        dx = xi[0, 1] - xi[0, 0]

        r_aber = r / r_aper

        coef = wvl0 * 1e3 * 15  # 15 waves of defocus
        phs = polynomials.hopkins(0, 2, 0, r_aber, t, 1) * coef

        amp = geometry.circle(r_aper, r)

        # monochromatic psf
        wf = propagation.Wavefront.from_amp_and_phase(amp, phs, wvl0, dx)
        psf_mono_prysm = wf.focus_fixed_sampling(efl, res_el, res).intensity

        # polychromatic psf
        halfbw = 0.2
        wvls = np.linspace(
            wvl0 * (1 - halfbw), wvl0 * (1 + halfbw), 11
        )  # 11 discrete wavelengths
        spectral_weights = np.ones_like(wvls)

        components = []
        for wvl in wvls:
            wf = propagation.Wavefront.from_amp_and_phase(amp, phs, wvl, dx)
            focused = wf.focus_fixed_sampling(
                efl, res_el, res
            )  # 512 samples in the output domain
            components.append(
                focused.intensity.data
            )  # sum of intensities, wvls are incoherent to each other

        # psf is just an array
        psf_poly_prysm = polynomials.sum_of_2d_modes(components, spectral_weights)
        # until we enrich it
        psf_poly_prysm = _richdata.RichData(psf_poly_prysm, res_el, wvl0)

        # opticks
        # -------
        # opticks set up
        yaml_text = """
        name: Prysm Test data
        focal_length: 150 mm
        aperture_diameter: 37.5 mm
        image_diam_on_focal_plane: 10 mm # unknown"""

        optics = Optics.from_yaml_text(yaml_text)

        aperture = Aperture.circle_aperture(
            optics.params.aperture_diameter, ap_samples, with_units=True
        )
        optics.set_aperture_model(aperture)

        # compute opd
        opd = OptPathDiff(phs * u.nm)

        # generate and add the Wavefront or Pupil function (mono)
        optics.add_pupil_func(wvl0 * u.um, opd)

        # compute monochromatic PSF
        psf_mono = optics.psf(
            wvl0 * u.um, res_el * u.um, res, spectral_weights=None, with_units=True
        )

        # polychromatic psf

        # reset pupil functs
        optics.pupils = []

        # reuse the same opd
        for wvl in wvls:
            optics.add_pupil_func(wvl * u.um, opd)

        # compute polychromatic PSF
        psf_poly = optics.psf(
            wvl0 * u.um, res_el * u.um, res, spectral_weights=None, with_units=True
        )

        # verification
        np.testing.assert_allclose(psf_mono.data, psf_mono_prysm.data, rtol=1e-10)
        assert_allclose(psf_mono.dx, psf_mono_prysm.dx * u.um, rtol=1e-10)
        np.testing.assert_allclose(psf_poly.data, psf_poly_prysm.data, rtol=1e-10)
        assert_allclose(psf_poly.dx, psf_poly_prysm.dx * u.um, rtol=1e-10)
