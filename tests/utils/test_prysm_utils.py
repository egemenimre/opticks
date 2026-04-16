# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import numpy as np
import pytest
from prysm.coordinates import cart_to_polar, make_xy_grid
from prysm.polynomials import (
    ansi_j_to_nm,
    fringe_to_nm,
    noll_to_nm,
    sum_of_2d_modes,
    zernike_nm_seq,
)

from opticks import u
from opticks.utils.prysm_utils import Grid, OptPathDiff


class TestPrysmUtils:
    def test_zernike(self):
        """Replicate the 'Q2d_sequence' example in
        the prysm docs 'Image Simulation' with Zernikes."""

        # prysm
        # ------
        # prysm set up
        samples = 1800
        fno = 2.8  # f number
        efl = 100  # effective focal length in mm
        epd = efl / fno  # aperture diameter
        r_aper = epd / 2  # aperture radius

        xi, eta = make_xy_grid(samples, diameter=epd)  # type: ignore[arg-type]
        r, t = cart_to_polar(xi, eta)

        r_aber = r / r_aper

        nms = [ansi_j_to_nm(j) for j in range(0, 13)]
        basis = np.asarray(list(zernike_nm_seq(nms, r_aber, t)))

        # set a "reproducible" random engine
        seed = 1
        rnd_gen = np.random.default_rng(seed)

        # generate WFE coeffs (in nm)
        phs_coefs = rnd_gen.uniform(size=len(basis)) * 2000

        phs = sum_of_2d_modes(basis, phs_coefs)

        # opticks
        # -------

        # compute grid
        ap_diam = efl * u.mm / fno

        grid = Grid.from_size(samples, ap_diam)

        # compute opd
        wfe_rms = phs_coefs * u.nm
        opd = OptPathDiff.from_zernike(wfe_rms, ap_diam, grid)

        # verification
        np.testing.assert_allclose(opd.data, phs * u.nm, rtol=1e-12)

    # ---- Zernike ordering support ----

    samples = 256
    ap_diam: u.Quantity = 35 * u.mm

    @pytest.fixture(scope="class")
    def grid(self) -> Grid:
        return Grid.from_size(self.samples, self.ap_diam)

    @staticmethod
    def _reference_opd(nms, coefs, samples, ap_diam):
        """Build the reference OPD from raw prysm calls."""
        epd = ap_diam.to_value(u.mm)
        xi, eta = make_xy_grid(samples, diameter=epd)  # type: ignore[arg-type]
        r, t = cart_to_polar(xi, eta)
        r_aber = r / (epd / 2)
        basis = np.asarray(list(zernike_nm_seq(nms, r_aber, t)))
        return sum_of_2d_modes(basis, coefs)

    def test_zernike_ordering_default_is_ansi(self, grid):
        """Default ordering is ANSI — equivalent to explicit 'ansi'."""
        coefs = np.linspace(10.0, 120.0, 12)
        wfe_rms = coefs * u.nm

        opd_default = OptPathDiff.from_zernike(wfe_rms, self.ap_diam, grid)
        opd_ansi = OptPathDiff.from_zernike(
            wfe_rms, self.ap_diam, grid, ordering="ansi"
        )

        np.testing.assert_allclose(opd_default.data, opd_ansi.data, rtol=1e-12)

    def test_zernike_ordering_noll(self, grid):
        """Noll ordering uses 1-based noll_to_nm mapping."""
        coefs = np.linspace(5.0, 90.0, 10)
        wfe_rms = coefs * u.nm

        # reference: Noll is 1-based
        nms = [noll_to_nm(j) for j in range(1, len(coefs) + 1)]
        ref = self._reference_opd(nms, coefs, self.samples, self.ap_diam)

        opd = OptPathDiff.from_zernike(wfe_rms, self.ap_diam, grid, ordering="noll")

        np.testing.assert_allclose(opd.data, ref * u.nm, rtol=1e-12)

    def test_zernike_ordering_fringe(self, grid):
        """Fringe ordering uses 1-based fringe_to_nm mapping."""
        coefs = np.linspace(2.0, 80.0, 9)
        wfe_rms = coefs * u.nm

        # reference: Fringe is 1-based
        nms = [fringe_to_nm(j) for j in range(1, len(coefs) + 1)]
        ref = self._reference_opd(nms, coefs, self.samples, self.ap_diam)

        opd = OptPathDiff.from_zernike(wfe_rms, self.ap_diam, grid, ordering="fringe")

        np.testing.assert_allclose(opd.data, ref * u.nm, rtol=1e-12)

    def test_zernike_ordering_case_insensitive(self, grid):
        """The ordering argument is case-insensitive."""
        coefs = np.linspace(1.0, 50.0, 6)
        wfe_rms = coefs * u.nm

        opd_lower = OptPathDiff.from_zernike(
            wfe_rms, self.ap_diam, grid, ordering="fringe"
        )
        opd_upper = OptPathDiff.from_zernike(
            wfe_rms, self.ap_diam, grid, ordering="FRINGE"
        )

        np.testing.assert_allclose(opd_lower.data, opd_upper.data, rtol=1e-12)

    def test_zernike_ordering_invalid_raises(self, grid):
        """Unknown ordering names raise ValueError."""
        wfe_rms = np.ones(4) * u.nm

        with pytest.raises(ValueError, match="Unknown Zernike ordering"):
            OptPathDiff.from_zernike(wfe_rms, self.ap_diam, grid, ordering="zemax")
