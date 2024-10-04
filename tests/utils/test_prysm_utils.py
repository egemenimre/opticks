# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import numpy as np
from prysm.coordinates import cart_to_polar, make_xy_grid
from prysm.polynomials import ansi_j_to_nm, sum_of_2d_modes, zernike_nm_sequence

from opticks import u
from opticks.utils.prysm_utils import Grid, OptPathDiff
from opticks.utils.testing_utils import assert_allclose


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

        xi, eta = make_xy_grid(samples, diameter=epd)
        r, t = cart_to_polar(xi, eta)

        r_aber = r / r_aper

        nms = [ansi_j_to_nm(j) for j in range(0, 13)]
        basis = list(zernike_nm_sequence(nms, r_aber, t))

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
        assert_allclose(opd.data, phs * u.nm, rtol=1e-12)
