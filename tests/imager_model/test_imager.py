# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import pytest
from astropy.units import Unit
from numpy.testing import assert_allclose

from opticks import u
from opticks.imager_model.imager import Imager
from tests import process_paths


class TestImager:

    file_directory = Path("sat_pushbroom_data")
    alt_file_directory = Path("tests", "imager_model", "sat_pushbroom_data")
    optics_file_path = Path("optics.yaml")
    detector_file_path = Path("pan_detector.yaml")
    rw_electronics_file_path = Path("rw_electronics.yaml")

    band_id = "pan"

    @pytest.fixture(scope="class")
    def imager(self) -> Imager:
        """Loads the yaml file and inits the Imager."""

        # different test environments work with different paths
        optics_file_path = process_paths(
            self.optics_file_path, self.file_directory, self.alt_file_directory
        )
        detector_file_path = process_paths(
            self.detector_file_path, self.file_directory, self.alt_file_directory
        )
        rw_electronics_file_path = process_paths(
            self.rw_electronics_file_path, self.file_directory, self.alt_file_directory
        )

        # checks only whether the file can be opened and whether it fits the schema
        return Imager.from_yaml_file(
            optics_file_path,
            detector_file_path,
            rw_electronics_file_path,
        )

    def test_q_factor(self, imager: Imager):
        """Tests the Q factor."""

        # check values
        truth = 0.7636094674556213

        # select the PAN channel
        ref_wavelength = 500 * u.nm

        # Generate the Q Factor
        q = imager.q_factor(ref_wavelength, self.band_id, with_binning=True)

        # verification
        assert q == pytest.approx(truth, abs=1e-9)

    def test_fov(self, imager: Imager):
        """Tests the horiz and vert FoV."""

        # set up
        horiz_fov_truth = 1.73139509 * u.deg
        vert_fov_truth = 0.0577131695 * u.mdeg

        # computation
        horiz_fov = imager.horizontal_fov(self.band_id)
        vert_fov = imager.vertical_fov(self.band_id)

        # verification
        assert_allclose(horiz_fov, horiz_fov_truth, atol=0.001 * u.mdeg)
        assert_allclose(vert_fov, vert_fov_truth, atol=0.00001 * u.mdeg)

    def test_ifov(self, imager: Imager):
        """Tests the instantaneous FoV."""

        # set up
        truth = 0.0577131695 * u.mdeg

        # computation
        ifov = imager.ifov(self.band_id, with_binning=False)

        # verification
        assert_allclose(ifov, truth, atol=0.00000001 * u.mdeg)

    def test_ssd_centre(self, imager: Imager):
        """Tests the Spatial Sample Distance."""

        # set up - centre
        truth_horiz = 0.699055672 * u.m
        truth_vert = 0.699055672 * u.m

        sat_altitude = 694.0 * u.km

        # computation
        ssd = imager.spatial_sample_distance(
            sat_altitude, self.band_id, False, "centre"
        )

        # verification
        assert_allclose(ssd.horiz, truth_horiz, atol=0.00000001 * u.m)
        assert_allclose(ssd.vert, truth_vert, atol=0.00000001 * u.m)

        # set up - centre right
        truth_horiz = 0.699215273 * u.m
        truth_vert = 0.699135473 * u.m

        # computation
        ssd = imager.spatial_sample_distance(
            sat_altitude, self.band_id, False, "centre right"
        )

        # verification
        assert_allclose(ssd.horiz, truth_horiz, atol=0.00000001 * u.m)
        assert_allclose(ssd.vert, truth_vert, atol=0.00000001 * u.m)

    def test_rw(self, imager: Imager):
        """Tests the read/write rates."""

        # set up
        truth_pix_read_rate_no_tdi = 294.348 * Unit("megapixel / second")
        truth_pix_read_rate_tdi = 2943.48 * Unit("megapixel / second")
        truth_data_write_rate_uncomp = 3638.14128 * Unit("megabit / second")
        truth_data_write_rate_comp = 909.53532 * Unit("megabit / second")

        # computation

        #  pixel read rate (without TDI)
        pix_read_rate_no_tdi = imager.detector.pix_read_rate(self.band_id, False, False)

        # pixel read rate (with TDI)
        pix_read_rate_tdi = imager.detector.pix_read_rate(self.band_id, False, True)

        # data write rate (uncompressed, incl. overheads)
        data_write_rate_uncomp = imager.data_write_rate(self.band_id, False, False)

        # data write rate (compressed, incl. overheads)
        data_write_rate_comp = imager.data_write_rate(self.band_id, False, True)

        # verification
        assert_allclose(
            pix_read_rate_no_tdi,
            truth_pix_read_rate_no_tdi,
            atol=0.00000001 * Unit("megapixel / second"),
        )
        assert_allclose(
            pix_read_rate_tdi,
            truth_pix_read_rate_tdi,
            atol=0.00000001 * Unit("megapixel / second"),
        )
        assert_allclose(
            data_write_rate_uncomp,
            truth_data_write_rate_uncomp,
            atol=0.00000001 * Unit("megabit / second"),
        )
        assert_allclose(
            data_write_rate_comp,
            truth_data_write_rate_comp,
            atol=0.00000001 * Unit("megabit / second"),
        )
