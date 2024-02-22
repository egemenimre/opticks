# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
from pint import Quantity
from strictyaml import Map, Str, Enum, Int, Optional

from opticks import u
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.yaml_helpers import Qty


detector_schema = {
    "name": Str(),
    "detector_type": Enum(["pushbroom", "full frame"]),
    "pixel_pitch": Qty(),
    "horizontal_pixels": Int(),
    "vertical_pixels": Int(),
    "horizontal_pixels_used": Int(),
    "vertical_pixels_used": Int(),
    Optional("binning", default=1): Int(),
    # TDI ignored for full frame detectors
    Optional("tdi_stages", default=1): Int(),
    Optional("full_well_capacity"): Qty(),
    Optional("noise"): Map(
        {
            Optional("dark_current"): Qty(),
            Optional("temporal_dark_noise"): Qty(),
        }
    ),
    "timings": Map(
        {
            # TODO revalid: required for frame-rate imager only
            Optional("frame_rate", default=None): Qty(),
            "integration_duration": Qty(),
            Optional("frame_overhead_duration", default=0 * u.ms): Qty(),
            Optional("frame_overlap_duration", default=0 * u.ms): Qty(),
        }
    ),
}
"""Schema containing detector parameters."""

# TODO schema: what to do with zero Optionals and None Optionals


class Detector(ImagerComponent):
    """
    Class containing generic Detector parameters.
    """

    @classmethod
    def schema(cls) -> Map:
        return Map(detector_schema)

    def _params_class_name(cls) -> str:
        return "DetectorParams"

    # ---------- begin modelling functions ----------

    def pix_pitch(self, with_binning: bool = True) -> Quantity:
        """
        Returns pixel pitch parameter with or without binning.

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not

        Returns
        -------
        Quantity
            Pixel pitch with or without binning
        """
        binning = self.params.binning if with_binning else 1

        return self.params.pixel_pitch * binning

    def nyquist_freq(self, with_binning: bool = True) -> Quantity:
        """
        Returns Nyquist Frequency / Limit parameter with or without binning.

        Nyquist frequency is defined as: $\frac{1}{2 \times \text{pix pitch}}$,
        which translates to one line pair per two pixels (e.g., one pixel white,
        next one black). This is the absolute maximum limiting resolution of the sensor.

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not

        Returns
        -------
        Quantity
            Nyquist limit in lp/mm
        """

        u.define("lp = 1 * dimensionless = lp")

        return 1 * u.lp / (2 * self.pix_pitch(with_binning)).to(u.mm)

    def pixel_count(self, with_binning: bool = True, used: bool = True) -> Quantity:
        """
        Computes the total number of pixels in the detector with/without binning
        or total pixels/used pixels.

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not
        used : bool
            Return the value for used or total number of pixels
        Returns
        -------
        Quantity
            Total number of pixels for the requested configuration
        """

        binning = self.params.binning if with_binning else 1

        if used:
            # number of used pixels in the frame
            return (
                self.params.horizontal_pixels_used
                / binning
                * self.params.vertical_pixels_used
                / binning
                * u.pixel
            ).to("Mpixel")
        else:
            # total number of pixels in the frame
            return (
                self.params.horizontal_pixels
                / binning
                * self.params.vertical_pixels
                / binning
                * u.pixel
            ).to("Mpixel")

    def pix_area(self, with_binning: bool = True) -> Quantity:
        """
        Pixel physical area.

        Computed as $text{pix pitch} \times \text{pix pitch}$

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not

        Returns
        -------
        Quantity
            Pixel area with or without binning
        """

        return self.pix_pitch(with_binning) ** 2
