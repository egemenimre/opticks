# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
from pint import Quantity
from strictyaml import YAML, Enum, Int, Map, Optional, Str

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

    def __init__(self, yaml: YAML):
        super().__init__(yaml)

        # init some useful parameters
        self._init_useful_params()

    def _init_useful_params(self):
        # init some useful parameters

        # frame duration: inverse of the frame rate
        self.params.timings["frame_duration"] = (1 / self.params.timings.frame_rate).to(
            u.ms
        )
        self.params["is_binned"] = False if self.params.binning == 1 else True
        self.params.timings["max_integration_duration"] = _max_integration_duration(
            self.params.timings
        )
        self.params.timings["total_tdi_col_duration"] = _total_tdi_col_duration(
            self.params
        )

    @classmethod
    def schema(cls) -> Map:
        return Map(detector_schema)

    @classmethod
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

    def pix_area(self, with_binning: bool = True) -> Quantity:
        """
        Pixel physical area.

        Computed as the square of pixel pitch. The effective pixel area is smaller,
        even for a pixel with a single sensing element.

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

    def nyquist_freq(self, with_binning: bool = True) -> Quantity:
        r"""
        Returns Nyquist Frequency / Limit parameter with or without binning.

        Nyquist frequency is defined as: `1 / (2 x pix pitch)`,
        which translates to one line pair per two pixels (e.g., one pixel white,
        next one black). This is the absolute maximum limiting resolution of the sensor.

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not

        Returns
        -------
        nyquist_limit: Quantity
            Nyquist limit in lp/mm
        """

        u.define("lp = 1 * dimensionless = lp")

        return 1 * u.lp / (2 * self.pix_pitch(with_binning)).to(u.mm)

    def pixel_count_frame(
        self, with_binning: bool = True, used: bool = True
    ) -> Quantity:
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

    def pixel_count_line(
        self, with_binning: bool = True, used: bool = True
    ) -> Quantity:
        """
        Computes the total number of pixels in the detector line with/without binning
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
            return (self.params.horizontal_pixels_used / binning * u.pixel).to("Mpixel")
        else:
            # total number of pixels in the frame
            return (self.params.horizontal_pixels / binning * u.pixel).to("Mpixel")

    def pix_read_rate(
        self, with_binning: bool = True, with_tdi: bool = True
    ) -> Quantity:
        r"""
        Pixel read rate.

        Computed as:
         - Pushbroom type: `horiz pix (binned) x TDI stages x line rate`
         - Full frame type: `horiz pix (binned) x vert pix (binned) x frame rate`


        Note that the unused pixels are also read, this assumes that the
        detector does not have ROI functionality.

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not
        with_tdi : bool
            Return the value with TDI or not (valid for pushbroom only)

        Returns
        -------
        Quantity
            Pixel read rate with or without binning (Mpixel/s)
        """
        tdi = self.params.tdi_stages if with_tdi else 1

        if self.params.detector_type == "pushbroom":
            pix_read_rate = (
                self.pixel_count_line(with_binning, False)
                * tdi
                / self.params.timings.frame_duration
            )
        elif self.params.detector_type == "full frame":
            pix_read_rate = (
                self.pixel_count_frame(with_binning, False)
                / self.params.timings.frame_duration
            )
        else:
            raise ValueError(f"Invalid detector type: {self.params.detector_type}")

        return pix_read_rate.to("Mpixel/s")


def _max_integration_duration(timings) -> Quantity:
    """
    Computes the maximum integration duration, where:

    `max integ duration = line or frame duration (w/o bin) - frame overhead duration
    + frame overlap duration`
    """
    return (
        timings.frame_duration
        - timings.frame_overhead_duration
        + timings.frame_overlap_duration
    )


def _total_tdi_col_duration(params) -> Quantity:
    """
    Computes the total TDI column duration (if applicable).

    `total tdi col duration = line or frame duration (w/o bin) x Number of TDI stages`
    """
    return (params.timings.frame_duration * params.tdi_stages).to(u.ms)
