# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
from collections.abc import Iterable

import numpy as np
from imager_model.optics import Optics
from pint import Quantity
from strictyaml import YAML, Enum, Int, Map, MapPattern, Optional, Str

from opticks import u
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.yaml_helpers import Qty

channel_schema = {
    "name": Str(),
    "horizontal_pixels": Int(),
    "vertical_pixels": Int(),
    Optional("binning", default=1): Int(),
    # TDI ignored for full frame detectors
    Optional("tdi_stages", default=1): Int(),
}
"""Schema containing per-channel parameters."""

detector_schema = {
    "name": Str(),
    "detector_type": Enum(["pushbroom", "full frame"]),
    "pixel_pitch": Qty(),
    "horizontal_pixels": Int(),
    "vertical_pixels": Int(),
    "channels": MapPattern(Str(), Map(channel_schema)),
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


class Channel:

    def __init__(self):
        self.tdi_stages = None
        self._detector_type = None
        self.vertical_pixels = None
        self.horizontal_pixels = None
        self.binning = None
        self._det_pixel_pitch = None

    def pixel_pitch(self, with_binning: bool = True) -> Quantity:
        """
        Returns pixel pitch parameter with or without binning.

        This corresponds to an effective pixel pitch.

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not

        Returns
        -------
        Quantity
            Pixel pitch with or without binning
        """
        binning = self.binning if with_binning else 1

        return self._det_pixel_pitch * binning

    def pixel_area(self, with_binning: bool = True) -> Quantity:
        """
        Pixel geometric area.

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

        return self.pixel_pitch(with_binning) ** 2

    def ifov(self, optics: Optics, with_binning: bool = True) -> Quantity:
        """
        Computes the Instantaneous field of view (works in vertical and horizontal).

        Assumes constant IFOV per pixel.

        Parameters
        ----------
        optics : Optics
            Optics in front of the detector
        with_binning : bool
            Return the value with binning or not

        Returns
        -------
        Quantity
            IFOV angle

        """
        return 2 * np.arctan(
            (self.pixel_pitch(with_binning) / 2.0) / optics.params.focal_length
        ).to(u.mdeg)

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

        return 1 * u.lp / (2 * self.pixel_pitch(with_binning).to(u.mm))

    def pixel_count_frame(self, with_binning: bool = True) -> Quantity:
        """
        Computes the total number of pixels in the channel with/without binning.

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not
        Returns
        -------
        Quantity
            Total number of pixels for the requested configuration
        """

        binning = self.binning if with_binning else 1

        # total number of pixels in the frame
        return (
            self.horizontal_pixels / binning * self.vertical_pixels / binning * u.pixel
        ).to("Mpixel")

    def pixel_count_line(self, with_binning: bool = True) -> Quantity:
        """
        Computes the total number of pixels in the detector line with/without binning.

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not
        Returns
        -------
        Quantity
            Total number of pixels for the requested configuration
        """

        binning = self.binning if with_binning else 1

        # total number of pixels in the frame
        return (self.horizontal_pixels / binning * u.pixel).to("Mpixel")

    def pix_solid_angle(self, optics: Optics, with_binning=True) -> Quantity:
        """
        Pixel solid angle (of a pyramid).

        Parameters
        ----------
        optics : Optics
            Optics in front of the detector
        with_binning : bool
            Return the value with binning or not

        Returns
        -------
        Quantity
            Pixel solid angle in steradians
        """
        #

        pix_solid_angle = 4 * np.arcsin(
            np.sin(self.ifov(optics, with_binning) / 2.0)
            * np.sin(self.ifov(optics, with_binning) / 2.0)
        )

        # correct the unit from rad to sr (or rad**2)
        return (pix_solid_angle * u.rad).to(u.steradian)

    def horizontal_fov(self, optics: Optics) -> Quantity:
        """
        Computes the full field of view in the horizontal direction.

        Assumes constant IFOV per pixel. Used pixels only.

        Parameters
        ----------
        optics : Optics
            Optics in front of the detector

        Returns
        -------
        Quantity
            Horizontal FOV angle

        """
        return 2 * np.tan(self.ifov(optics, False) * self.horizontal_pixels / 2.0).to(
            u.deg
        )

    def vertical_fov(self, optics: Optics) -> Quantity:
        """
        Computes the full field of view in the vertical direction.

        Assumes constant IFOV per pixel. Used pixels only.

        Parameters
        ----------
        optics : Optics
            Optics in front of the detector

        Returns
        -------
        Quantity
            Vertical FOV angle

        """
        return 2 * np.tan(self.ifov(optics, False) * self.vertical_pixels / 2.0).to(
            u.deg
        )

    def pix_read_rate(
        self, frame_rate: Quantity, with_binning: bool = True, with_tdi: bool = True
    ) -> Quantity:
        r"""
        Pixel read rate.

        Computed as:
         - Pushbroom type: `horiz pix (binned) x TDI stages x line rate`
         - Full frame type: `horiz pix (binned) x vert pix (binned) x frame rate

        Parameters
        ----------
        frame_rate : Quantity
            Frame or line rate (in Hz)
        with_binning : bool
            Return the value with binning or not
        with_tdi : bool
            Return the value with TDI or not (valid for pushbroom only)

        Returns
        -------
        Quantity
            Pixel read rate with or without binning (Mpixel/s)
        """
        tdi = self.tdi_stages if with_tdi else 1

        if self._detector_type == "pushbroom":
            pix_read_rate = self.pixel_count_line(with_binning) * tdi * frame_rate
        elif self._detector_type == "full frame":
            pix_read_rate = self.pixel_count_frame(with_binning) * frame_rate
        else:
            raise ValueError(f"Invalid detector type: {self._detector_type}")

        return pix_read_rate.to("Mpixel/s")


class Detector(ImagerComponent):
    """
    Class containing generic Detector parameters.
    """

    def __init__(self, yaml: YAML):
        super().__init__(yaml)

        # modify some internal auto generated classes
        self._prepare_internal_classes()

        # init some useful parameters
        self._init_useful_params()

    def _prepare_internal_classes(self):
        # make the internal dict of the channels accessible as a dict
        self.params.channels.all = {
            key: value for key, value in self.params.channels.__dict__.items()
        }

        # cast all channels into a common Channel object
        for channel in self.params.channels.all.values():
            channel.__class__ = Channel

    def _init_useful_params(self):
        # init some useful parameters

        # shorthand for timings
        timings = self.params.timings

        # frame duration: inverse of the frame rate
        timings["frame_duration"] = (1 / timings.frame_rate).to(u.ms)
        # max integration duration possible
        self.params.timings["max_integration_duration"] = (
            timings.frame_duration
            - timings.frame_overhead_duration
            + timings.frame_overlap_duration
        )

        # per channel params
        for channel in self.params.channels.all.values():
            # physical pixel pitch size
            channel._detector_type = self.params.detector_type
            # physical pixel pitch size
            channel._det_pixel_pitch = self.params.pixel_pitch
            # is binned
            channel.is_binned = False if channel.binning == 1 else True
            # total TDI column duration
            # total tdi col duration = line/frame duration (w/o bin) x Nr of TDI stages
            channel.total_tdi_col_duration = (
                timings.frame_duration * channel.tdi_stages
            ).to(u.ms)

    @classmethod
    def schema(cls) -> Map:
        return Map(detector_schema)

    @classmethod
    def _params_class_name(cls) -> str:
        return "DetectorParams"

    # ---------- begin modelling functions ----------

    @property
    def pixel_count(self) -> Quantity:
        """
        Computes the total number of pixels in the detector.

        Parameters
        ----------

        Returns
        -------
        Quantity
            Total number of pixels for the requested configuration
        """

        # total number of pixels in the frame
        return (
            self.params.horizontal_pixels * self.params.vertical_pixels * u.pixel
        ).to("Mpixel")

    def pix_read_rate(
        self,
        channel: Channel | Iterable[Channel],
        with_binning: bool = True,
        with_tdi: bool = True,
    ) -> Quantity:
        r"""
        Pixel read rate.

        Computed as:
         - Pushbroom type: `horiz pix (binned) x TDI stages x line rate`
         - Full frame type: `horiz pix (binned) x vert pix (binned) x frame rate`

        If a list of channels is given, then  the returned result is the sum of all
        the requested channels.

        Parameters
        ----------
        channel : Channel or Iterable[Channel]
            Channels to compute the readout
        with_binning : bool
            Return the value with binning or not
        with_tdi : bool
            Return the value with TDI or not (valid for pushbroom only)

        Returns
        -------
        Quantity
            Pixel read rate with or without binning (Mpixel/s)
        """

        pix_read_rate = 0 * u.Mpixel / u.s

        if isinstance(channel, Iterable):
            # there are multiple channels, sum the values
            for single_channel in channel:
                pix_read_rate_single_chan = single_channel.pix_read_rate(
                    self.params.timings.frame_rate, with_binning, with_tdi
                )
                pix_read_rate += pix_read_rate_single_chan
        else:
            # there is a single channel
            pix_read_rate = channel.pix_read_rate(
                self.params.timings.frame_rate, with_binning, with_tdi
            )

        return pix_read_rate.to("Mpixel/s")
