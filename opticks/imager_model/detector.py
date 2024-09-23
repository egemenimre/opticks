# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
from collections.abc import Iterable
from typing import List

import numpy as np
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
    "cuton_wvl": Qty(),
    "cutoff_wvl": Qty(),
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
        self.tdi_stages: int = None
        self._detector_type = None
        self.vertical_pixels: int = None
        self.horizontal_pixels: int = None
        self.binning: int = None
        self._det_pixel_pitch: float = None

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

    def nyquist_freq(self, with_binning: bool = True) -> Quantity:
        """
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
            Nyquist limit in cycles/mm
        """

        return 1 * u.cy / (2 * self.pixel_pitch(with_binning).to(u.mm))

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
            Total number of pixels for the requested configuration (in Mpixels)
        """

        binning = self.binning if with_binning else 1

        # total number of pixels in the frame
        return (self.horizontal_pixels / binning * u.pixel).to("Mpixel")

    def pix_read_rate(
        self, frame_rate: Quantity, with_binning: bool = True, with_tdi: bool = False
    ) -> Quantity:
        r"""
        Pixel read rate.

        Computed as:
            - Pushbroom type: `horiz pix (binned) x TDI stages x line rate`
            - Full frame type: `horiz pix (binned) x vert pix (binned) x frame rate`

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

        binning = self.binning if with_binning else 1

        if self._detector_type == "pushbroom":
            pix_read_rate = (
                self.pixel_count_line(with_binning) * tdi * (frame_rate / binning)
            )
        elif self._detector_type == "full frame":
            pix_read_rate = self.pixel_count_frame(with_binning) * (
                frame_rate / binning
            )
        else:
            raise ValueError(f"Invalid detector type: {self._detector_type}")

        return pix_read_rate.to("Mpixel/s")

    @property
    def centre_wavelength(self) -> Quantity:
        """
        Computes the centre wavelength of the channel.

        Returns
        -------
        Quantity
            Centre wavelength of the channel
        """

        return (self.cuton_wvl + self.cutoff_wvl) * 0.5

    @property
    def bandwidth(self) -> Quantity:
        """
        Computes the bandwidth of the channel
        (the difference between the cut-off and cut-on wavelengths).

        Returns
        -------
        Quantity
            Bandwidth of the channel (in wavelengths)
        """

        return np.abs(self.cutoff_wvl - self.cuton_wvl)


class Detector(ImagerComponent):
    """
    Class containing generic Detector parameters.

    A 'Detector' has one or more 'Channel's.
    """

    def __init__(self, yaml: YAML):
        super().__init__(yaml)

        # modify some internal auto generated classes
        self._prepare_internal_classes()

        # init some useful parameters
        self._init_useful_params()

    def _prepare_internal_classes(self):
        """
        Prepares the internal classes upon initialisation.
        """
        # make the internal dict of the channels accessible as a dict
        self.params.channels.all = {
            key: value for key, value in self.params.channels.__dict__.items()
        }

        # cast all channels into a common Channel object
        for channel in self.params.channels.all.values():
            channel.__class__ = Channel

    def _init_useful_params(self):
        """
        Initialises some useful parameters.
        """

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
            # line/frame duration (with binning)
            channel.frame_duration = timings.frame_duration * channel.binning
            # line/frame rate (with binning)
            channel.frame_rate = timings.frame_rate / channel.binning
            # total exposure duration (with binning)
            channel.integration_duration = (
                timings.integration_duration * channel.binning
            )
            # total TDI column duration
            # total tdi col duration = line/frame duration (w/bin) x Nr of TDI stages
            channel.total_tdi_col_duration = (
                timings.frame_duration * channel.binning * channel.tdi_stages
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
            Total number of pixels for the requested configuration (in Mpixels)
        """

        # total number of pixels in the frame
        return (
            self.params.horizontal_pixels * self.params.vertical_pixels * u.pixel
        ).to("Mpixel")

    def pix_read_rate(
        self,
        band_id: str | Iterable[str],
        with_binning: bool = True,
        with_tdi: bool = True,
    ) -> Quantity:
        """
        Pixel read rate.

        Computed as:
         - Pushbroom type: `horiz pix (binned) x TDI stages x line rate`
         - Full frame type: `horiz pix (binned) x vert pix (binned) x frame rate`

        If a list of channels is given, then  the returned result is the sum of all
        the requested channels.

        Parameters
        ----------
        band_id : str or Iterable[str]
            band ID to compute the readout
        with_binning : bool
            Return the value with binning or not
        with_tdi : bool
            Return the value with TDI or not (valid for pushbroom only)

        Returns
        -------
        Quantity
            Pixel read rate with or without binning (Mpixels)
        """

        pix_read_rate = 0 * u.Mpixel / u.s

        # shorthand
        timings = self.params.timings

        if isinstance(band_id, str):

            # there is a single channel
            channel = self.get_channel(band_id)

            pix_read_rate = channel.pix_read_rate(
                timings.frame_rate, with_binning, with_tdi
            )
        else:
            # there are multiple channels, sum the values
            channels = self.get_channels(band_id)

            for channel in channels:
                pix_read_rate_single_chan = channel.pix_read_rate(
                    timings.frame_rate, with_binning, with_tdi
                )
                pix_read_rate += pix_read_rate_single_chan

        return pix_read_rate.to("Mpixel/s")

    def get_channel(self, band_id: str) -> Channel:
        """
        Gets the channel with the 'band_id'.

        Alias for 'self.params.channels.all[band_id]'.

        Parameters
        ----------
        band_id : str
            band ID corresponding to the requested channel

        Returns
        -------
        Channel
            Requested channel with the 'band_id'
        """
        return self.params.channels.all[band_id]

    def get_channels(self, band_ids: Iterable[str]) -> List[Channel]:
        """
        Gets the list channel with the 'band_id'.

        Alias for 'self.params.channels.all[band_id]'.

        Parameters
        ----------
        band_ids : Iterable[str]
            band IDs corresponding to the requested channels

        Returns
        -------
        Channel
            Requested channels with the 'band_id'
        """

        return [self.params.channels.all[band_id] for band_id in band_ids]
