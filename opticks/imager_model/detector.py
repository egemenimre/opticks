# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from collections.abc import Iterable
from enum import StrEnum

from astropy.units import Quantity
from pydantic import BaseModel, ConfigDict, Field, PositiveInt, PrivateAttr

from opticks import u
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.parser_helpers import PositivePydanticQty, PydanticQty


class DetectorType(StrEnum):
    """Detector type enumeration."""

    PUSHBROOM = "pushbroom"
    FULL_FRAME = "full frame"


class Channel(BaseModel):

    model_config = ConfigDict(arbitrary_types_allowed=True)

    name: str
    horizontal_pixels: PositiveInt
    vertical_pixels: PositiveInt
    binning: int = 1
    tdi_stages: int = 1
    cuton_wvl: PositivePydanticQty
    cutoff_wvl: PositivePydanticQty

    # Derived attributes (set by Detector in model_post_init, not serialized)
    _detector_type: DetectorType | None = PrivateAttr(default=None)
    _det_pixel_pitch: Quantity | None = PrivateAttr(default=None)
    is_binned: bool = Field(default=False, exclude=True)
    frame_duration: Quantity | None = Field(default=None, exclude=True)
    frame_rate: Quantity | None = Field(default=None, exclude=True)
    integration_duration: Quantity | None = Field(default=None, exclude=True)
    total_tdi_col_duration: Quantity | None = Field(default=None, exclude=True)

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

        return self._det_pixel_pitch * binning  # type: ignore[operator]

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
        ).to(
            "Mpixel"
        )  # type: ignore[return-value]

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
        return (self.horizontal_pixels / binning * u.pixel).to("Mpixel")  # type: ignore[return-value]

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

        if self._detector_type == DetectorType.PUSHBROOM:
            pix_read_rate = (
                self.pixel_count_line(with_binning) * tdi * (frame_rate / binning)
            )
        elif self._detector_type == DetectorType.FULL_FRAME:
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

        return (self.cuton_wvl + self.cutoff_wvl) * 0.5  # type: ignore[return-value]

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

        return abs(self.cutoff_wvl - self.cuton_wvl)  # type: ignore[return-value]


class Timings(BaseModel):
    """Timing parameters for a Detector."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    frame_rate: PydanticQty | None = None
    integration_duration: PydanticQty
    frame_overhead_duration: PydanticQty = Field(
        default_factory=lambda: Quantity(0, u.ms)
    )
    frame_overlap_duration: PydanticQty = Field(
        default_factory=lambda: Quantity(0, u.ms)
    )

    # Derived (set by Detector in model_post_init)
    frame_duration: Quantity | None = Field(default=None, exclude=True)
    max_integration_duration: Quantity | None = Field(default=None, exclude=True)


class Noise(BaseModel):
    """Noise parameters for a Detector."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    dark_current: PydanticQty | None = None
    temporal_dark_noise: PydanticQty | None = None


class Detector(ImagerComponent):
    """
    Class containing generic Detector parameters.

    A `Detector` has one or more `Channel` objects.
    """

    name: str
    detector_type: DetectorType
    pixel_pitch: PositivePydanticQty
    horizontal_pixels: PositiveInt
    vertical_pixels: PositiveInt
    channels: dict[str, Channel]
    full_well_capacity: PydanticQty | None = None
    noise: Noise | None = None
    timings: Timings

    def model_post_init(self, __context):
        """Initialises derived parameters after Pydantic validation.

        Overrides the abstract pydantic method to perform additional initialization
        after `__init__` and `model_construct`."""
        self._init_useful_params()

    def _init_useful_params(self):
        """
        Initialises some useful parameters.
        """

        # shorthand for timings
        timings = self.timings

        # frame duration: inverse of the frame rate
        timings.frame_duration = (1 / timings.frame_rate).to(u.ms)  # type: ignore[operator]
        # max integration duration possible
        timings.max_integration_duration = (  # type: ignore[misc]
            timings.frame_duration
            - timings.frame_overhead_duration
            + timings.frame_overlap_duration
        )

        # per channel params
        for channel in self.channels.values():
            # physical pixel pitch size
            channel._detector_type = self.detector_type
            # physical pixel pitch size
            channel._det_pixel_pitch = self.pixel_pitch
            # is binned
            channel.is_binned = False if channel.binning == 1 else True
            # line/frame duration (with binning)
            channel.frame_duration = timings.frame_duration * channel.binning  # type: ignore[operator]
            # line/frame rate (with binning)
            channel.frame_rate = timings.frame_rate / channel.binning  # type: ignore[operator]
            # total exposure duration (with binning)
            channel.integration_duration = (
                timings.integration_duration * channel.binning  # type: ignore[operator]
            )
            # total TDI column duration
            # total tdi col duration = line/frame duration (w/bin) x Nr of TDI stages
            channel.total_tdi_col_duration = (
                timings.frame_duration * channel.binning * channel.tdi_stages  # type: ignore[operator]
            ).to(u.ms)

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
        return (self.horizontal_pixels * self.vertical_pixels * u.pixel).to("Mpixel")  # type: ignore[return-value]

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
        timings = self.timings

        if isinstance(band_id, str):

            # there is a single channel
            channel = self.get_channel(band_id)

            pix_read_rate = channel.pix_read_rate(
                timings.frame_rate, with_binning, with_tdi  # type: ignore[arg-type]
            )
        else:
            # there are multiple channels, sum the values
            channels = self.get_channels(band_id)

            for channel in channels:
                pix_read_rate_single_chan = channel.pix_read_rate(
                    timings.frame_rate, with_binning, with_tdi  # type: ignore[arg-type]
                )
                pix_read_rate += pix_read_rate_single_chan

        return pix_read_rate.to("Mpixel/s")  # type: ignore[return-value]

    def get_channel(self, band_id: str) -> Channel:
        """
        Gets the channel with the 'band_id'.

        Parameters
        ----------
        band_id : str
            band ID corresponding to the requested channel

        Returns
        -------
        Channel
            Requested channel with the 'band_id'
        """
        return self.channels[band_id]

    def get_channels(self, band_ids: Iterable[str]) -> list[Channel]:
        """
        Gets the list of channels with the 'band_id'.

        Parameters
        ----------
        band_ids : Iterable[str]
            band IDs corresponding to the requested channels

        Returns
        -------
        list[Channel]
            Requested channels with the 'band_id'
        """

        return [self.channels[band_id] for band_id in band_ids]
