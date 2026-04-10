# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from collections.abc import Iterable
from enum import StrEnum
from pathlib import Path
from typing import Self

import numpy as np
from astropy.units import Quantity, Unit
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    PositiveInt,
    PrivateAttr,
    model_validator,
)

from opticks import u
from opticks.contrast_model.detector_mtf import (
    DetectorDiffusionModel,
    DetectorDiffusionPreset,
)
from opticks.contrast_model.mtf import MTF_Model_1D
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.math_utils import InterpolatorWithUnits, InterpolatorWithUnitTypes
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
    read_blocks: PositiveInt = 1
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

        return self.pixel_pitch(with_binning) ** 2  # type: ignore[return-value]

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
        ).to("Mpixel")  # type: ignore[return-value]

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
        self,
        frame_rate: Quantity,
        with_binning: bool = True,
        with_tdi: bool = False,
        with_read_blocks: bool = True,
    ) -> Quantity:
        r"""
        Pixel read rate.

        Computed as:
            - Pushbroom type: `horiz pix (binned) x TDI stages x read blocks x line rate`
            - Full frame type: `horiz pix (binned) x vert pix (binned) x frame rate`

        Parameters
        ----------
        frame_rate : Quantity
            Frame or line rate (in Hz)
        with_binning : bool
            Return the value with binning or not
        with_tdi : bool
            Return the value with TDI or not (valid for pushbroom only)
        with_read_blocks : bool
            Return the value with read blocks or not (valid for pushbroom only)

        Returns
        -------
        Quantity
            Pixel read rate with or without binning (Mpixel/s)
        """
        tdi = self.tdi_stages if with_tdi else 1

        binning = self.binning if with_binning else 1

        read_blocks = self.read_blocks if with_read_blocks else 1

        if self._detector_type == DetectorType.PUSHBROOM:
            pix_read_rate = (
                self.pixel_count_line(with_binning)
                * tdi
                * read_blocks
                * (frame_rate / binning)
            )
        elif self._detector_type == DetectorType.FULL_FRAME:
            pix_read_rate = self.pixel_count_frame(with_binning) * (
                frame_rate / binning
            )
        else:
            raise ValueError(f"Invalid detector type: {self._detector_type}")

        return pix_read_rate.to("Mpixel/s")  # type: ignore[union-attr]

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


class AbsorptionData(BaseModel):
    """Reference to an external absorption coefficient data file."""

    file: str
    csv_separator: str | None = None


class SensorParams(BaseModel):
    """Image sensor physical parameters (diffusion MTF, CTE, crosstalk)."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    # Diffusion model selection (mutually exclusive)
    diffusion_model: str | None = None
    diffusion_preset: str | None = None

    # Diffusion geometry params (required when using diffusion_model,
    # overrides when using preset)
    diffusion_length: PositivePydanticQty | None = None
    field_free_depth: PositivePydanticQty | None = None
    depletion_depth: PositivePydanticQty | None = None
    surface_recomb_velocity: PydanticQty | None = None
    diffusion_coeff: PydanticQty | None = None

    # Crosstalk coefficients (dimensionless fractions)
    crosstalk_xs: float | None = None
    crosstalk_xd: float | None = None

    # Absorption data file (relative to cwd)
    absorption_data: AbsorptionData | None = None

    # Loaded interpolator (private, not serialized)
    _absorption_ipol: InterpolatorWithUnits | None = PrivateAttr(default=None)

    @model_validator(mode="after")
    def _check_mutually_exclusive(self) -> Self:
        if self.diffusion_model is not None and self.diffusion_preset is not None:
            raise ValueError(
                "diffusion_model and diffusion_preset are mutually exclusive. "
                "Supply only one."
            )
        return self

    def load_absorption_data(self) -> None:
        """Load absorption coefficient data from the external file.

        Resolves the file path relative to the current working directory.
        The file is self-describing: comment lines start with ``#``,
        the first non-comment line is the header declaring column names
        and units (e.g. ``wavelength (nm)  alpha (1/um)``),
        and subsequent lines are data rows.
        """
        if self.absorption_data is None:
            raise ValueError("No absorption_data configured in sensor_params.")

        file_path = Path.cwd() / self.absorption_data.file
        sep = self.absorption_data.csv_separator

        wavelengths: list[float] = []
        alphas: list[float] = []
        x_unit: Unit | None = None
        y_unit: Unit | None = None

        with open(file_path) as f:
            for line in f:
                stripped = line.strip()
                # skip empty lines and comments
                if not stripped or stripped.startswith("#"):
                    continue
                if x_unit is None:
                    # first non-comment line is the header
                    x_unit, y_unit = _parse_header_units(stripped, sep)
                else:
                    # data row
                    parts = stripped.split(sep) if sep else stripped.split()
                    wavelengths.append(float(parts[0].strip()))
                    alphas.append(float(parts[1].strip()))

        x_qty = np.array(wavelengths) * x_unit
        y_qty = np.array(alphas) * y_unit

        self._absorption_ipol = InterpolatorWithUnits.from_ipol_method(
            InterpolatorWithUnitTypes.PCHIP, x_qty, y_qty
        )

    def get_absorption_coeff(self, wavelength: Quantity) -> Quantity:
        """Return the absorption coefficient at the given wavelength.

        Parameters
        ----------
        wavelength : Quantity
            Wavelength at which to interpolate α.

        Returns
        -------
        Quantity
            Absorption coefficient with units from the data file.
        """
        if self._absorption_ipol is None:
            raise ValueError(
                "Absorption data not loaded. Call load_absorption_data() first."
            )
        return self._absorption_ipol(wavelength)


def _parse_header_units(header_line: str, sep: str | None) -> tuple[Unit, Unit]:
    """Extract units from a self-describing header line.

    Expects two columns with units in parentheses,
    e.g. ``wavelength (nm)  alpha (1/um)`` or ``wavelength (nm), alpha (1/um)``.
    """
    parts = header_line.split(sep) if sep else header_line.split()

    # Reassemble tokens into two column groups based on parenthesised units.
    # Each column ends when we find a token containing ')'.
    columns: list[str] = []
    current_tokens: list[str] = []
    for token in parts:
        current_tokens.append(token.strip())
        if ")" in token:
            columns.append(" ".join(current_tokens))
            current_tokens = []
    if current_tokens:
        columns.append(" ".join(current_tokens))

    if len(columns) < 2:
        raise ValueError(f"Cannot parse two column groups from header: '{header_line}'")

    def _extract_unit(col: str) -> Unit:
        start = col.index("(")
        end = col.index(")")
        return Unit(col[start + 1 : end])

    return _extract_unit(columns[0]), _extract_unit(columns[1])


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
    sensor_params: SensorParams | None = None

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
            # line/frame duration (with binning and read blocks)
            channel.frame_duration = (
                timings.frame_duration * channel.binning / channel.read_blocks  # type: ignore[operator]
            )  # type: ignore[operator]
            # line/frame rate (with binning read blocks)
            channel.frame_rate = (
                timings.frame_rate / channel.binning * channel.read_blocks  # type: ignore[operator]
            )  # type: ignore[operator]
            # total exposure duration (with binning and blocks)
            channel.integration_duration = (
                timings.integration_duration * channel.binning  # type: ignore[operator]
            )
            # total TDI column duration
            # total tdi col duration = line/frame duration (w/bin) x Nr of TDI stages
            channel.total_tdi_col_duration = (
                timings.frame_duration * channel.binning * channel.tdi_stages  # type: ignore[operator]
            ).to(u.ms)  # type: ignore[union-attr]

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
        with_read_blocks: bool = True,
    ) -> Quantity:
        """
        Pixel read rate.

        Computed as:
         - Pushbroom type: `horiz pix (binned) x TDI stages x read blocks x line rate`
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
        with_read_blocks : bool
            Return the value with read blocks or not (valid for pushbroom only)

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
                timings.frame_rate,  # type: ignore[arg-type]
                with_binning,
                with_tdi,
                with_read_blocks,
            )
        else:
            # there are multiple channels, sum the values
            channels = self.get_channels(band_id)

            for channel in channels:
                pix_read_rate_single_chan = channel.pix_read_rate(
                    timings.frame_rate,  # type: ignore[arg-type]
                    with_binning,
                    with_tdi,
                    with_read_blocks,
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

    def get_diffusion_mtf_1d(self, wavelength: Quantity) -> MTF_Model_1D:
        """Generate diffusion MTF model for this detector at the given wavelength.

        The diffusion MTF is isotropic (same in ALT and ACT directions).

        Looks up the absorption coefficient from the sensor_params absorption
        table, then delegates to ``MTF_Model_1D.detector_diffusion()`` or
        ``MTF_Model_1D.detector_diffusion_preset()``.

        Parameters
        ----------
        wavelength : Quantity
            Wavelength at which to compute the diffusion MTF.

        Returns
        -------
        MTF_Model_1D
            Diffusion MTF model.

        Raises
        ------
        ValueError
            If ``sensor_params`` or ``absorption_data`` is not configured.
        """
        from opticks.contrast_model.mtf import MTF_Model_1D

        sp = self.sensor_params
        if sp is None:
            raise ValueError(
                "sensor_params is not configured on this Detector. "
                "Add a sensor_params block to the detector YAML."
            )

        alpha = sp.get_absorption_coeff(wavelength)

        if sp.diffusion_model is not None:
            model = next(
                (m for m in DetectorDiffusionModel if m.label == sp.diffusion_model),
                None,
            )
            if model is None:
                raise ValueError(
                    f"Unknown diffusion_model '{sp.diffusion_model}'. "
                    f"Valid values: {[m.label for m in DetectorDiffusionModel]}"
                )
            if sp.diffusion_length is None:
                raise ValueError(
                    "diffusion_length is required when using diffusion_model."
                )
            return MTF_Model_1D.detector_diffusion(
                model,
                absorption_coeff=alpha,
                diffusion_length=sp.diffusion_length,
                field_free_depth=sp.field_free_depth,
                depletion_depth=sp.depletion_depth,
                surface_recomb_velocity=sp.surface_recomb_velocity,
                diffusion_coeff=sp.diffusion_coeff,
            )
        elif sp.diffusion_preset is not None:
            preset = next(
                (p for p in DetectorDiffusionPreset if p.label == sp.diffusion_preset),
                None,
            )
            if preset is None:
                raise ValueError(
                    f"Unknown diffusion_preset '{sp.diffusion_preset}'. "
                    f"Valid values: {[p.label for p in DetectorDiffusionPreset]}"
                )
            return MTF_Model_1D.detector_diffusion_preset(
                preset,
                absorption_coeff=alpha,
                diffusion_length=sp.diffusion_length,
                field_free_depth=sp.field_free_depth,
                depletion_depth=sp.depletion_depth,
                surface_recomb_velocity=sp.surface_recomb_velocity,
                diffusion_coeff=sp.diffusion_coeff,
            )
        else:
            raise ValueError(
                "Neither diffusion_model nor diffusion_preset is set in sensor_params."
            )

    def get_det_sampling_mtf_1d(self, channel_id: str | None = None) -> MTF_Model_1D:
        """Generate detector sampling MTF model for this detector.

        The sampling MTF is a sinc function of the effective pixel pitch.
        If a channel ID is supplied, the channel's effective pitch (native pitch
        × binning factor) is used. If no channel is supplied, the native
        detector pixel pitch is used.

        Parameters
        ----------
        channel_id : str, optional
            Channel identifier. If provided, the channel's effective pixel pitch
            (accounting for binning) is used. If ``None``, the native detector
            pixel pitch is used.

        Returns
        -------
        MTF_Model_1D
            Detector sampling MTF model.
        """

        if channel_id is not None:
            pixel_pitch = self.get_channel(channel_id).pixel_pitch(with_binning=True)
        else:
            pixel_pitch = self.pixel_pitch

        return MTF_Model_1D.detector_sampling(pixel_pitch)

    def get_crosstalk_mtf_1d(self, channel_id: str | None = None) -> MTF_Model_1D:
        """Generate crosstalk MTF model for this detector.

        Uses the crosstalk coefficients from ``sensor_params`` and the pixel
        pitch (optionally adjusted for binning via *channel_id*).

        Parameters
        ----------
        channel_id : str, optional
            Channel identifier. If provided, the channel's effective pixel
            pitch (accounting for binning) is used. If ``None``, the native
            detector pixel pitch is used.

        Returns
        -------
        MTF_Model_1D
            Crosstalk MTF model.

        Raises
        ------
        ValueError
            If ``sensor_params`` or ``crosstalk_xs`` is not configured.
        """
        from opticks.contrast_model.mtf import MTF_Model_1D

        sp = self.sensor_params
        if sp is None:
            raise ValueError(
                "sensor_params is not configured on this Detector. "
                "Add a sensor_params block to the detector YAML."
            )
        if sp.crosstalk_xs is None:
            raise ValueError(
                "crosstalk_xs is not set in sensor_params. "
                "Add crosstalk_xs to the sensor_params block."
            )

        if channel_id is not None:
            pixel_pitch = self.get_channel(channel_id).pixel_pitch(with_binning=True)
        else:
            pixel_pitch = self.pixel_pitch

        xd = sp.crosstalk_xd if sp.crosstalk_xd is not None else 0.0

        return MTF_Model_1D.detector_crosstalk(pixel_pitch, sp.crosstalk_xs, xd)
