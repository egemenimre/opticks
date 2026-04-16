# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Package for Modulation Transfer Function (MTF) related classes and functions.

"""

from typing import Self

import numpy as np
from astropy.units import Quantity, imperial
from matplotlib import pyplot as plt
from numpy.typing import NDArray
from prysm._richdata import RichData

from opticks import u
from opticks.contrast_model.detector_mtf import (
    detector_crosstalk_mtf_1d,
    detector_diffusion_mtf,
    validate_crosstalk_params,
    validate_diffusion_params,
)
from opticks.contrast_model.optics_mtf import (
    _aberrated_optical_mtf,
    _ideal_optical_mtf,
)
from opticks.utils.math_utils import InterpolatorWithUnits, InterpolatorWithUnitTypes


class MTF_Model_1D:
    def __init__(self, id: str, mtf_value_func) -> None:
        """
        Modulation Transfer Function (MTF) Model in 1D.

        Usually this is not called by the User. There are other constructors instead.

        Parameters
        ----------
        id : str
            Model explanation or identifier
        mtf_value_func
            Function that returns the MTF value for a given input line frequency
        """

        self.id = id
        self._value_func = mtf_value_func

    def mtf_value(self, input_line_freq: Quantity) -> float | NDArray[np.float64]:
        """
        Gets the MTF value for the given input line frequency.

        Parameters
        ----------
        input_line_freq : Quantity | ArrayLike[Quantity]
            Input line frequency (in cycles/mm)

        Returns
        -------
        float | NDArray[np.float64]
            MTF value (usually between 0 and 1, though can be negative)
        """
        return self._value_func(input_line_freq)

    def __str__(self) -> str:
        return self.id

    @staticmethod
    def external_data(
        freq_values: Quantity,
        mtf_values: NDArray[np.float64],
        id: str | None = None,
    ) -> "MTF_Model_1D":
        """
        MTF Model from an external data.

        This model is mainly to ingest actual measurements or the results
        from complex simulations (such as Zemax).

        Uses the scipy 'Akima1DInterpolator' (default Akima method)
        in the background, though withy unit support via the
        `InterpolatorWithUnits` wrapper.


        Parameters
        ----------
        freq_data : ArrayLike[Quantity]
            List of spatial frequencies (in cy/mm or lp/mm)
        spatial_freqs : ArrayLike[np.float64]
            List of corresponding MTF values (should be less than or equal to 1)
        id : str
            id text associated with the data. Default text if `None`.

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the id if None
        if not id:
            id = "External MTF Data"

        # prepare interpolator
        interpolator = InterpolatorWithUnits.from_ipol_method(
            InterpolatorWithUnitTypes.AKIMA,
            freq_values,
            mtf_values,  # type: ignore[arg-type]
        )

        # set the value function (with the interpolator)
        def value_func(input_line_freq):
            return _external_data_mtf(input_line_freq, interpolator)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def ideal_optics(
        spatial_cutoff_freq: Quantity,
        wavelength: Quantity | None = None,
    ) -> "MTF_Model_1D":
        """
        Ideal optical MTF model.

        Assumes uniformly illuminated circular aperture, no aberrations.

        Parameters
        ----------
        spatial_cutoff_freq : Quantity
            Spatial cutoff frequency (in cycles/mm).
            Can be computed via `Optics.spatial_cutoff_freq(wavelength)`.
        wavelength : Quantity, optional
            Wavelength (used only for the model id string)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the id
        id = (
            f"Ideal optical MTF at {wavelength}"
            if wavelength is not None
            else "Ideal optical MTF"
        )

        # set the value function (with the fixed spatial cutoff frequency)
        def value_func(input_line_freq):
            return _ideal_optical_mtf(input_line_freq, spatial_cutoff_freq)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def emp_model_aberrated_optics(
        spatial_cutoff_freq: Quantity,
        w_rms: float | np.ndarray,
        wavelength: Quantity | None = None,
    ) -> "MTF_Model_1D":
        """
        Aberrated optical MTF model with an empirical model.

        Uses an empirical model for the optical aberrations, such that:
        MTF_true = MTF_ideal x ATF. See Shannon's The Art and Science of Optical
        Design for more information.

        The 'w_rms' value corresponds to the RMS of the total wavefront error,
        or how much the actual wavefront deviates from the ideal wavefront.
        The unit of this deviation is the multiple wavelengths
        (such as 0.15 x lambda).

        Parameters
        ----------
        spatial_cutoff_freq : Quantity
            Spatial cutoff frequency (in cycles/mm).
            Can be computed via `Optics.spatial_cutoff_freq(wavelength)`.
        w_rms : float or ndarray
            RMS of the total wavefront error (in wavelengths)
        wavelength : Quantity, optional
            Wavelength (used only for the model id string)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the id
        id = (
            f"Aberrated optical MTF at {wavelength} (W_RMS = {w_rms})"
            if wavelength is not None
            else f"Aberrated optical MTF (W_RMS = {w_rms})"
        )

        # set the value function (with the fixed spatial cutoff frequency and w_rms)
        def value_func(input_line_freq):
            return _aberrated_optical_mtf(input_line_freq, spatial_cutoff_freq, w_rms)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def from_mtf_2d(mtf_2d: RichData, slice: str) -> "MTF_Model_1D":
        """Converts the 2D MTF into a 1D MTF Model object.

        The `prysm` MTF object produces the data in 2D, as the contrast transfer
        may differ in different directions. This method extracts a slice (usually
        in x or y directions) and generates an `MTF_Model_1D` object.

        Possible slice strings are `x`, `y`, `azavg`, `azavmedian`, `azmin`, `azpv`,
        `azvar`, `azstd`. The first two are simply slices in the x and y axes.
        The remaining are different ways of sampling the data in the azimuthal direction.

        Parameters
        ----------
        mtf_2d : RichData
            2D MTF data (line freq should be in cy/mm)
        slice : str
            slice type (e.g., "x" or "y" direction)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # get the slice
        fx, mtf_values = getattr(mtf_2d.slices(twosided=False), slice)

        if not isinstance(mtf_2d.dx, Quantity):
            # MTF 2D has no units, add units to fx
            fx = fx * u.cy / u.mm

        # generate the model
        return MTF_Model_1D.external_data(fx, mtf_values, id=f"MTF in {slice}")

    @staticmethod
    def detector_sampling(pixel_pitch: Quantity) -> "MTF_Model_1D":
        """
        Detector sampling MTF model.

        MTF value is usually between 0 and 1, though contrast reversal
        may result in negative values.

        Parameters
        ----------
        pixel_pitch : Quantity
            Pixel pitch (with or without binning)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the id
        id = f"Detector MTF with pixel pitch {pixel_pitch}"

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return _detector_sampling_mtf(input_line_freq, pixel_pitch)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def motion_blur(pixel_pitch: Quantity, blur_extent: float) -> "MTF_Model_1D":
        """
        Motion blur MTF model.

        Parameters
        ----------
        pixel_pitch : Quantity
            Pixel pitch (with or without binning)
        blur_extent : float
            Blur extent on the image, given in pixels (e.g. 0.1 or 10% pix)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the id
        id = (
            f"Motion Blur MTF with pixel pitch {pixel_pitch}"
            f"and blur extent {blur_extent:.6f}"
        )

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return _smear_mtf(input_line_freq, pixel_pitch, blur_extent)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def smear(pixel_pitch: Quantity, blur_extent: float) -> "MTF_Model_1D":
        """
        Drift/Smear MTF model.

        The drift/smear in the along-track and across-track dimensions
        should be computed separately.

        Parameters
        ----------
        pixel_pitch : Quantity
            Pixel pitch (with or without binning)
        blur_extent : float
            Blur extent on the image, given in pixels (e.g. 0.1 or 10% pix)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the id
        id = (
            f"Drift/Smear MTF with pixel pitch {pixel_pitch} "
            f"and blur extent {blur_extent:.6f}"
        )

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return _smear_mtf(input_line_freq, pixel_pitch, blur_extent)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def jitter(
        pixel_pitch: Quantity,
        jitter_stdev: float,
    ) -> "MTF_Model_1D":
        """
        Jitter MTF model.

        The `jitter_stdev` value is defined as the 1 sigma value of the jitter
        amplitude, defined in pixels (e.g. 10% of the pixel).

        Jitter is defined with respect to the relevant frequency of the imaging problem.
        For example an imaging system with a 10 msec integration time will have the
        disturbances at about 100 Hz or higher registered as jitter. Slower ones will be
        classified as drift/smear.

        Parameters
        ----------
        pixel_pitch : Quantity
            Pixel pitch (with or without binning)
        jitter_stdev : float
            Standard deviation of the jitter value in pixels (e.g. 0.1 or 10% pix)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the id
        id = (
            f"Jitter MTF with pixel pitch {pixel_pitch} and"
            f"jitter standard deviation {jitter_stdev:.6f}"
        )

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return _jitter_mtf(input_line_freq, pixel_pitch, jitter_stdev)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def combined(*mtf_models: "MTF_Model_1D") -> "MTF_Model_1D":
        """
        Combination MTF models.

        The combined MTF is useful when describing a combination of
        multiple MTF Models. For example, the Imager MTF is a combination
        of Optical MTF, Detector Sampling MTF and Detector Diffusion MTF.

        Parameters
        ----------
        mtf_models : tuple containing multiple "MTF_Model_1D" objects
            list of MTF Models to be combined

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the id
        id = "A combination of multiple MTF Models."

        # build list of value functions
        value_funcs = [mtf_model._value_func for mtf_model in mtf_models]  # type: ignore[union-attr]

        # set the value function
        def value_func(input_line_freq):
            return _combined_mtf(input_line_freq, value_funcs)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def fixed(mtf_value: float) -> "MTF_Model_1D":
        """
        Fixed value MTF model.

        Used for disturbances and imperfections.

        Parameters
        ----------
        mtf_value : float
            MTF value over the spatial frequency domain (between 0 and 1, inclusive)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """
        # check mtf value
        if mtf_value < 0 or mtf_value > 1:
            raise ValueError(
                f"Fixed value MTF model - input value should be between 0 and 1: {mtf_value}"
            )

        # set the id
        id = f"Fixed MTF with value {mtf_value:.6f}"

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return mtf_value

        return MTF_Model_1D(id, value_func)

    @u.quantity_input(
        diffusion_length="length",
        field_free_depth="length",
        depletion_depth="length",
        surface_recomb_velocity="speed",
    )
    @staticmethod
    def detector_diffusion(
        model,
        absorption_coeff: Quantity,
        diffusion_length: Quantity,
        field_free_depth: Quantity | None = None,
        depletion_depth: Quantity | None = None,
        surface_recomb_velocity: Quantity | None = None,
        diffusion_coeff: Quantity | None = None,
    ):
        """
        Detector diffusion MTF model (Crowell & Labuda 1969).

        Computes the diffusion-limited MTF for a given detector geometry and
        boundary conditions. Use one of the five simplified model variants
        (BSI-1/2/3, FSI-1/2) matching your detector type.

        Parameters
        ----------
        model : DetectorDiffusionModel
            Diffusion model variant (geometry + boundary conditions)
        absorption_coeff : Quantity["1/length"]
            Optical absorption coefficient (e.g., in 1/µm or 1/cm)
        diffusion_length : Quantity["length"]
            Minority carrier diffusion length L_o
        field_free_depth : Quantity["length"], optional
            Field-free (bulk) region depth L_a. Required for BSI-1/2/3, FSI-1.
        depletion_depth : Quantity["length"], optional
            Depletion region depth L_b (BSI) or L_D (FSI). Required for
            BSI-1, BSI-2, FSI-1, FSI-2.
        surface_recomb_velocity : Quantity["speed"], optional
            Surface recombination velocity S. Required for BSI-3 only.
        diffusion_coeff : Quantity, optional
            Minority carrier diffusion coefficient D (provide in cm²/s).
            Required for BSI-3 only.

        Returns
        -------
        MTF_Model_1D
            Diffusion MTF model
        """
        # validate forbidden/required parameters per model
        validate_diffusion_params(
            model,
            field_free_depth,
            depletion_depth,
            surface_recomb_velocity,
            diffusion_coeff,
        )

        # convert to internal float units (µm, µm⁻¹)
        alpha = absorption_coeff.to(1 / u.um).value
        L_o = diffusion_length.to(u.um).value
        L_a = field_free_depth.to(u.um).value if field_free_depth is not None else None
        depth = depletion_depth.to(u.um).value if depletion_depth is not None else None
        s_over_d = (
            (surface_recomb_velocity / diffusion_coeff).decompose().to(1 / u.um).value
            if surface_recomb_velocity is not None
            else None
        )

        id_str = (
            f"Detector Diffusion MTF ({model}): "
            f"alpha={alpha:.4g}/µm, L_o={L_o:.4g}µm"
            + (f", L_a={L_a:.4g}µm" if L_a is not None else "")
            + (f", depth={depth:.4g}µm" if depth is not None else "")
            + (f", s/d={s_over_d:.4g}/µm" if s_over_d is not None else "")
        )

        def value_func(input_line_freq):
            return detector_diffusion_mtf(
                input_line_freq.to(u.cy / u.mm).value,
                model,
                alpha,
                L_o,
                L_a,
                depth,
                s_over_d,
            )

        return MTF_Model_1D(id_str, value_func)

    @u.quantity_input(
        diffusion_length="length",
        field_free_depth="length",
        depletion_depth="length",
        surface_recomb_velocity="speed",
    )
    @staticmethod
    def detector_diffusion_preset(
        preset,
        absorption_coeff: Quantity,
        diffusion_length: Quantity | None = None,
        field_free_depth: Quantity | None = None,
        depletion_depth: Quantity | None = None,
        surface_recomb_velocity: Quantity | None = None,
        diffusion_coeff: Quantity | None = None,
    ):
        """
        Detector diffusion MTF model from a predefined detector category preset.

        Looks up default parameters for the preset, then applies any supplied
        overrides. Raises ``ValueError`` if an override is incompatible with
        the preset's model (e.g., supplying ``surface_recomb_velocity`` for a
        BSI-1 preset).

        Parameters
        ----------
        preset : DetectorDiffusionPreset
            Detector category preset
        absorption_coeff : Quantity["1/length"]
            Optical absorption coefficient
        diffusion_length : Quantity["length"], optional
            Override for the preset diffusion length L_o
        field_free_depth : Quantity["length"], optional
            Override for the preset field-free depth L_a
        depletion_depth : Quantity["length"], optional
            Override for the preset depletion/substrate depth
        surface_recomb_velocity : Quantity["speed"], optional
            Override for the preset surface recombination velocity (BSI-3 only)
        diffusion_coeff : Quantity, optional
            Override for the preset diffusion coefficient (BSI-3 only)

        Returns
        -------
        MTF_Model_1D
            Diffusion MTF model
        """
        model = preset.model

        # validate overrides against the preset's model
        overrides = {
            "field_free_depth": field_free_depth,
            "depletion_depth": depletion_depth,
            "surface_recomb_velocity": surface_recomb_velocity,
            "diffusion_coeff": diffusion_coeff,
        }
        for param_name, value in overrides.items():
            if value is not None and param_name not in model.required_params:
                raise ValueError(
                    f"Parameter '{param_name}' is not used by model {model} "
                    f"(preset {preset}). Remove it or choose a different preset."
                )

        # merge overrides into defaults (user-supplied wins)
        merged = dict(preset.params)
        if diffusion_length is not None:
            merged["diffusion_length"] = diffusion_length
        if field_free_depth is not None:
            merged["field_free_depth"] = field_free_depth
        if depletion_depth is not None:
            merged["depletion_depth"] = depletion_depth
        if surface_recomb_velocity is not None:
            merged["surface_recomb_velocity"] = surface_recomb_velocity
        if diffusion_coeff is not None:
            merged["diffusion_coeff"] = diffusion_coeff

        return MTF_Model_1D.detector_diffusion(
            model=model,
            absorption_coeff=absorption_coeff,
            diffusion_length=merged["diffusion_length"],
            field_free_depth=merged.get("field_free_depth"),
            depletion_depth=merged.get("depletion_depth"),
            surface_recomb_velocity=merged.get("surface_recomb_velocity"),
            diffusion_coeff=merged.get("diffusion_coeff"),
        )

    @staticmethod
    def detector_crosstalk(
        pixel_pitch: Quantity,
        crosstalk_xs: float,
        crosstalk_xd: float = 0.0,
    ) -> "MTF_Model_1D":
        """
        Detector crosstalk MTF model (center pixel, 8 neighbours).

        Computes the MTF due to nearest-neighbour charge sharing using
        the discrete kernel model with separate side and diagonal
        crosstalk coefficients.

        The 1D MTF (fy=0 slice) is:
        ``MTF(f) = 1 - 2(Xs + 2Xd)(1 - cos(2π f P))``

        When ``crosstalk_xd = 0``, this reduces to the classical
        nearest-neighbour formula ``MTF(f) = 1 - 2Xs(1 - cos(2π f P))``.

        Parameters
        ----------
        pixel_pitch : Quantity["length"]
            Pixel pitch
        crosstalk_xs : float
            Side-neighbour crosstalk coefficient (dimensionless fraction,
            e.g. 0.03 for 3%).
        crosstalk_xd : float, optional
            Diagonal-neighbour crosstalk coefficient (dimensionless fraction,
            e.g. 0.005 for 0.5%). Default is 0.0 (no diagonal crosstalk).

        Returns
        -------
        MTF_Model_1D
            Crosstalk MTF model
        """
        validate_crosstalk_params(crosstalk_xs, crosstalk_xd)

        id_str = (
            f"Detector Crosstalk MTF: "
            f"Xs={crosstalk_xs:.4g}, Xd={crosstalk_xd:.4g}, "
            f"pitch={pixel_pitch}"
        )

        def value_func(input_line_freq):
            return detector_crosstalk_mtf_1d(
                input_line_freq,
                crosstalk_xs,
                crosstalk_xd,
                pixel_pitch,
            )

        return MTF_Model_1D(id_str, value_func)


def _force_return_float(mtf_value):
    if isinstance(mtf_value, Quantity):
        return mtf_value.decompose().value
    else:
        return mtf_value


def _combined_mtf(
    input_line_freq: Quantity, value_funcs
) -> float | NDArray[np.float64]:
    """
    Combination of multiple MTF Model MTF data for the given
    input line frequency.

    Parameters
    ----------
    input_line_freq: Quantity | ArrayLike[Quantity]
        Input line frequency (in cycles/mm)
    value_funcs
        list of value functions for all the MTF Models

    Returns
    -------
    float | NDArray[np.float64]
        MTF value between 0 and 1
    """

    mtf_value = 1
    for value_func in value_funcs:
        mtf_value = mtf_value * value_func(input_line_freq)

    # force return float
    return _force_return_float(mtf_value)


def _external_data_mtf(
    input_line_freq: Quantity,
    interpolator: InterpolatorWithUnits,
) -> float | NDArray[np.float64]:
    """
    Interpolated MTF data for the given input line frequency.

    Note that the interpolator deletes the units and
    therefore has no units support.

    Parameters
    ----------
    input_line_freq: Quantity | ArrayLike[Quantity]
        Input line frequency (in cycles/mm)
    interpolator: InterpolatorWithUnits
        Interpolator

    Returns
    -------
    float | NDArray[np.float64]
        MTF value between 0 and 1
    """

    # Compute the interpolated value
    mtf_value = interpolator(input_line_freq)

    # force return float
    return _force_return_float(mtf_value)


@u.quantity_input(pixel_pitch="length")
def _detector_sampling_mtf(
    input_line_freq: Quantity, pixel_pitch: Quantity
) -> float | NDArray[np.float64]:
    """
    Detector sampling MTF for the given input line frequency.

    MTF value is usually between 0 and 1, though contrast reversal
    may result in negative values.

    Parameters
    ----------
    input_line_freq: Quantity | ArrayLike[Quantity]
        Input line frequency (in cycles/mm)
    pixel_pitch : Quantity
        Pixel pitch (with or without binning)

    Returns
    -------
    float | NDArray[np.float64]
        MTF value (usually between 0 and 1, though can be negative)
    """

    # pixel pitch (um) x input line freq (cycles/mm)
    a_fx = (pixel_pitch * input_line_freq / u.cy).decompose()

    # This is the alternative formulation
    # a_fx = (input_line_freq / self.nyquist_freq).to_reduced_units()/2

    # This is the alternative formulation (negative values possible)
    # return np.sin(np.pi * a_fx) / (np.pi * a_fx)

    with u.set_enabled_equivalencies(u.dimensionless_angles()):
        mtf_value = np.sinc(a_fx)

    # force return float
    return _force_return_float(mtf_value)


@u.quantity_input(pixel_pitch="length")
def _smear_mtf(
    input_line_freq: Quantity,
    pixel_pitch: Quantity,
    blur_extent: float | Quantity,
) -> float | NDArray[np.float64]:
    """
    Drift/Smear MTF for the given input line frequency.

    MTF value is usually between 0 and 1, though contrast reversal
    may result in negative values.

    Parameters
    ----------
    input_line_freq: Quantity | ArrayLike[Quantity]
        Input line frequency (in cycles/mm)
    pixel_pitch : Quantity
        Pixel pitch (with or without binning)
    blur_extent : float | Quantity
        Blur extent on the image, given in pixels (e.g. 0.1 or 10% pix)

    Returns
    -------
    float | NDArray[np.float64]
        MTF value (usually between 0 and 1, though can be negative)
    """
    # pixel pitch (um) x input line freq (cycles/mm)
    a_fx = pixel_pitch * input_line_freq / u.cy

    mtf_value = np.sinc((blur_extent * a_fx).decompose().value)  # type: ignore[union-attr]

    # force return float
    return _force_return_float(mtf_value)


@u.quantity_input(pixel_pitch="length")
def _jitter_mtf(
    input_line_freq: Quantity,
    pixel_pitch: Quantity,
    jitter_stdev: float | Quantity,
) -> float | NDArray[np.float64]:
    """
    Jitter MTF for the given input line frequency.

    The `jitter_stdev` value is defined as the 1 sigma value of the jitter amplitude,
    defined in pixels (e.g. 10% of the pixel).

    Jitter is defined with respect to the relevant frequency of the imaging problem.
    For example an imaging system with a 10 msec integration time will have the
    disturbances at about 100 Hz or higher registered as jitter. Slower ones will be
    classified as drift/smear.

    Returns the MTF value between 0 and 1.

    Parameters
    ----------
    input_line_freq: Quantity | ArrayLike[Quantity]
        Input line frequency (in cycles/mm)
    pixel_pitch : Quantity
        Pixel pitch (with or without binning)
    jitter_stdev : float | Quantity
        Standard deviation of the jitter value in pixels

    Returns
    -------
    float | NDArray[np.float64]
        MTF value between 0 and 1
    """

    # pixel pitch (um) x input line freq (cycles/mm)
    a_fx = pixel_pitch * input_line_freq / u.cy

    mtf_value = np.exp(-2 * ((np.pi * jitter_stdev * a_fx).decompose() ** 2))  # type: ignore[union-attr]

    # force return float
    return _force_return_float(mtf_value)


class MTF_Plot_1D:  # pragma: no cover
    """Generates an MTF Plot.

    Each MTF Model is used to generate the plot y values,
    using the `freq_list` as the discrete x axis values
    and the dict key as the label.

    Use the `set_plot_style` method to further detail the
    style options and decorators like title and axis labels.

    Parameters
    ----------
    freq_list : Arraylike
        list of spatial frequency values
    mtf_data : dict[str, MTF_Model_1D]
        list of MTF models and labels
    acceptable_limit : float, optional
        horizontal "acceptable limit" line value, by default 0.1
    nyq_limit : Quantity, optional
        spatial frequency corresponding to the Nyquist limit, by default None
    """

    def __init__(
        self,
        freq_list,
        mtf_data: dict[str, MTF_Model_1D],
        acceptable_limit: float = 0.1,
        nyq_limit: Quantity | None = None,
    ) -> None:

        self.fig, self.ax = plt.subplots()

        self.populate_plot(freq_list, mtf_data, acceptable_limit, nyq_limit)

    def populate_plot(
        self,
        freq_list,
        mtf_data: dict[str, MTF_Model_1D],
        acceptable_limit: float = 0.1,
        nyq_limit: Quantity | None = None,
    ) -> Self:
        """
        Populates the MTF plot lines using the MTF Models.

        This conveniently adds items in addition to those in the constructor.

        Each MTF Model is used to generate the plot y values,
        using the `freq_list` as the discrete x axis values
        and the dict key as the label.

        If the `freq_list` has a single column (`len(freq_list.shape) == 1`)
        then the same list is used for each MTF Model. If the `freq_list` has
        as many columns as the `mtf_data` (`len(freq_list) == len(mtf_data)`)
        then each column is used for the consecutive MTF Model in the dict.
        The order of the frequency list should therefore match the order of the
        models in the dict. The dict implementation keeps the insertion order
        since Python 3.6.

        Parameters
        ----------
        freq_list : Arraylike
            list of spatial frequency values
        mtf_data : dict[str, MTF_Model_1D]
            list of MTF models and labels
        acceptable_limit : float, optional
            horizontal "acceptable limit" line value, by default 0.1
        nyq_limit : Quantity, optional
            spatial frequency corresponding to the Nyquist limit, by default None

        Returns
        -------
        MTF_Plot_1D
            self object for convenience
        """

        if len(freq_list) == len(mtf_data):
            # one freq list for each

            # generate MTF data lines
            for label, freqs in zip(mtf_data, freq_list, strict=True):
                # generate values (y axis)
                mtf_values = mtf_data[label].mtf_value(freqs)

                # generate plot line
                self.ax.plot(freqs, mtf_values, label=label)

        elif len(freq_list.shape) == 1:
            # one freq list for all

            # generate MTF data lines
            for label, mtf_model in mtf_data.items():
                # generate values (y axis)
                mtf_values = mtf_model.mtf_value(freq_list)

                # generate plot line
                self.ax.plot(freq_list, mtf_values, label=label)
        else:
            raise ValueError(
                f"Columns of the frequency list ({len(freq_list.shape)}) does not match "
                f"the columns of the MTF data list ({len(mtf_data)})"
            )

        # -----------------

        if acceptable_limit:
            self.ax.axhline(
                acceptable_limit,
                label="Acceptable MTF Limit",
                linestyle="-.",
            )

        # ax.axhline(26400 *ureg.feet, color='tab:red')
        # ax.axvline(120* ureg.minutes, color='tab:green')

        if nyq_limit is not None:
            self.ax.axvline(
                nyq_limit.value,
                label="Detector Nyq Limit",
                linestyle="--",
            )

        return self

    def set_plot_style(
        self,
        x_max=None,
        y_min=0,
        title=None,
        xlabel="input line freq (cycles/mm)",
        ylabel="MTF",
        height: Quantity | float = 10 * u.cm,
        width: Quantity | float = 15 * u.cm,
    ) -> Self:
        """
        Sets some default style parameters for MTF plots.

        Parameters
        ----------
        x_max : float, optional
            max value of x axis, by default None
        y_min : float, optional
            minimum value of y axis, by default 0
        title : str, optional
            title of the plot, by default None
        xlabel : str, optional
            x-axis label, by default "input line freq (cycles/mm)"
        ylabel : str, optional
            y-axis label, by default "MTF"
        height : int, optional
            height of the figure (in cm), by default 10
        width : int, optional
            width of the figure (in cm), by default 15

        Returns
        -------
        MTF_Plot_1D
            self object for convenience
        """

        # set decorators
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        if title:
            self.ax.set_title(title)

        self.fig.legend(bbox_to_anchor=(1, 1), loc="upper left")

        # set plot formatting
        self.ax.xaxis.grid(True)
        self.ax.yaxis.grid(True)
        if x_max is not None:
            self.ax.set_xlim(0, x_max)
        self.ax.set_ylim(y_min, 1)

        self.fig.tight_layout()

        with imperial.enable():
            if isinstance(height, Quantity):
                height = height.to_value(imperial.inch)  # type: ignore[assignment]
            if isinstance(width, Quantity):
                width = width.to_value(imperial.inch)  # type: ignore[assignment]

        self.fig.set_figheight(height)  # type: ignore[arg-type]
        self.fig.set_figwidth(width)  # type: ignore[arg-type]

        return self

    # def show_plot(self):
    #     plt.show()
