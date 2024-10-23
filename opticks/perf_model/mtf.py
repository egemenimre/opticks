# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Package for Modulation Transfer Function (MTF) related classes and functions.

"""

from typing import Self

import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import NDArray
from pint import Quantity
from prysm._richdata import RichData
from prysm.otf import mtf_from_psf

from opticks import u
from opticks.imager_model.optics import Optics
from opticks.utils.math_utils import InterpolatorWithUnits, InterpolatorWithUnitTypes
from opticks.utils.prysm_utils import richdata_with_units


class MTF_Model_1D:

    def __init__(self, id: str, mtf_value_func) -> None:
        """
        Modulation Transfer Function (MTF) Model in 1D.

        Parameters
        ----------
        id : str
            Model explanation or identifier
        mtf_value_func
            Function that returns the MTF value for a given input line frequency
        """

        self.id = id
        self._value_func = mtf_value_func

    def mtf_value(
        self, input_line_freq: Quantity | np.ndarray[Quantity]
    ) -> float | NDArray[np.float64]:
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
        freq_values: np.ndarray[Quantity],
        mtf_values: np.ndarray[np.float64],
        id: str = None,
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
            InterpolatorWithUnitTypes.AKIMA, freq_values, mtf_values
        )

        # set the value function (with the interpolator)
        def value_func(input_line_freq):
            return _external_data_mtf(input_line_freq, interpolator)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    @u.check("[length]", None)
    def ideal_optics(
        wavelength: Quantity | np.ndarray[Quantity], optics: Optics
    ) -> "MTF_Model_1D":
        """
        Ideal optical MTF model.

        Assumes uniformly illuminated circular aperture, no aberrations.

        Parameters
        ----------
        wavelength : Quantity | ArrayLike[Quantity]
            Wavelength at which MTF is computed
        optics: Optics
            Optics model (to compute the spatial cut-off frequency)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the spatial cutoff frequency
        spatial_cutoff_freq = optics.spatial_cutoff_freq(wavelength)

        # set the id
        id = f"Ideal optical MTF at {wavelength:~P}"

        # set the value function (with the fixed spatial cutoff frequency)
        def value_func(input_line_freq):
            return _ideal_optical_mtf(input_line_freq, spatial_cutoff_freq)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    @u.check("[length]", None, None)
    def emp_model_aberrated_optics(
        wavelength: Quantity | np.ndarray[Quantity],
        w_rms: float,
        optics: Optics,
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
        wavelength : Quantity | ArrayLike[Quantity]
            Wavelength at which MTF is computed
        w_rms : float
            RMS of the total wavefront error (in wavelengths)
        optics: Optics
            Optics model (to compute the spatial cutoff frequency)

        Returns
        -------
        MTF_Model_1D
            MTF model
        """

        # set the spatial cutoff frequency
        spatial_cutoff_freq = optics.spatial_cutoff_freq(wavelength)

        # set the id
        id = f"Aberrated optical MTF at {wavelength:~P} (W_RMS = {w_rms})"

        # set the value function (with the fixed spatial cutoff frequency and w_rms)
        def value_func(input_line_freq):
            return _aberrated_optical_mtf(input_line_freq, spatial_cutoff_freq, w_rms)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def from_mtf_2d(mtf_2d: RichData, slice: str) -> "MTF_Model_1D":
        """Converts the 2D MTF into a 1D MTF Model object.

        The `prysm` MTF object produces the data in 2D,
        as the contrast transfer may differ in different directions.
        This method extracts a slice (usually in x or y directions)
        and generates an `MTF_Model_1D` object.

        Possible slice strings are `x`, `y`, `azavg`, `azavmedian`,
        `azmin`, `azpv`, `azvar`, `azstd`. The first two are simply
        slices in the x and y axes. The remaining are different
        ways of sampling the data in the azimuthal direction.

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
    @u.check("[length]")
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
        id = f"Detector MTF with pixel pitch {pixel_pitch:~P}"

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return _detector_sampling_mtf(input_line_freq, pixel_pitch)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    @u.check("[length]", None)
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
            f"Motion Blur MTF with pixel pitch {pixel_pitch:~P}"
            f"and blur extent {blur_extent:.6f}"
        )

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return _smear_mtf(input_line_freq, pixel_pitch, blur_extent)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    @u.check("[length]", None)
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
            f"Drift/Smear MTF with pixel pitch {pixel_pitch:~P}"
            f"and blur extent {blur_extent:.6f}"
        )

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return _smear_mtf(input_line_freq, pixel_pitch, blur_extent)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    @u.check("[length]", None)
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
            f"Jitter MTF with pixel pitch {pixel_pitch:~P} and"
            f"jitter standard deviation {jitter_stdev:.6f}"
        )

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return _jitter_mtf(input_line_freq, pixel_pitch, jitter_stdev)

        return MTF_Model_1D(id, value_func)

    @staticmethod
    def combined(*mtf_models: tuple["MTF_Model_1D", ...]) -> "MTF_Model_1D":
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
        value_funcs = [mtf_model._value_func for mtf_model in mtf_models]

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


def _combined_mtf(
    input_line_freq: Quantity | np.ndarray[Quantity], value_funcs
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

    return mtf_value


def _external_data_mtf(
    input_line_freq: Quantity | np.ndarray[Quantity],
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

    # Return the interpolated value
    return interpolator(input_line_freq)


def _ideal_optical_mtf(
    input_line_freq: Quantity | np.ndarray[Quantity], spatial_cutoff_freq: Quantity
) -> float | NDArray[np.float64]:
    """
    Ideal optical MTF for the given input line frequency.

    Assumes uniformly illuminated circular aperture, no significant aberrations.

    Returns the MTF value between 0 and 1.

    Parameters
    ----------
    input_line_freq: Quantity | ArrayLike[Quantity]
        Input line frequency (in cycles/mm)
    spatial_cutoff_freq: Quantity
        Spatial cutoff frequency (in cycles/mm)

    Returns
    -------
    float | NDArray[np.float64]
        MTF value between 0 and 1
    """

    # normalised optical frequency
    nu = input_line_freq / spatial_cutoff_freq

    # This is the alternative formulation
    # psi = np.arccos(f_ov_fc)
    # mtf_ideal_optical = 2/np.pi * (psi-np.cos(psi)*np.sin(psi))

    return 2 / np.pi * (np.arccos(nu) - nu * np.sqrt(1 - nu**2)).m


def _aberrated_optical_mtf(
    input_line_freq: Quantity | np.ndarray[Quantity],
    spatial_cutoff_freq: Quantity,
    w_rms: float,
) -> float | NDArray[np.float64]:
    """
    Aberrated optical MTF for the given input line frequency.

    Assumes an empirical aberration model for the overall
    RMS Wavefront Error.

    Returns the MTF value (usually between 0 and 1, though can be negative).

    Parameters
    ----------
    input_line_freq: Quantity | ArrayLike[Quantity]
        Input line frequency (in cycles/mm)
    spatial_cutoff_freq: Quantity
        Spatial cutoff frequency (in cycles/mm)
    w_rms : float
        RMS of the total wavefront error (in wavelengths)

    Returns
    -------
    float | NDArray[np.float64]
        MTF value (usually between 0 and 1, though can be negative).
    """

    # ideal mtf
    mtf_ideal = _ideal_optical_mtf(input_line_freq, spatial_cutoff_freq)

    # aberration transfer factor
    atf = _aberration_transfer_factor(input_line_freq, spatial_cutoff_freq, w_rms)

    return mtf_ideal * atf


def _aberration_transfer_factor(
    input_line_freq: Quantity | np.ndarray[Quantity],
    spatial_cutoff_freq: Quantity,
    w_rms: float,
) -> float | NDArray[np.float64]:
    """
    Aberration Transfer Factor (ATF) for the given input line frequency.

    Computes an empirical model for the optical aberrations, such that:
    MTF_true = MTF_ideal x ATF. See Shannon's The Art and Science of Optical
    Design for more information.

    The `w_rms` value corresponds to the RMS of the total wavefront error,
    or how much the actual wavefront deviates from the ideal wavefront.
    The unit of this deviation is the multiple wavelengths (such as 0.15 x lambda).

    Returns the ATF value (usually between 0 and 1, though can be negative).

    Parameters
    ----------
    input_line_freq: Quantity | ArrayLike[Quantity]
        Input line frequency (in cycles/mm)
    spatial_cutoff_freq: Quantity
        Spatial cutoff frequency (in cycles/mm)
    w_rms : float
        RMS of the total wavefront error (in wavelengths)

    Returns
    -------
    float | NDArray[np.float64]
        ATF value (<1)
    """

    # normalised optical frequency
    nu = input_line_freq / spatial_cutoff_freq

    # return ATF
    return (1 - ((w_rms / 0.18) ** 2 * (1 - 4 * (nu - 0.5) ** 2))).m


@u.check(None, "[length]")
def _detector_sampling_mtf(
    input_line_freq: Quantity | np.ndarray[Quantity], pixel_pitch: Quantity
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
    a_fx = (pixel_pitch * input_line_freq / u.cy).to_reduced_units()

    # This is the alternative formulation
    # a_fx = (input_line_freq / self.nyquist_freq).to_reduced_units()/2

    # This is the alternative formulation (negative values possible)
    # return np.sin(np.pi * a_fx) / (np.pi * a_fx)

    # sinc does not receive Quantity input.
    return np.sinc(a_fx.m)


@u.check(None, "[length]", None)
def _smear_mtf(
    input_line_freq: Quantity | np.ndarray[Quantity],
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
    a_fx = (pixel_pitch * input_line_freq / u.cy).to_reduced_units()

    # sinc does not receive Quantity input.
    if isinstance(blur_extent, Quantity):
        # blur extent is Quantity
        return np.sinc(blur_extent.to_reduced_units().m * a_fx.m)
    else:
        # blur extent is float
        return np.sinc(blur_extent * a_fx.m)


@u.check(None, "[length]", None)
def _jitter_mtf(
    input_line_freq: Quantity | np.ndarray[Quantity],
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
    a_fx = (pixel_pitch * input_line_freq / u.cy).to_reduced_units()

    if isinstance(jitter_stdev, Quantity):
        # jitter std dev is Quantity
        return np.exp(-2 * ((np.pi * jitter_stdev.to_reduced_units() * a_fx) ** 2))
    else:
        # jitter std dev is float
        return np.exp(-2 * ((np.pi * jitter_stdev * a_fx) ** 2))


def psf_to_mtf(psf: RichData, with_units=False) -> RichData:
    """
    Computes the MTF.

    This is the Modulation Transfer Function for the PSF.

    `prysm` does not work with units, but MTF with units
    can be generated when the method is called `with_units=True`.

    Parameters
    ----------
    psf: RichData
        PSF data
    with_units : bool, optional
        create the MTF with units

    Returns
    -------
    RichData
        Modulation Transfer Function (MTF) with spacing in cy/mm
    """

    if isinstance(psf.dx, Quantity):
        # psf has units, create one without for safety
        psf_copy = RichData(psf.data, psf.dx.m_as(u.um), None)

        mtf = mtf_from_psf(psf_copy)

        if psf.wavelength:
            mtf.wavelength = psf.wavelength.m_as(u.um)
    else:
        # psf does not have units, create mtf directly
        mtf = mtf_from_psf(psf)

        if psf.wavelength:
            mtf.wavelength = psf.wavelength

    # the resulting mtf is guaranteed to be without units
    # output dx in cy/mm

    if with_units:
        # add units to MTF
        return richdata_with_units(mtf, dx_units=u.cy / u.mm)
    else:
        # mtf returned without units
        return mtf


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
        nyq_limit: Quantity = None,
    ) -> None:

        self.fig, self.ax = plt.subplots()

        self.populate_plot(freq_list, mtf_data, acceptable_limit, nyq_limit)

    def populate_plot(
        self,
        freq_list,
        mtf_data: dict[str, MTF_Model_1D],
        acceptable_limit: float = 0.1,
        nyq_limit: Quantity = None,
    ) -> Self:
        """
        Populates the MTF plot lines using the MTF Models.

        This conveniently adds items in addition to those
        in the constructor.

        Each MTF Model is used to generate the plot y values,
        using the `freq_list` as the discrete x axis values
        and the dict key as the label.

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

        # generate MTF data lines
        for label, mtf_model in mtf_data.items():

            # generate values (y axis)
            mtf_values = mtf_model.mtf_value(freq_list)

            # generate plot line
            self.ax.plot(freq_list, mtf_values, label=label)

        # -----------------

        if acceptable_limit:
            self.ax.axhline(
                acceptable_limit,
                label="Acceptable MTF Limit",
                linestyle="-.",
            )

        # ax.axhline(26400 *ureg.feet, color='tab:red')
        # ax.axvline(120* ureg.minutes, color='tab:green')

        if nyq_limit:
            self.ax.axvline(
                nyq_limit.m,
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
        height=4,
        width=8,
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
            height of the figure (in inches), by default 4
        width : int, optional
            width of the figure (in inches), by default 6

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

        self.fig.legend()

        # set plot formatting
        self.ax.xaxis.grid(True)
        self.ax.yaxis.grid(True)
        if x_max:
            self.ax.set_xlim(0, x_max)
        self.ax.set_ylim(y_min, 1)

        self.fig.tight_layout()

        self.fig.set_figheight(height)
        self.fig.set_figwidth(width)

        return self

    # def show_plot(self):
    #     plt.show()
