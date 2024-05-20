# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


import numpy as np
from numpy.typing import ArrayLike, NDArray
from pint import Quantity

from opticks import u
from opticks.imager_model.optics import Optics


class MTF_Model:

    def __init__(self, id: str, mtf_value_func) -> None:
        """
        Modulation Transfer Function (MTF) Model.

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
        Gets the MTF value (between 0 and 1) for the given
        input line frequency.

        Parameters
        ----------
        input_line_freq : Quantity | ArrayLike[Quantity]
            Input line frequency (in lp/mm)

        Returns
        -------
        float | NDArray[np.float64]
            MTF value (usually between 0 and 1, though can be negative)
        """
        return self._value_func(input_line_freq)

    def __str__(self) -> str:
        return self.id

    @staticmethod
    @u.check("[length]", None)
    def ideal_optics(
        wavelength: Quantity | np.ndarray[Quantity], optics: Optics
    ) -> "MTF_Model":
        """
        Ideal optical MTF model (for the given input line frequency).

        Assumes uniformly illuminated circular aperture, no significant aberrations.

        Parameters
        ----------
        wavelength : Quantity | ArrayLike[Quantity]
            Wavelength at which MTF is computed
        optics: Optics
            Optics model (to compute the spatial cutoff frequency)

        Returns
        -------
        MTF_Model
            MTF model
        """

        # set the spatial cutoff frequency
        spatial_cutoff_freq = optics.spatial_cutoff_freq(wavelength)

        # set the id
        id = f"Ideal optical MTF at {wavelength:~P}"

        # set the value function (with the fixed spatial cutoff frequency)
        def value_func(input_line_freq):
            return _ideal_optical_mtf(input_line_freq, spatial_cutoff_freq)

        return MTF_Model(id, value_func)

    @staticmethod
    @u.check("[length]")
    def detector_sampling(pixel_pitch: Quantity) -> "MTF_Model":
        """
        Detector sampling MTF model (for the given input line frequency).

        Parameters
        ----------
        pixel_pitch : Quantity
            Pixel pitch (with or without binning)

        Returns
        -------
        Quantity
            MTF value between 0 and 1
        """

        # set the id
        id = f"Detector MTF with pixel pitch {pixel_pitch:~P}"

        # set the value function (with the fixed pixel pitch)
        def value_func(input_line_freq):
            return _detector_sampling_mtf(input_line_freq, pixel_pitch)

        return MTF_Model(id, value_func)


@u.check(None, "[length]")
def _detector_sampling_mtf(
    input_line_freq: Quantity | np.ndarray[Quantity], pixel_pitch: Quantity
) -> float | NDArray[np.float64]:
    """
    Detector sampling MTF for the given input line frequency.

    Returns the MTF value between 0 and 1.

    Parameters
    ----------
    input_line_freq: Quantity | ArrayLike[Quantity]
        Input line frequency (in lp/mm)
    pixel_pitch : Quantity
        Pixel pitch (with or without binning)

    Returns
    -------
    Quantity
        MTF value between 0 and 1
    """

    # pixel pitch (um) x input line freq (lp/mm)
    a_fx = (pixel_pitch * input_line_freq / u.lp).to_reduced_units()

    # This is the alternative formulation
    # a_fx = (input_line_freq / self.nyquist_freq).to_reduced_units()/2

    # This is the alternative formulation
    # mtf_det_sampling = np.abs(np.sin(np.pi*a_fx)/(np.pi*a_fx))

    # sinc does not receive Quantity input.
    return np.sinc(a_fx.m)


def _ideal_optical_mtf(
    input_line_freq: Quantity | np.ndarray[Quantity], spatial_cutoff_freq: Quantity
) -> float | NDArray[np.float64]:
    """
    Ideal optical MTF for the given input line frequency.

    Assumes uniformly illuminated circular aperture, no significant aberrations.

    Returns the MTF value between 0 and 1.

    Parameters
    ----------
    input_line_freq: QuanQuantity | ArrayLike[Quantity]tity
        Input line frequency (in lp/mm)
    spatial_cutoff_freq: Quantity
        Spatial cutoff frequency (in lp/mm)

    Returns
    -------
    Quantity
        MTF value between 0 and 1
    """

    # normalised optical frequency
    f_ov_fc = input_line_freq / spatial_cutoff_freq

    # This is the alternative formulation
    # psi = np.arccos(f_ov_fc)
    # mtf_ideal_optical = 2/np.pi * (psi-np.cos(psi)*np.sin(psi))

    return 2 / np.pi * (np.arccos(f_ov_fc) - f_ov_fc * np.sqrt(1 - f_ov_fc**2)).m
