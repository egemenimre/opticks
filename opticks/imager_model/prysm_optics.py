# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


import copy

from numpy import ndarray
from pint import Quantity
from prysm._richdata import RichData
from prysm.coordinates import cart_to_polar, make_xy_grid
from prysm.geometry import circle
from prysm.otf import mtf_from_psf
from prysm.propagation import Wavefront

from opticks import u


class ApertureFactory:
    """Generates oft used aperture models for `prysm`."""

    @classmethod
    @u.wraps(None, (None, u.mm, u.mm), False)
    def circle_aperture(cls, aperture_diam: Quantity | float, samples):
        """Circle aperture model for `prysm`.

        This is valid for many if not most refractive optics.

        The output is an `ndarray` of bools, therefore the unit input
        is not important for the aperture generation.

        Parameters
        ----------
        aperture_diam : Quantity | float
            aperture diameter in mm
        samples : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions.  Order is numpy axis convention, (row, col)

        Returns
        -------
        numpy.ndarray
            binary ndarray representation of the mask
        """

        # TODO docstring return eksik, r,t ekle

        # aperture radius
        aperture_radius = aperture_diam / 2.0

        # generate the square grid in polar coords
        r, t = _generate_grid(aperture_diam, samples)

        # generate aperture (circle mask applied to the square grid)
        return circle(aperture_radius, r), r, t

    @classmethod
    @u.wraps(None, (None, u.mm, None, u.mm), False)
    def circle_aperture_with_obscuration(
        cls, aperture_diam: Quantity | float, obscuration_ratio: float, samples
    ):
        """Circle aperture model for `prysm` with circular centre obscuration.

        This is valid for many reflective telescopes. The obscuration
        is usually the secondary mirror.

        The output is an `ndarray` of bools, therefore the unit input
        is not important for the aperture generation.

        Parameters
        ----------
        aperture_diam : Quantity | float
            aperture diameter in mm
        obscuration_ratio: float
            obscuration ratio (between 0 and 1)
        samples : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions.  Order is numpy axis convention, (row, col)

        Returns
        -------
        numpy.ndarray
            binary ndarray representation of the mask
        """

        # aperture radius
        aperture_radius = aperture_diam / 2.0

        # generate the square grid in polar coords
        r, t = _generate_grid(aperture_diam, samples)

        # generate full aperture (circle mask applied to the square grid)
        full_aperture = circle(aperture_radius, r)

        # generate circular centre obscuration
        obscuration = circle(aperture_radius * obscuration_ratio, r)

        # apply full aperture minus obscuration as the mask
        aperture = full_aperture ^ obscuration  # or full_aperture & ~obscuration

        return aperture, r, t


class PrysmOpticsModel:
    """Creates a `prysm` wavefront model from amplitude and phase.

    This essentially creates an optics model (or more precisely,
    a Pupil Function) with a Diffraction plus Wavefront Error Model.

    Parameters
    ----------
    focal_length : float
        effective focal length (or focussing distance) in mm
    amplitude : numpy.ndarray
        array containing the amplitude
    phase : numpy.ndarray, optional
        array containing the optical path error in nm
        if None, assumed zero
    wavelength : float
        wavelength of light with units of microns
    dx : float
        sample spacing with units of mm

    """

    _psf: RichData = None
    _mtf: RichData = None

    @u.wraps(None, (None, u.mm, None, u.nm, u.um, u.mm), False)
    def __init__(self, focal_length, aperture: ndarray, phase, wavelength, dx) -> None:

        # this is monochromatic

        # aperture is an ndarray of True and False
        self.aperture = aperture

        # set focal length
        self.focal_length = focal_length

        # pupil function
        self.pupil = Wavefront.from_amp_and_phase(
            aperture, phase=phase, wavelength=wavelength, dx=dx
        )

    def psf(self, with_units=False) -> RichData:
        """
        Gets the PSF of the Wavefront.

        This is the Point Spread Function of the Wavefront (or pupil function)
        focussed by effective focal length onto the image plane.

        Lazily computes the PSF the first time it is requested.

        `prysm` does not work with units, but a deepcopy of the PSF
        can be generated when the method is called `with_units=True`.

        Parameters
        ----------
        with_units : bool, optional
            a deepcopy of the PSF with units

        Returns
        -------
        RichData
            _description_
        """

        # TODO this is currently hardcoded but shouldn't
        Q = 2

        # lazy init PSF
        if not self._psf:
            # complex field in the plane of the PSF
            coherent_psf = self.pupil.focus(self.focal_length, Q=Q)
            self._psf: RichData = coherent_psf.intensity

        if with_units:
            # psf returned with units
            return _add_units(self._psf)
        else:
            # psf returned without units
            return self._psf

    def mtf(self, with_units=False) -> RichData:

        # lazy init MTF
        if not self._mtf:
            self._mtf = mtf_from_psf(self.psf())

        if with_units:
            # mtf returned with units
            return _add_units(self._mtf)
        else:
            # mtf returned without units
            return self._mtf


# deep copy of the rich data object with units
def _add_units(rich_data: RichData):

    # data is a simple ndarray without units
    data = copy.deepcopy(rich_data.data)

    # inter-sample spacing, mm
    dx = rich_data.dx * u.mm

    # wavelength of light, um
    if rich_data.wavelength:
        wavelength = rich_data.wavelength * u.um
    else:
        wavelength = None

    return RichData(data, dx, wavelength)


def _generate_grid(diameter, samples):
    """Generates a polar grid with the given internal sample points.

    Parameters
    ----------
    diameter : Quantity | float
        centre diameter
    samples : int or tuple of int
        number of samples per dimension.  If a scalar value, broadcast to
        both dimensions.  Order is numpy axis convention, (row, col)

    Returns
    -------
    rho : numpy.ndarray or number
        radial coordinate
    phi : numpy.ndarray or number
        azimuthal coordinate
    """

    # Cartesian grid
    x, y = make_xy_grid(samples, diameter=diameter)
    # radial and azimuthal coords
    r, t = cart_to_polar(x, y)

    return r, t
