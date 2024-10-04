# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


import numpy as np
from numpy import ndarray
from pint import Quantity
from prysm._richdata import RichData
from prysm.geometry import circle
from prysm.polynomials import sum_of_2d_modes
from prysm.propagation import Wavefront
from strictyaml import Map, Str

from opticks import u
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.prysm_utils import Grid, OptPathDiff
from opticks.utils.yaml_helpers import Qty

optics_schema = {
    "name": Str(),
    "focal_length": Qty(),
    "aperture_diameter": Qty(),
    "image_diam_on_focal_plane": Qty(),
}
"""Schema containing optical parameters."""


class ApertureFactory:
    """Generates oft used aperture models for `prysm`."""

    @classmethod
    def circle_aperture(
        cls,
        aperture_diam: Quantity,
        samples,
        with_units=True,
    ) -> tuple[ndarray, Grid]:
        """Circle aperture model for `prysm`.

        This is valid for many if not most refractive optics.

        The output is an `ndarray` of bools, therefore the unit input
        is not important for the aperture generation.

        Parameters
        ----------
        aperture_diam : Quantity
            aperture diameter in mm
        samples : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions.  Order is numpy axis convention, (row, col)
        with_units : bool, optional
            Grid returned with units (in mm and radians)

        Returns
        -------
        aperture: numpy.ndarray
            ndarray representation of the mask
        grid : Grid
            Grid object associated with the aperture
        """

        # aperture radius
        aperture_radius = aperture_diam / 2.0

        # generate the square grid in polar coords
        grid = Grid.from_size(samples, aperture_diam)
        r, t = grid.polar()

        # generate aperture (circle mask applied to the square grid)
        aperture = circle(aperture_radius, r)

        # strip units if needed
        if not with_units:
            grid = grid.strip_units(u.mm)

        return aperture, grid

    @classmethod
    def circle_aperture_with_obscuration(
        cls,
        aperture_diam: Quantity,
        obscuration_ratio: float,
        samples,
        with_units=True,
    ) -> tuple[ndarray, Grid]:
        """Circle aperture model for `prysm` with circular centre obscuration.

        This is valid for many reflective telescopes. The obscuration
        is usually the secondary mirror.

        The output is an `ndarray` of bools, therefore the unit input
        is not important for the aperture generation.

        Parameters
        ----------
        aperture_diam : Quantity
            aperture diameter in mm
        obscuration_ratio: float
            obscuration ratio (between 0 and 1)
        samples : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions.  Order is numpy axis convention, (row, col)
        with_units : bool, optional
            Grid returned with units (in mm and radians)

        Returns
        -------
        aperture: numpy.ndarray of booleans
            ndarray representation of the mask
        grid : Grid
            Grid object associated with the aperture
        """

        # aperture radius
        aperture_radius = aperture_diam / 2.0

        # generate the square grid in polar coords
        grid = Grid.from_size(samples, aperture_diam)
        r, t = grid.polar()

        # generate full aperture (circle mask applied to the square grid)
        full_aperture = circle(aperture_radius, r)

        # generate circular centre obscuration
        obscuration = circle(aperture_radius * obscuration_ratio, r)

        # apply full aperture minus obscuration as the mask
        aperture = full_aperture ^ obscuration  # or full_aperture & ~obscuration

        # can also use floats of 0 and 1 instead of bool
        # this also enables varying the amount of light through the aperture
        # aperture = aperture.astype(float)

        # strip units if needed
        if not with_units:
            grid = grid.strip_units(u.mm)

        return aperture, grid


class Optics(ImagerComponent):
    """
    Class containing generic Optics parameters.
    """

    # Empty list of pupil functions
    pupils: list[Wavefront] = []

    @classmethod
    def schema(cls) -> Map:
        return Map(optics_schema)

    @classmethod
    def _params_class_name(cls) -> str:
        return "OpticsParams"

    # ---------- begin modelling functions ----------

    aperture_dx: Quantity = None
    """Aperture sample distance (in mm).
    """

    def set_aperture_model(self, aperture: ndarray = None, samples: int = 400) -> None:
        """Sets the aperture model for the optics.

        The aperture model is as defined by `prysm`, it is an `ndarray` grid of
        `True` and `False` data, where `True` allows light through. It can also be
        an `ndarray` grid of 0 and 1 (and possibly anything in between).

        If `None`, a new circular aperture is defined using the aperture diameter
        of the optics and the `samples` parameter (for one side of the square grid).
        If an aperture is defined, `samples` is ignored.

        Parameters
        ----------
        aperture : ndarray, optional
            aperture grid
        samples : int, optional
            sample size of the default circular aperture,
            ignored if aperture is user defined
        """

        # use units for the grid (currently not used)
        with_units = True

        if aperture is None:
            # aperture function not defined, use circular aperture
            # grid currently ignored
            self.aperture, grid = ApertureFactory.circle_aperture(
                self.params.aperture_diameter, samples, with_units
            )
        else:
            # aperture function defined, use it
            self.aperture = aperture

        # reset the sample size to the actual sample size
        samples = len(self.aperture)

        # compute aperture sample distance
        self.aperture_dx = (self.params.aperture_diameter / samples).to(u.mm)

    def add_pupil_func(self, wavelength, opd: OptPathDiff = None):
        """Adds a pupil function.

        A pupil function combines amplitude with phase information
        to form a wavefront.

        Note that the pupil function is valid for a single point
        on the image frame and for a single wavelength (monochromatic).
        One such functions for each wavelength should be evaluated and
        then properly summed to yield the polychromatic PSF.

        Parameters
        ----------
        wavelength : Quantity
            wavelength of light with units of microns
        opd : OptPathDiff, optional
            array containing the optical path error in nm
            if None, assumed zero
        """

        # strip units for the Wavefront as it cannot handle units well
        opd_data = opd.strip_units(u.nm).data

        # Generate the pupil function (no units)
        pupil = Wavefront.from_amp_and_phase(
            self.aperture,  # amplitudes
            phase=opd_data,  # float ndarray in nm
            wavelength=wavelength.m_as(u.um),
            dx=self.aperture_dx.m_as(u.mm),
        )

        # add to the list of pupil functions
        self.pupils.append(pupil)

    def psf(
        self,
        wvl_ref: Quantity,
        psf_dx: Quantity,
        psf_samples: int = 512,
        spectral_weights: np.ndarray = None,
        with_units=True,
    ) -> RichData:
        """Computes the PSF for a single point on the Image Plane.

        The function can handle monochromatic or polychromatic PSF
        computations. The PSF or Image Plane is resampled to
        the user defined grid.

        The spectral weights array should have the same number of elements
        as the number of Pupil Functions. If no weighting is defined,
        a uniform weighting is assumed.

        The operation can be fairly expensive, therefore it is advised
        to keep the output PSF stored.

        Parameters
        ----------
        wvl_ref : Quantity
            reference wavelength (in microns)
        psf_dx : Quantity
            sample distance of the output PSF Plane grid (in microns)
        psf_samples : int, optional
            number of samples in the output plane.
            If int, interpreted as square else interpreted as (x,y),
            which is the reverse of numpy's (y, x) row major ordering.
            By default 512
        spectral_weights : np.ndarray, optional
            spectral weight of each wavelength, by default None
        with_units : bool, optional
            output the PSF with or without units, by default True

        Returns
        -------
        RichData
            PSF model with or without units
        """

        focal_length = self.params.focal_length

        if not spectral_weights:
            # set uniform weights if no weights defined
            spectral_weights = np.ones_like(self.pupils)

        # compute the PSF (no units via wraps)
        psf = _compute_psf(
            self.pupils, focal_length, wvl_ref, psf_dx, psf_samples, spectral_weights
        )

        # return the PSF with or without units
        if with_units:
            # pupil to PSF plane means dx is switched from mm to um
            return RichData(psf.data, psf.dx * u.um, psf.wavelength * u.um)
        else:
            return psf

    @property
    def f_number(self) -> float:
        """
        F-number.

        Computed as: focal length / aperture diameter

        """
        return (
            (self.params.focal_length / self.params.aperture_diameter)
            .to_reduced_units()
            .m
        )

    @property
    def full_optical_fov(self) -> Quantity:
        """
        Full optical Field of View.

        Actual FoV depends on the detector size, but cannot be wider than this value.

        """
        return 2 * np.arctan(
            (self.params.image_diam_on_focal_plane / 2.0) / self.params.focal_length
        ).to(u.deg)

    @property
    def aperture_area(self) -> Quantity:
        """
        Aperture area.

        Computed as: pi x (aperture diameter/2 )^2
        """
        return np.pi * (self.params.aperture_diameter / 2.0) ** 2

    @property
    def aperture_solid_angle(self) -> Quantity:
        r"""
        Aperture solid angle in steradians.
        """
        return (
            np.pi
            / (self.params.focal_length / (self.params.aperture_diameter / 2.0)) ** 2
            * u.steradian
        )

    @u.check(None, "[length]")
    def spatial_cutoff_freq(self, ref_wavelength: Quantity) -> Quantity:
        """
        Spatial cutoff frequency, assumes perfect incoherent optics.

        Determines the theoretical limit of the optical resolution, or the smallest
        object resolvable by the optical system.

        Computed as: `1/(wavelength x F#)` in cycles per mm

        Parameters
        ----------
        ref_wavelength : Quantity
            Reference wavelength

        Returns
        -------
        Quantity
            Spatial cutoff frequency (in cycles/mm)
        """
        # perfect incoherent optics
        return (1.0 * u.cy) / (ref_wavelength * self.f_number).to(u.mm)


@u.wraps(None, (None, u.mm, u.um, u.um, None, None), False)
def _compute_psf(
    pupils: list[Wavefront],
    focal_length: float | Quantity,
    wvl: float | Quantity,
    psf_dx: float | Quantity,
    psf_samples: int,
    spectral_weights: ndarray,
) -> RichData:
    """Computes the PSF for a single point on the Image Plane.

    The function can handle monochromatic or polychromatic PSF
    computations. The PSF or Image Plane is resampled to
    the user defined grid.

    The spectral weights array should have the same number of elements
    as the number of Pupil Functions.

    The operation can be fairly expensive, therefore it is advised
    to keep the output PSF stored.

    Parameters
    ----------
    pupils : list[Wavefront]
        list of Pupil functions or Wavefronts
    focal_length : float | Quantity
        focal length in mm
    wvl : float | Quantity
        reference wavelength (in microns)
    psf_dx : float | Quantity
        sample distance of the output PSF Plane grid (in microns)
    psf_samples : int
        number of samples in the output plane.
        If int, interpreted as square else interpreted as (x,y),
        which is the reverse of numpy's (y, x) row major ordering.
    spectral_weights : np.ndarray
        spectral weight of each wavelength

    Returns
    -------
    RichData
        PSF model
    """

    psf_components = []

    # focus all WF objects and compute the monochromatic PSF
    for pupil in pupils:
        # complex field in the plane of the PSF (no unit support)
        # Note: focusing changes aperture sampling,
        # and is a function of wavelength. Therefore we use fixed sampling
        # for multiple wavelengths.
        coherent_psf = pupil.focus_fixed_sampling(focal_length, psf_dx, psf_samples)
        psf_data = coherent_psf.intensity.data
        # sum of intensities, wvls are incoherent to each other
        psf_components.append(psf_data)

    # Create psf array via summation (no unit support)
    psf_data = sum_of_2d_modes(psf_components, spectral_weights)

    # Add scaling and wavelength information
    # pupil to psf plane means dx is switched from mm to um
    psf = RichData(psf_data, psf_dx, wvl)

    return psf
