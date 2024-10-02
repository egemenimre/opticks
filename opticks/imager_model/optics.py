# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.


import numpy as np
from numpy import ndarray
from pint import Quantity
from prysm._richdata import RichData
from prysm.coordinates import cart_to_polar, make_xy_grid
from prysm.geometry import circle
from prysm.polynomials import ansi_j_to_nm, sum_of_2d_modes, zernike_nm_sequence
from prysm.propagation import Wavefront
from strictyaml import Map, Str

from opticks import u
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.prysm_utils import richdata_with_units
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
    @u.wraps(None, (None, u.mm, None, None), False)
    def circle_aperture(
        cls,
        aperture_diam: Quantity | float,
        samples,
        with_units=False,
    ) -> tuple[ndarray, ndarray, ndarray]:
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
        with_units : bool, optional
            `r` and `t` returned with units (in mm and radians, respectively)

        Returns
        -------
        aperture: numpy.ndarray
            binary ndarray representation of the mask
        rho : numpy.ndarray or number
            radial coordinate
        phi : numpy.ndarray or number
            azimuthal coordinate
        """

        # aperture radius
        aperture_radius = aperture_diam / 2.0

        # generate the square grid in polar coords
        r, t = _generate_polar_grid(aperture_diam, samples)

        # generate aperture (circle mask applied to the square grid)
        if with_units:
            return circle(aperture_radius, r), r * u.mm, t * u.rad
        else:
            return circle(aperture_radius, r), r, t

    @classmethod
    @u.wraps(None, (None, u.mm, None, None, None), False)
    def circle_aperture_with_obscuration(
        cls,
        aperture_diam: Quantity | float,
        obscuration_ratio: float,
        samples,
        with_units=False,
    ) -> tuple[ndarray, ndarray, ndarray]:
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
        with_units : bool, optional
            `r` and `t` returned with units (in mm and radians, respectively)

        Returns
        -------
        aperture: numpy.ndarray of booleans
            binary ndarray representation of the mask
        rho : numpy.ndarray or number
            radial coordinate
        phi : numpy.ndarray or number
            azimuthal coordinate
        """

        # aperture radius
        aperture_radius = aperture_diam / 2.0

        # generate the square grid in polar coords
        r, t = _generate_polar_grid(aperture_diam, samples)

        # generate full aperture (circle mask applied to the square grid)
        full_aperture = circle(aperture_radius, r)

        # generate circular centre obscuration
        obscuration = circle(aperture_radius * obscuration_ratio, r)

        # apply full aperture minus obscuration as the mask
        aperture = full_aperture ^ obscuration  # or full_aperture & ~obscuration

        # can also use floats of 0 and 1 instead of bool
        # this also enables varying the amount of light through the aperture
        # aperture = aperture.astype(float)

        if with_units:
            return aperture, r * u.mm, t * u.rad
        else:
            return aperture, r, t


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

    aperture_dx = None
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

        # use units for the r and t values
        with_units = True

        if aperture is None:
            # aperture function not defined, use circular aperture
            # r and t (with units) currently ignored
            self.aperture, r, t = ApertureFactory.circle_aperture(
                self.params.aperture_diameter, samples, with_units
            )
        else:
            # aperture function defined, use it
            self.aperture = aperture

        # reset the sample size to the actual sample size
        samples = len(self.aperture)

        # compute aperture sample distance
        self.aperture_dx = (self.params.aperture_diameter / samples).to(u.mm)

    @u.wraps(None, (None, u.um, u.nm), False)
    def add_pupil_func(self, wavelength, opd=None):
        """Adds a pupil function.

        A pupil function combines amplitude with phase information
        to form a wavefront.

        Note that the pupil function is valid for a single point
        on the image frame and for a single wavelength (monochromatic).
        One such functions for each wavelength should be evaluated and
        then properly summed to yield the polychromatic PSF.

        Parameters
        ----------
        wavelength : float | Quantity
            wavelength of light with units of microns
        opd : numpy.ndarray | Quantity, optional
            array containing the optical path error in nm
            if None, assumed zero
        """

        # Generate the pupil function
        pupil = Wavefront.from_amp_and_phase(
            self.aperture,  # amplitudes
            phase=opd,  # float ndarray in nm with u.wraps
            wavelength=wavelength,  # float in um with u.wraps
            dx=self.aperture_dx.m_as(u.mm),
        )

        # add to the list of pupil functions
        self.pupils.append(pupil)

    def psf(
        self,
        wvl_ref: float | Quantity,
        psf_dx: float | Quantity,
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
        wvl_ref : float | Quantity
            reference wavelength (in microns)
        psf_dx : float | Quantity
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

        # compute the PSF
        psf = _compute_psf(
            self.pupils, focal_length, wvl_ref, psf_dx, psf_samples, spectral_weights
        )

        # return PSF with or without units
        if with_units:
            # pupil to PSF plane means dx is switched from mm to um
            return richdata_with_units(psf, dx_units=u.um)
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


@u.wraps(u.nm, (u.nm, u.mm, u.mm, u.rad), False)
def zernike_opd(
    wfe_rms: list, aperture_diameter: Quantity | float, r: ndarray, t: ndarray
) -> ndarray:
    """Computes the Optical Path Difference via Zernike Polynomials.

    Generates the Zernike Polynomials to the order corresponding
    to the number of coefficients (e.g., 9 coefficients = mode 8)
    sums them properly, adding the WFE RMS Zernike coefficients.

    The result is the monochromatic OPD for a single location
    on the PSF plane.

    Parameters
    ----------
    wfe_rms : list
        ANSI list aberration coefficients (WFE RMS) (in nm)
    diameter : Quantity | float
        aperture diameter in mm
    r : ndarray
        radial coordinate of the aperture grid (in mm)
    t : ndarray
        azimuthal coordinate of the aperture grid (in rad)

    Returns
    -------
    ndarray
        Optical Path Difference (OPD)
    """

    # mode n = (n+1) elements
    elems = len(wfe_rms)

    # Generate the (n,m) tuples in ANSI order
    nms = [ansi_j_to_nm(i) for i in range(0, elems)]

    # radial coords normalised by aperture radius
    # normalisation required by the polynomials
    ap_radius = aperture_diameter / 2.0
    rho = r / ap_radius

    # compute the polynomials (dimensionless)
    mode = list(zernike_nm_sequence(nms, rho, t))

    # monochromatic OPD with multiple aberrations
    # units are in nm as the wfe_rms is forced to be in nm (via wraps)
    opd = sum_of_2d_modes(mode, wfe_rms) * u.nm

    return opd


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


def _generate_polar_grid(diameter, samples) -> tuple[ndarray, ndarray]:
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
    rho, phi: (ndarray, ndarray)
        radial coordinate and azimuthal coordinate
    """

    # Cartesian grid
    x, y = make_xy_grid(samples, diameter=diameter)
    # radial and azimuthal coords
    r, t = cart_to_polar(x, y)

    return r, t
