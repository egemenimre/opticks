# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from enum import Enum

import numpy as np
from astropy.units import Quantity
from numpy import ndarray
from prysm._richdata import RichData
from prysm.polynomials import sum_of_2d_modes
from prysm.propagation import Wavefront

from opticks import u
from opticks.contrast_model.mtf import MTF_Model_1D, _psf_to_mtf
from opticks.imager_model.aperture import Aperture
from opticks.utils.prysm_utils import OptPathDiff
from opticks.utils.unit_utils import split_value_and_force_unit


class WvlRef(Enum):
    """Reference wavelength selector for PSF output metadata."""

    FIRST = "first"
    LAST = "last"
    MID = "mid"  # middle element of the array
    AVERAGE = "average"

    def resolve(self, wavelengths: Quantity) -> Quantity:
        """Select a wavelength from the Quantity array based on this selector.

        Parameters
        ----------
        wavelengths : Quantity
            wavelengths array to select from

        Returns
        -------
        Quantity
            the selected wavelength
        """
        if self == WvlRef.FIRST:
            return wavelengths[0]  # type: ignore[return-value]
        elif self == WvlRef.LAST:
            return wavelengths[-1]  # type: ignore[return-value]
        elif self == WvlRef.MID:
            return wavelengths[len(wavelengths) // 2]  # type: ignore[return-value]
        elif self == WvlRef.AVERAGE:
            return wavelengths.mean()  # type: ignore[return-value]
        else:
            raise ValueError(f"Unknown WvlRef value: {self}")


class PupilFunction:
    """Pupil function combining amplitude and phase information.

    A pupil function encapsulates one or more wavelength+OPD combinations
    (monochromatic or polychromatic) and can compute the resulting PSF.

    Users should not create PupilFunction objects directly.
    Use `Optics.add_mono_pupil_function()` or
    `Optics.add_poly_pupil_function()` instead.

    Parameters
    ----------
    wavelengths : Quantity
        wavelengths array (in microns)
    opds : list[OptPathDiff | None]
        list of optical path differences (in nm), one per wavelength.
        Use None for zero phase (perfect wavefront).
    aperture : Aperture
        aperture object
    aperture_dx : Quantity
        aperture sample distance (in mm)
    focal_length : Quantity
        focal length (in mm)
    spectral_weights : np.ndarray, optional
        spectral weight of each wavelength, by default uniform
    """

    def __init__(
        self,
        wavelengths: Quantity,
        opds: list[OptPathDiff | None],
        aperture: Aperture,
        aperture_dx: Quantity,
        focal_length: Quantity,
        spectral_weights: np.ndarray | None = None,
    ) -> None:

        # validate lengths
        if len(wavelengths) != len(opds):
            raise ValueError(
                f"wavelengths and opds must have the same length, "
                f"got {len(wavelengths)} and {len(opds)}."
            )

        if spectral_weights is not None and len(spectral_weights) != len(wavelengths):
            raise ValueError(
                f"spectral_weights must have the same length as wavelengths, "
                f"got {len(spectral_weights)} and {len(wavelengths)}."
            )

        self._wavelengths = wavelengths
        self._opds = opds
        self._aperture = aperture
        self._aperture_dx = aperture_dx
        self._focal_length = focal_length
        self._spectral_weights = (
            spectral_weights
            if spectral_weights is not None
            else np.ones(len(wavelengths))
        )

        # build wavefronts eagerly
        self._wavefronts = self._build_wavefronts()

        # PSF and MTF caches
        self._psf: RichData | None = None
        self._mtf: RichData | None = None

    @classmethod
    def monochromatic(
        cls,
        wavelength: Quantity,
        opd: OptPathDiff | None,
        aperture: Aperture,
        aperture_dx: Quantity,
        focal_length: Quantity,
    ) -> "PupilFunction":
        """Create a monochromatic PupilFunction (single wavelength).

        This is useful for sampling of a narrow beam.

        Parameters
        ----------
        wavelength : Quantity
            wavelength of light (in microns)
        opd : OptPathDiff | None
            optical path difference (in nm), or None for zero phase
        aperture : Aperture
            aperture object
        aperture_dx : Quantity
            aperture sample distance (in mm)
        focal_length : Quantity
            focal length (in mm)
        """
        return cls(
            wavelengths=Quantity([wavelength]),
            opds=[opd],
            aperture=aperture,
            aperture_dx=aperture_dx,
            focal_length=focal_length,
        )

    @classmethod
    def polychromatic(
        cls,
        wavelengths: Quantity,
        opds: list[OptPathDiff | None],
        aperture: Aperture,
        aperture_dx: Quantity,
        focal_length: Quantity,
        spectral_weights: np.ndarray | None = None,
    ) -> "PupilFunction":
        """Create a polychromatic PupilFunction (multiple wavelengths).

        This is useful for sampling of a broadband beam.

        Parameters
        ----------
        wavelengths : Quantity
            wavelengths array (in microns)
        opds : list[OptPathDiff | None]
            list of optical path differences (in nm), one per wavelength
        aperture : Aperture
            aperture object
        aperture_dx : Quantity
            aperture sample distance (in mm)
        focal_length : Quantity
            focal length (in mm)
        spectral_weights : np.ndarray, optional
            spectral weight of each wavelength, by default uniform
        """
        if len(wavelengths) == 1:
            return cls.monochromatic(
                wavelengths[0],  # type: ignore[arg-type]
                opds[0],
                aperture,
                aperture_dx,
                focal_length,
            )

        return cls(
            wavelengths=Quantity(wavelengths),
            opds=opds,
            aperture=aperture,
            aperture_dx=aperture_dx,
            focal_length=focal_length,
            spectral_weights=spectral_weights,
        )

    @property
    def is_monochromatic(self) -> bool:
        """True if this pupil function has a single wavelength."""
        return len(self._wavelengths) == 1

    @property
    def num_wavelengths(self) -> int:
        """Number of wavelengths in this pupil function."""
        return len(self._wavelengths)

    @property
    def wavelengths(self) -> Quantity:
        """Wavelengths array."""
        return self._wavelengths

    @property
    def spectral_weights(self) -> np.ndarray:
        """Spectral weights array."""
        return self._spectral_weights

    def compute_psf(
        self,
        psf_dx: Quantity,
        psf_samples: int = 512,
        wvl_ref: WvlRef | None = None,
        with_units: bool = True,
    ) -> RichData:
        """Compute the PSF, cache internally, and return it.

        Parameters
        ----------
        psf_dx : Quantity
            sample distance of the output PSF plane grid (in microns)
        psf_samples : int, optional
            number of samples in the output plane, by default 512
        wvl_ref : WvlRef, optional
            reference wavelength selector for output metadata.
            For monochromatic, defaults to the only wavelength.
            For polychromatic, defaults to WvlRef.AVERAGE.
        with_units : bool, optional
            output the PSF with or without units, by default True

        Returns
        -------
        RichData
            PSF model
        """
        # resolve reference wavelength
        if wvl_ref is None:
            wvl_ref = WvlRef.FIRST if self.is_monochromatic else WvlRef.AVERAGE
        ref_wvl = wvl_ref.resolve(self._wavelengths)

        # compute PSF
        psf = _compute_psf(
            self._wavefronts,
            self._focal_length,
            ref_wvl,
            psf_dx,
            psf_samples,
            self._spectral_weights,
        )

        # add units if requested
        if with_units:
            psf = RichData(psf.data, psf.dx * u.um, psf.wavelength * u.um)

        # cache and invalidate MTF cache
        self._psf = psf
        self._mtf = None

        return psf

    @property
    def psf(self) -> RichData:
        """Return the cached PSF.

        Raises
        ------
        ValueError
            if compute_psf() has not been called yet
        """
        if self._psf is None:
            raise ValueError("PSF has not been computed yet. Call compute_psf() first.")
        return self._psf

    @property
    def mtf(self) -> RichData:
        """Return the MTF, computing lazily from the cached PSF if needed.

        Raises
        ------
        ValueError
            if compute_psf() has not been called yet
        """
        if self._psf is None:
            raise ValueError("PSF has not been computed yet. Call compute_psf() first.")

        if self._mtf is None:
            self._mtf = _psf_to_mtf(self._psf, with_units=True)

        return self._mtf

    def to_MTF_Model_1D(self, slice: str) -> MTF_Model_1D:
        """Convert the cached 2D MTF to a 1D MTF model.

        Extracts a 1D slice from the 2D MTF and returns it as an
        ``MTF_Model_1D`` object. The PSF must have been computed
        before calling this method.

        Possible slice strings are ``x``, ``y``, ``azavg``, ``azavmedian``,
        ``azmin``, ``azpv``, ``azvar``, ``azstd``.

        Parameters
        ----------
        slice : str
            slice type (e.g., "x", "y", "azavg")

        Returns
        -------
        MTF_Model_1D
            1D MTF model

        Raises
        ------
        ValueError
            if compute_psf() has not been called yet
        """
        return MTF_Model_1D.from_mtf_2d(self.mtf, slice)

    def _build_wavefronts(self) -> list[Wavefront]:
        """Build prysm Wavefront objects from stored wavelengths/OPDs."""

        wavefronts = []
        for wvl, opd in zip(self._wavelengths, self._opds, strict=True):
            opd_data = opd.strip_units(u.nm).data if opd else None
            wf = Wavefront.from_amp_and_phase(
                self._aperture.data,
                phase=opd_data,
                wavelength=wvl.to_value(u.um),
                dx=self._aperture_dx.to_value(u.mm),
            )
            wavefronts.append(wf)
        return wavefronts


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
        PSF model (without units)
    """

    focal_length_val, _ = split_value_and_force_unit(focal_length, u.mm)
    wvl_val, _ = split_value_and_force_unit(wvl, u.um)
    psf_dx_val, _ = split_value_and_force_unit(psf_dx, u.um)

    psf_components = []

    # focus all WF objects and compute the monochromatic PSF
    for pupil in pupils:
        # complex field in the plane of the PSF (no unit support)
        # Note: focusing changes aperture sampling,
        # and is a function of wavelength. Therefore we use fixed sampling
        # for multiple wavelengths.
        coherent_psf = pupil.focus_fixed_sampling(
            focal_length_val, psf_dx_val, psf_samples
        )
        psf_data = coherent_psf.intensity.data
        # sum of intensities, wvls are incoherent to each other
        psf_components.append(psf_data)

    # Create psf array via summation (no unit support)
    psf_data = sum_of_2d_modes(np.asarray(psf_components), spectral_weights)

    # Add scaling and wavelength information
    # pupil to psf plane means dx is switched from mm to um
    psf = RichData(psf_data, psf_dx_val, wvl_val)

    return psf
