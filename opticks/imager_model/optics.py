# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import numpy as np
from astropy.units import Quantity
from prysm._richdata import RichData
from pydantic import Field

from opticks import u
from opticks.contrast_model.mtf import MTF_Model_1D
from opticks.contrast_model.optics_mtf import FieldAberrationModel
from opticks.imager_model.aperture import Aperture
from opticks.imager_model.imager_component import ImagerComponent
from opticks.imager_model.pupil import PupilFunction, WvlRef
from opticks.utils.parser_helpers import PositivePydanticQty
from opticks.utils.prysm_utils import OptPathDiff


class Optics(ImagerComponent):
    """
    Class containing generic Optics parameters.
    """

    name: str
    focal_length: PositivePydanticQty
    aperture_diameter: PositivePydanticQty
    image_diam_on_focal_plane: PositivePydanticQty

    # Non-serialized state (set after init, not from YAML)
    pupil_functs: dict[str, PupilFunction] = Field(default_factory=dict, exclude=True)
    aperture_dx: Quantity | None = Field(default=None, exclude=True)
    aperture: Aperture | None = Field(default=None, exclude=True)

    # ---------- begin modelling functions ----------

    def set_aperture_model(
        self, aperture: Aperture | None = None, samples: int = 400
    ) -> None:
        """Sets the aperture model for the optics.

        The aperture model is as defined by `prysm`, it contains an `ndarray` grid of
        `True` and `False` data, where `True` allows light through. It can also be
        an `ndarray` grid of 0 and 1 (and anything in between).

        If `None`, a new circular aperture is defined using the aperture diameter
        of the optics and the `samples` parameter (for one side of the square grid).
        If an aperture is defined, `samples` is ignored.

        Grid is optional and is added only when a new aperture is defined

        Parameters
        ----------
        aperture : aperture
            aperture object
        samples : int, optional
            sample size of the default circular aperture,
            ignored if aperture is user defined
        """

        # use units for the grid (currently not used)
        with_units = True

        if aperture is None:
            # aperture function not defined, use circular aperture
            # grid currently ignored
            self.aperture = Aperture.circle_aperture(
                self.aperture_diameter, samples, with_units
            )
        else:
            # aperture function defined by the User, set it
            self.aperture = aperture

        # reset the sample size to the actual sample size
        samples = len(self.aperture.data)

        # compute aperture sample distance
        self.aperture_dx = (self.aperture_diameter / samples).to(u.mm)

    def add_mono_pupil_function(
        self,
        name: str,
        wavelength: Quantity,
        opd: OptPathDiff | None = None,
    ) -> PupilFunction:
        """Create a monochromatic PupilFunction and register it.

        Parameters
        ----------
        name : str
            identifier for this pupil function
        wavelength : Quantity
            wavelength of light (in microns)
        opd : OptPathDiff, optional
            optical path difference (in nm), or None for zero phase

        Returns
        -------
        PupilFunction
            the created PupilFunction for convenience
        """
        if self.aperture is None or self.aperture_dx is None:
            raise ValueError(
                "Aperture model must be set before adding pupil functions. "
                "Call set_aperture_model() first."
            )

        pf = PupilFunction.monochromatic(
            wavelength,
            opd,
            self.aperture,
            self.aperture_dx,
            self.focal_length,
        )
        self.pupil_functs[name] = pf
        return pf

    def add_poly_pupil_function(
        self,
        name: str,
        wavelengths: Quantity,
        opds: list[OptPathDiff | None],
        spectral_weights: np.ndarray | None = None,
    ) -> PupilFunction:
        """Create a polychromatic PupilFunction and register it.

        Parameters
        ----------
        name : str
            identifier for this pupil function
        wavelengths : Quantity
            wavelengths array (in microns)
        opds : list[OptPathDiff | None]
            list of optical path differences (in nm), one per wavelength
        spectral_weights : np.ndarray, optional
            spectral weight of each wavelength, by default uniform

        Returns
        -------
        PupilFunction
            the created PupilFunction for convenience
        """
        if self.aperture is None or self.aperture_dx is None:
            raise ValueError(
                "Aperture model must be set before adding pupil functions. "
                "Call set_aperture_model() first."
            )

        pf = PupilFunction.polychromatic(
            wavelengths,
            opds,
            self.aperture,
            self.aperture_dx,
            self.focal_length,
            spectral_weights,
        )
        self.pupil_functs[name] = pf
        return pf

    def remove_pupil_function(self, name: str) -> None:
        """Remove a named PupilFunction.

        Parameters
        ----------
        name : str
            identifier of the pupil function to remove
        """
        del self.pupil_functs[name]

    def clear_pupil_functs(self) -> None:
        """Remove all PupilFunctions."""
        self.pupil_functs.clear()

    def get_pupil_function(self, name: str) -> PupilFunction:
        """Retrieve a named PupilFunction.

        Parameters
        ----------
        name : str
            identifier of the pupil function

        Returns
        -------
        PupilFunction
            the requested PupilFunction
        """
        return self.pupil_functs[name]

    def compute_psf(
        self,
        pupil_name: str,
        psf_dx: Quantity,
        psf_samples: int = 512,
        wvl_ref: WvlRef | None = None,
        with_units: bool = True,
    ) -> RichData:
        """Compute PSF for the named PupilFunction (caches + returns).

        Parameters
        ----------
        pupil_name : str
            identifier of the pupil function
        psf_dx : Quantity
            sample distance of the output PSF plane grid (in microns)
        psf_samples : int, optional
            number of samples in the output plane, by default 512
        wvl_ref : WvlRef, optional
            reference wavelength selector for output metadata
        with_units : bool, optional
            output the PSF with or without units, by default True

        Returns
        -------
        RichData
            PSF model
        """
        return self.pupil_functs[pupil_name].compute_psf(
            psf_dx, psf_samples, wvl_ref, with_units
        )

    def psf(self, pupil_name: str) -> RichData:
        """Return cached PSF for the named PupilFunction.

        Parameters
        ----------
        pupil_name : str
            identifier of the pupil function

        Returns
        -------
        RichData
            cached PSF model

        Raises
        ------
        ValueError
            if compute_psf() has not been called yet
        """
        pf = self.pupil_functs[pupil_name]

        if pf._psf is None:
            raise ValueError(
                f"PSF has not been computed yet for {pupil_name}."
                "Call compute_psf() with the desired parameters first."
            )
        return pf.psf

    def mtf(self, pupil_name: str) -> RichData:
        """Return MTF derived from cached PSF for the named PupilFunction.

        Lazily computed on first access; cached until PSF is recomputed.

        Parameters
        ----------
        pupil_name : str
            identifier of the pupil function

        Returns
        -------
        RichData
            MTF model

        Raises
        ------
        ValueError
            if compute_psf() has not been called yet
        """
        pf = self.pupil_functs[pupil_name]

        if pf._psf is None:
            raise ValueError(
                f"PSF has not been computed yet for {pupil_name}. "
                "Call compute_psf() with the desired parameters first and then call mtf()."
            )
        return pf.mtf

    def to_MTF_Model_1D(self, pupil_name: str, slice: str) -> MTF_Model_1D:
        """Convert the cached 2D MTF to a 1D MTF model for the named PupilFunction.

        Parameters
        ----------
        pupil_name : str
            identifier of the pupil function
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
        pf = self.pupil_functs[pupil_name]

        if pf._psf is None:
            raise ValueError(
                f"PSF has not been computed yet for {pupil_name}. "
                "Call compute_psf() with the desired parameters first and then call mtf()."
            )

        return pf.to_MTF_Model_1D(slice)

    @u.quantity_input(wavelength="length")
    def field_mtf_model_1d(
        self,
        field_model: FieldAberrationModel,
        h: float,
        wavelength: Quantity,
    ) -> MTF_Model_1D:
        """Aberrated 1D MTF at normalised field position *h* (Tier B path).

        Uses the Seidel model to compute the total RMS WFE at the given
        field position and feeds it into the empirical ATF model.

        Parameters
        ----------
        field_model : FieldAberrationModel
            Seidel coefficient set
        h : float
            Normalised radial field coordinate (0 = on-axis, 1 = edge)
        wavelength : Quantity["length"]
            Wavelength for the MTF computation

        Returns
        -------
        MTF_Model_1D
            Aberrated 1D MTF model at the given field position
        """
        w_rms_waves = field_model.w_rms_waves(h, wavelength)
        spatial_cutoff = self.spatial_cutoff_freq(wavelength)
        return MTF_Model_1D.emp_model_aberrated_optics(
            spatial_cutoff,
            w_rms=w_rms_waves,
            wavelength=wavelength,
        )

    @property
    def f_number(self) -> float:
        """
        F-number.

        Computed as: focal length / aperture diameter

        """
        return (self.focal_length / self.aperture_diameter).decompose().value

    @property
    def full_optical_fov(self) -> Quantity:
        """
        Full optical Field of View.

        Actual FoV depends on the detector size, but cannot be wider than this value.

        """
        return 2 * np.arctan(
            (self.image_diam_on_focal_plane / 2.0) / self.focal_length
        ).to(u.deg, copy=False)

    @property
    def aperture_area(self) -> Quantity:
        """
        Aperture area.

        Computed as: pi x (aperture diameter/2 )^2
        """
        return np.pi * (self.aperture_diameter / 2.0) ** 2

    @property
    def aperture_solid_angle(self) -> Quantity:
        r"""
        Aperture solid angle in steradians.
        """
        return (
            np.pi
            / (self.focal_length / (self.aperture_diameter / 2.0)) ** 2
            * u.steradian
        )

    @u.quantity_input(ref_wavelength="length")
    def spatial_cutoff_freq(self, ref_wavelength: Quantity) -> Quantity:
        """
        Spatial cut-off frequency, assumes perfect incoherent optics.

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
        return (1.0 * u.cy) / (ref_wavelength * self.f_number).to(u.mm, copy=False)  # type: ignore[union-attr]
