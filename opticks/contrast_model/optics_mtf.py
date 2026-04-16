# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Optics-level MTF models (ideal diffraction, aberration transfer factor,
field-dependent Seidel model, PSF-to-MTF conversion).

Mirrors the pattern of ``detector_mtf.py`` — houses the backend
functions that ``MTF_Model_1D`` classmethods in ``mtf.py`` delegate to.
"""

import numpy as np
from astropy.units import Quantity
from numpy.typing import NDArray
from prysm._richdata import RichData
from prysm.otf import mtf_from_psf

from opticks import u
from opticks.utils.prysm_utils import richdata_with_units

# ---------------------------------------------------------------------------
# Backend helpers (moved from mtf.py)
# ---------------------------------------------------------------------------


def _ideal_optical_mtf(
    input_line_freq: Quantity, spatial_cutoff_freq: Quantity
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

    with u.set_enabled_equivalencies(u.dimensionless_angles()):
        mtf_value = 2 / np.pi * (np.arccos(nu) - nu * np.sqrt(1 - nu**2))

    # force return float
    return _force_return_float(mtf_value)


def _aberrated_optical_mtf(
    input_line_freq: Quantity,
    spatial_cutoff_freq: Quantity,
    w_rms: float | NDArray[np.float64],
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
    w_rms : float or NDArray[np.float64]
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

    mtf_value = mtf_ideal * atf

    # force return float
    return _force_return_float(mtf_value)


def _aberration_transfer_factor(
    input_line_freq: Quantity,
    spatial_cutoff_freq: Quantity,
    w_rms: float | NDArray[np.float64],
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
    w_rms : float or NDArray[np.float64]
        RMS of the total wavefront error (in wavelengths)

    Returns
    -------
    float | NDArray[np.float64]
        ATF value (<1)
    """

    # normalised optical frequency
    nu = input_line_freq / spatial_cutoff_freq

    # compute ATF
    mtf_value = 1 - ((w_rms / 0.18) ** 2 * (1 - 4 * (nu - 0.5) ** 2))

    # force return float
    return _force_return_float(mtf_value)


def _psf_to_mtf(psf: RichData, with_units=False) -> RichData:
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
        psf_copy = RichData(psf.data, psf.dx.to_value(u.um), None)

        mtf = mtf_from_psf(psf_copy)

        if psf.wavelength is not None:
            mtf.wavelength = psf.wavelength.to_value(u.um)
    else:
        # psf does not have units, create mtf directly
        mtf = mtf_from_psf(psf)

        if psf.wavelength is not None:
            mtf.wavelength = psf.wavelength

    # the resulting mtf is guaranteed to be without units
    # output dx in cy/mm

    if with_units:
        # add units to MTF
        return richdata_with_units(mtf, dx_units=u.cy / u.mm)
    else:
        # mtf returned without units
        return mtf


def _force_return_float(mtf_value):
    if isinstance(mtf_value, Quantity):
        return mtf_value.decompose().value
    else:
        return mtf_value


# ---------------------------------------------------------------------------
# Seidel / Hopkins field-dependent aberration model
# ---------------------------------------------------------------------------

# Mahajan RMS-per-term variance coefficients (Table 7.3, Aberration Theory
# Made Simple, 2nd Ed., SPIE Press, 2011).  Each maps a peak Seidel
# coefficient to its RMS contribution via  sigma_i = c_i * |W_ijk|.
_RMS_COEFF_SPHERICAL = 1.0 / (6.0 * np.sqrt(5.0))  # W_040
_RMS_COEFF_COMA = 1.0 / (2.0 * np.sqrt(2.0) * np.sqrt(2.0))  # W_131
_RMS_COEFF_ASTIGMATISM = 1.0 / (2.0 * np.sqrt(6.0))  # W_222
_RMS_COEFF_FIELD_CURVATURE = 1.0 / (2.0 * np.sqrt(3.0))  # W_220


class FieldAberrationModel:
    """Third-order (Seidel) wave-aberration model with field dependence.

    Stores the four MTF-affecting Seidel coefficients evaluated at the
    edge of the field (H = 1) as astropy length Quantities (nm, µm, …).
    Distortion (W_311) is omitted because it does not affect the MTF.

    The field parameter *h* is the normalised radial distance from the
    optical axis: 0 at centre, 1 at the edge of the field.

    Parameters
    ----------
    w040 : Quantity["length"]
        Spherical aberration (field-independent peak coefficient)
    w131_edge : Quantity["length"]
        Coma peak coefficient at H = 1
    w222_edge : Quantity["length"]
        Astigmatism peak coefficient at H = 1
    w220_edge : Quantity["length"]
        Field curvature / Petzval peak coefficient at H = 1
    reference_wavelength : Quantity["length"]
        Wavelength used for "in-waves" conversions
    """

    def __init__(
        self,
        w040: Quantity,
        w131_edge: Quantity,
        w222_edge: Quantity,
        w220_edge: Quantity,
        reference_wavelength: Quantity,
    ) -> None:
        self.w040: Quantity = w040.to(u.nm)
        self.w131_edge: Quantity = w131_edge.to(u.nm)
        self.w222_edge: Quantity = w222_edge.to(u.nm)
        self.w220_edge: Quantity = w220_edge.to(u.nm)
        self.reference_wavelength: Quantity = reference_wavelength.to(u.nm)

    # ---- constructors ----

    @classmethod
    def from_waves(
        cls,
        w040: float,
        w131_edge: float,
        w222_edge: float,
        w220_edge: float,
        reference_wavelength: Quantity,
    ) -> "FieldAberrationModel":
        """Create from peak coefficients expressed in waves.

        Parameters
        ----------
        w040 : float
            Spherical aberration in waves (e.g. 0.05 = λ/20)
        w131_edge : float
            Coma at H = 1 in waves
        w222_edge : float
            Astigmatism at H = 1 in waves
        w220_edge : float
            Field curvature at H = 1 in waves
        reference_wavelength : Quantity["length"]
            Wavelength that defines "one wave"

        Returns
        -------
        FieldAberrationModel
        """
        ref = reference_wavelength.to(u.nm)
        return cls(
            w040=w040 * ref,
            w131_edge=w131_edge * ref,
            w222_edge=w222_edge * ref,
            w220_edge=w220_edge * ref,
            reference_wavelength=reference_wavelength,
        )

    @classmethod
    def from_quantity(
        cls,
        w040: Quantity,
        w131_edge: Quantity,
        w222_edge: Quantity,
        w220_edge: Quantity,
        reference_wavelength: Quantity,
    ) -> "FieldAberrationModel":
        """Create from peak coefficients expressed as physical lengths.

        Parameters
        ----------
        w040 : Quantity["length"]
            Spherical aberration (e.g. 27.5 nm)
        w131_edge : Quantity["length"]
            Coma at H = 1
        w222_edge : Quantity["length"]
            Astigmatism at H = 1
        w220_edge : Quantity["length"]
            Field curvature at H = 1
        reference_wavelength : Quantity["length"]
            Reference wavelength for ``w_rms_waves``

        Returns
        -------
        FieldAberrationModel
        """
        return cls(
            w040=w040,
            w131_edge=w131_edge,
            w222_edge=w222_edge,
            w220_edge=w220_edge,
            reference_wavelength=reference_wavelength,
        )

    # ---- core methods ----

    def w_rms(self, h: float | np.ndarray) -> Quantity:
        """Total RMS wavefront error at normalised field position *h*.

        Parameters
        ----------
        h : float or array
            Normalised radial field coordinate, 0 ≤ h ≤ 1

        Returns
        -------
        Quantity["length"]
            RMS WFE (in nm)
        """
        # Per-term variance (Mahajan balanced):
        #   sigma_i^2 = (c_i * W_ijk * h^p)^2
        # where p is the field exponent (0 for spherical, 1 for coma, 2 for
        # astig and field curvature).
        h = np.asarray(h, dtype=float)

        var_sph = (_RMS_COEFF_SPHERICAL * self.w040) ** 2
        var_coma = (_RMS_COEFF_COMA * self.w131_edge * h) ** 2
        var_astig = (_RMS_COEFF_ASTIGMATISM * self.w222_edge * h**2) ** 2
        var_fc = (_RMS_COEFF_FIELD_CURVATURE * self.w220_edge * h**2) ** 2

        total_var = var_sph + var_coma + var_astig + var_fc

        return np.sqrt(total_var).to(u.nm)

    def w_rms_waves(
        self, h: float | np.ndarray, wavelength: Quantity | None = None
    ) -> float | np.ndarray:
        """Total RMS wavefront error in waves.

        Parameters
        ----------
        h : float or array
            Normalised radial field coordinate
        wavelength : Quantity["length"], optional
            Wavelength to normalise by (defaults to ``reference_wavelength``)

        Returns
        -------
        float or array
            Dimensionless RMS WFE in waves
        """
        wvl = (wavelength if wavelength is not None else self.reference_wavelength).to(
            u.nm
        )
        return (self.w_rms(h) / wvl).decompose().value

    def to_zernikes(self, h_x: float, h_y: float, n_terms: int = 15) -> Quantity:
        """ANSI-ordered Zernike coefficient vector at field point (h_x, h_y).

        The mapping follows Mahajan (Ch. 7) and Wyant & Creath.  Tilt
        balance for coma is omitted (tilt only shifts the PSF centroid
        and does not affect MTF).

        Parameters
        ----------
        h_x, h_y : float
            Normalised field coordinates, each in [-1, 1]
        n_terms : int, optional
            Length of the returned vector (default 15, up to ANSI j = 14)

        Returns
        -------
        Quantity
            Zernike coefficient array (in nm), ANSI-ordered
        """
        h = np.hypot(h_x, h_y)
        phi_h = np.arctan2(h_y, h_x)  # field azimuth

        z = np.zeros(n_terms) * u.nm

        # --- Spherical: W_040 → Z11 (ANSI j=11, n=4,m=0) with Z4 balance ---
        # Peak-to-RMS balanced spherical: Z11 coeff = W_040 / (6*sqrt(5))
        # but we store the *peak* coefficient and the Zernike expansion
        # uses peak coefficients directly in the Noll/ANSI expansion.
        # Standard mapping (Mahajan Table 7.2):
        #   Z11 = W_040                (primary spherical)
        #   Z4  += -W_040              (defocus balance for spherical)
        if n_terms > 11:
            z[11] = self.w040  # type: ignore[index]
        if n_terms > 4:
            z[4] = z[4] - self.w040  # type: ignore[index]  # defocus balance

        # --- Coma: W_131 → Z7/Z8 (ANSI j=7,8; n=3,m=±1) ---
        # At field (h_x, h_y), peak coma = W_131 * h, rotated by phi_h.
        #   Z7 = W_131 * h * cos(phi_h)    (coma-x)
        #   Z8 = W_131 * h * sin(phi_h)    (coma-y)
        # Tilt balance (Z1/Z2) omitted — no MTF effect.
        w131_at_h = self.w131_edge * h
        if n_terms > 8:
            z[7] = w131_at_h * np.cos(phi_h)  # type: ignore[index]
            z[8] = w131_at_h * np.sin(phi_h)  # type: ignore[index]

        # --- Astigmatism: W_222 → Z5/Z6 (ANSI j=5,6; n=2,m=±2) ---
        # At field (h_x, h_y), peak astigmatism = W_222 * h^2, at azimuth 2*phi_h.
        #   Z5 = W_222 * h^2 * cos(2*phi_h)
        #   Z6 = W_222 * h^2 * sin(2*phi_h)
        #   Z4 += -W_222 * h^2 / 2         (defocus balance for astigmatism)
        w222_at_h = self.w222_edge * h**2
        if n_terms > 6:
            z[5] = w222_at_h * np.cos(2 * phi_h)  # type: ignore[index]
            z[6] = w222_at_h * np.sin(2 * phi_h)  # type: ignore[index]
        if n_terms > 4:
            z[4] = z[4] - w222_at_h / 2  # type: ignore[index]  # defocus balance

        # --- Field curvature: W_220 → Z4 (ANSI j=4, n=2,m=0) ---
        w220_at_h = self.w220_edge * h**2
        if n_terms > 4:
            z[4] = z[4] + w220_at_h  # type: ignore[index]

        return z
