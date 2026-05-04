# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Detector-level MTF models (diffusion, CTE, crosstalk, etc.).

Currently includes:

- Diffusion MTF based on Crowell & Labuda (1969) with five simplified cases
  covering BSI and FSI detector geometries, plus preset parameter sets for
  common detector categories.
- Crosstalk MTF for the center-pixel 8-neighbour model.
- CTE (Charge Transfer Efficiency) MTF for CCD detectors.
"""

from enum import Enum
from typing import cast

import numpy as np
from astropy.units import Quantity
from numpy.typing import NDArray

from opticks import u


class DetectorDiffusionModel(Enum):
    """Diffusion MTF model variant (geometry + boundary conditions).

    Each member carries its ``required_params`` - the set of physical
    parameter names the model needs.  Any optional parameter *not* in
    this set is forbidden (triggers ``ValueError`` if supplied).
    """

    BSI_1 = ("bsi_1", {"diffusion_length", "field_free_depth", "depletion_depth"})
    BSI_2 = ("bsi_2", {"diffusion_length", "field_free_depth", "depletion_depth"})
    BSI_3 = (
        "bsi_3",
        {
            "diffusion_length",
            "field_free_depth",
            "surface_recomb_velocity",
            "diffusion_coeff",
        },
    )
    FSI_1 = ("fsi_1", {"diffusion_length", "field_free_depth", "depletion_depth"})
    FSI_2 = ("fsi_2", {"diffusion_length", "depletion_depth"})

    def __init__(self, label: str, required_params: set[str]):
        self.label = label
        self.required_params = required_params

    def __str__(self):
        return self.label


class DetectorDiffusionPreset(Enum):
    """Generic detector category presets for the diffusion MTF model.

    Each member carries its associated ``model`` and default ``params``
    dict (Quantity values keyed by parameter name).
    """

    SCIENTIFIC_CCD = (
        "scientific_ccd",
        DetectorDiffusionModel.BSI_1,
        {
            "diffusion_length": 200 * u.um,
            "field_free_depth": 17 * u.um,
            "depletion_depth": 15 * u.um,
        },
    )
    SCIENTIFIC_CMOS = (
        "scientific_cmos",
        DetectorDiffusionModel.BSI_2,
        {
            "diffusion_length": 20 * u.um,
            "field_free_depth": 5 * u.um,
            "depletion_depth": 1.5 * u.um,
        },
    )
    MCT_MWIR = (
        "mct_mwir",
        DetectorDiffusionModel.BSI_3,
        {
            "diffusion_length": 25 * u.um,
            "field_free_depth": 7 * u.um,
            "surface_recomb_velocity": 500 * u.cm / u.s,
            "diffusion_coeff": 15 * u.cm**2 / u.s,
        },
    )
    MCT_LWIR = (
        "mct_lwir",
        DetectorDiffusionModel.BSI_3,
        {
            "diffusion_length": 35 * u.um,
            "field_free_depth": 10 * u.um,
            "surface_recomb_velocity": 1000 * u.cm / u.s,
            "diffusion_coeff": 40 * u.cm**2 / u.s,
        },
    )
    CONSUMER_CMOS = (
        "consumer_cmos",
        DetectorDiffusionModel.FSI_1,
        {
            "diffusion_length": 15 * u.um,
            "field_free_depth": 3 * u.um,
            "depletion_depth": 2 * u.um,
        },
    )
    CONSUMER_CCD = (
        "consumer_ccd",
        DetectorDiffusionModel.FSI_2,
        {
            "diffusion_length": 50 * u.um,
            "depletion_depth": 3 * u.um,
        },
    )

    def __init__(self, label: str, model: DetectorDiffusionModel, params: dict):
        self.label = label
        self.model = model
        self.params = params

    def __str__(self):
        return self.label


# ---------- validation helpers ----------


def validate_diffusion_params(
    model: DetectorDiffusionModel,
    field_free_depth,
    depletion_depth,
    surface_recomb_velocity,
    diffusion_coeff,
) -> None:
    """Validate that supplied/missing parameters match the selected diffusion model."""
    supplied = {
        "field_free_depth": field_free_depth,
        "depletion_depth": depletion_depth,
        "surface_recomb_velocity": surface_recomb_velocity,
        "diffusion_coeff": diffusion_coeff,
    }
    required = model.required_params

    for param_name, value in supplied.items():
        if value is not None and param_name not in required:
            raise ValueError(
                f"Parameter '{param_name}' is not used by model {model}. "
                f"Remove it or choose a different model."
            )

    for param_name in required:
        if param_name != "diffusion_length" and supplied.get(param_name) is None:
            raise ValueError(
                f"Parameter '{param_name}' is required for model {model} "
                f"but was not supplied."
            )


def _check_mtf_range(mtf_array: NDArray[np.float64], model_id: str) -> None:
    """Raise ValueError if any MTF value is outside [0, 1]."""
    if np.any(mtf_array > 1) or np.any(mtf_array < 0):
        max_val = np.max(mtf_array)
        min_val = np.min(mtf_array)
        raise ValueError(
            f"MTF out of range [0, 1] for model '{model_id}': "
            f"min={min_val:.6g}, max={max_val:.6g}. "
            f"Check physical parameters."
        )


# ---------- collection efficiency (eta) functions ----------
# Pure numpy, all lengths in µm, frequency in cy/µm.


def _bsi_1_eta(
    f: NDArray[np.float64],
    alpha: float,
    L_o: float,
    L_a: float,
    L_b: float,
) -> NDArray[np.float64]:
    """BSI collection efficiency: dead back surface (S→∞), no reflection."""
    k = 2 * np.pi * f
    L = L_o / np.sqrt(1 + k**2 * L_o**2)
    xi = L_a / L
    aL = alpha * L
    return aL / (aL**2 - 1) * (
        (1 - np.cosh(xi) * np.exp(-alpha * L_a)) / np.sinh(xi)
        - np.exp(-alpha * L_a) / aL
    ) - np.exp(-alpha * L_b)


def _bsi_2_eta(
    f: NDArray[np.float64],
    alpha: float,
    L_o: float,
    L_a: float,
    L_b: float,
) -> NDArray[np.float64]:
    """BSI collection efficiency: passivated back surface (S=0), no reflection."""
    k = 2 * np.pi * f
    L = L_o / np.sqrt(1 + k**2 * L_o**2)
    xi = L_a / L
    aL = alpha * L
    return aL / (aL**2 - 1) * (
        aL / np.cosh(xi)
        - np.tanh(xi) * np.exp(-alpha * L_a)
        - np.exp(-alpha * L_a) / aL
    ) - np.exp(-alpha * L_b)


def _bsi_3_eta(
    f: NDArray[np.float64],
    alpha: float,
    L_o: float,
    L_a: float,
    s_over_d: float,
) -> NDArray[np.float64]:
    """BSI collection efficiency: thick substrate (L_b→∞), general surface S.

    Parameters
    ----------
    s_over_d : float
        S/D ratio in µm⁻¹ (computed as ``(S/D).decompose().to(1/u.um).value``)
    """
    k = 2 * np.pi * f
    L = L_o / np.sqrt(1 + k**2 * L_o**2)
    xi = L_a / L
    aL = alpha * L
    sLD = s_over_d * L
    cxi, sxi = np.cosh(xi), np.sinh(xi)
    denom = cxi + sLD * sxi
    return (
        aL
        / (aL**2 - 1)
        * (
            (aL + sLD) / denom
            - (sxi + sLD * cxi) / denom * np.exp(-alpha * L_a)
            - np.exp(-alpha * L_a) / aL
        )
    )


def _fsi_1_eta(
    f: NDArray[np.float64],
    alpha: float,
    L_o: float,
    L_a: float,
    L_D: float,
) -> NDArray[np.float64]:
    """FSI collection efficiency: dead far surface, finite field-free bulk."""
    k = 2 * np.pi * f
    L = L_o / np.sqrt(1 + k**2 * L_o**2)
    xi = L_a / L
    aL = alpha * L
    eta_dep = 1 - np.exp(-alpha * L_D)
    P = aL * L / (aL**2 - 1) * np.exp(-alpha * L_D)
    A = P * (np.cosh(xi) - np.exp(-alpha * L_a)) / np.sinh(xi)
    eta_ff = -(A / L - alpha * P)
    return eta_dep + eta_ff


def _fsi_2_eta(
    f: NDArray[np.float64],
    alpha: float,
    L_o: float,
    L_D: float,
) -> NDArray[np.float64]:
    """FSI collection efficiency: dead far surface, semi-infinite bulk (≡ Fiete-Sieb)."""
    k = 2 * np.pi * f
    L = L_o / np.sqrt(1 + k**2 * L_o**2)
    return 1 - np.exp(-alpha * L_D) / (1 + alpha * L)


# ---------- dispatcher ----------


def detector_diffusion_mtf(
    input_line_freq_cy_mm: float | NDArray[np.float64],
    model: DetectorDiffusionModel,
    alpha: float,
    L_o: float,
    L_a: float | None,
    depth: float | None,
    s_over_d: float | None,
) -> float | NDArray[np.float64]:
    """Dispatcher for the detector diffusion MTF.

    All length parameters are in µm, frequency in cy/mm (converted internally).
    Returns normalised MTF = η(f) / η(f≈0).
    """
    f = np.atleast_1d(np.asarray(input_line_freq_cy_mm, dtype=float)) / 1000  # cy/um
    f0 = np.array([1e-13])  # cy/um, approximation for f→0

    if model == DetectorDiffusionModel.BSI_1:
        eta = _bsi_1_eta(f, alpha, L_o, cast(float, L_a), cast(float, depth))
        eta0 = _bsi_1_eta(f0, alpha, L_o, cast(float, L_a), cast(float, depth))[0]
    elif model == DetectorDiffusionModel.BSI_2:
        eta = _bsi_2_eta(f, alpha, L_o, cast(float, L_a), cast(float, depth))
        eta0 = _bsi_2_eta(f0, alpha, L_o, cast(float, L_a), cast(float, depth))[0]
    elif model == DetectorDiffusionModel.BSI_3:
        eta = _bsi_3_eta(f, alpha, L_o, cast(float, L_a), cast(float, s_over_d))
        eta0 = _bsi_3_eta(f0, alpha, L_o, cast(float, L_a), cast(float, s_over_d))[0]
    elif model == DetectorDiffusionModel.FSI_1:
        eta = _fsi_1_eta(f, alpha, L_o, cast(float, L_a), cast(float, depth))
        eta0 = _fsi_1_eta(f0, alpha, L_o, cast(float, L_a), cast(float, depth))[0]
    elif model == DetectorDiffusionModel.FSI_2:
        eta = _fsi_2_eta(f, alpha, L_o, cast(float, depth))
        eta0 = _fsi_2_eta(f0, alpha, L_o, cast(float, depth))[0]
    else:
        raise ValueError(f"Unknown diffusion model: {model}")

    mtf = eta / eta0

    _check_mtf_range(mtf, str(model))

    scalar_input = np.ndim(input_line_freq_cy_mm) == 0
    return float(mtf[0]) if scalar_input else mtf


# ---------- crosstalk MTF (center pixel, 8 neighbours) ----------


def validate_crosstalk_params(xs: float, xd: float) -> None:
    """Validate crosstalk coefficients for the center-pixel model.

    Parameters
    ----------
    xs : float
        Side-neighbour crosstalk coefficient (dimensionless fraction).
    xd : float
        Diagonal-neighbour crosstalk coefficient (dimensionless fraction).
    """
    if xs < 0:
        raise ValueError(f"crosstalk_xs must be >= 0, got {xs}.")
    if xd < 0:
        raise ValueError(f"crosstalk_xd must be >= 0, got {xd}.")
    if 4 * xs + 4 * xd >= 1:
        raise ValueError(
            f"Kernel center weight (1 - 4*xs - 4*xd) must be positive. "
            f"Got xs={xs}, xd={xd} => 4*xs + 4*xd = {4 * xs + 4 * xd:.6g} >= 1."
        )


def detector_crosstalk_mtf(
    fx: Quantity,
    fy: Quantity,
    xs: float,
    xd: float,
    pixel_pitch: Quantity,
) -> float | NDArray[np.float64]:
    """Compute center-pixel (8-neighbour) crosstalk MTF at arbitrary 2D
    frequencies.

    The transfer function is real-valued for the symmetric center kernel::

        H(fx, fy) = (1 - 4Xs - 4Xd)
                   + 2Xs cos(2π fx P) + 2Xs cos(2π fy P)
                   + 4Xd cos(2π fx P) cos(2π fy P)

    Parameters
    ----------
    fx, fy : Quantity
        Spatial frequencies (e.g. in cy/mm).
    xs : float
        Side-neighbour crosstalk coefficient (dimensionless fraction).
    xd : float
        Diagonal-neighbour crosstalk coefficient (dimensionless fraction).
    pixel_pitch : Quantity["length"]
        Pixel pitch.

    Returns
    -------
    float | NDArray[np.float64]
        Crosstalk MTF value(s).
    """
    f_nyq = u.cy / (2 * pixel_pitch)
    cos_x = np.cos(np.pi * (fx / f_nyq).decompose().value * u.rad)
    cos_y = np.cos(np.pi * (fy / f_nyq).decompose().value * u.rad)
    mtf = (1 - 4 * xs - 4 * xd) + 2 * xs * (cos_x + cos_y) + 4 * xd * cos_x * cos_y

    _check_mtf_range(mtf, f"crosstalk (xs={xs}, xd={xd})")

    scalar_input = np.ndim(fx) == 0 and np.ndim(fy) == 0
    return float(mtf) if scalar_input else np.asarray(mtf, dtype=float)


def detector_crosstalk_mtf_1d(
    input_line_freq: Quantity,
    xs: float,
    xd: float,
    pixel_pitch: Quantity,
) -> float | NDArray[np.float64]:
    """1D slice of the center-pixel crosstalk MTF (fy = 0).

    Equivalent to ``1 - 2(Xs + 2Xd)(1 - cos(2π f P))``.
    When ``xd = 0`` this reduces to the classical nearest-neighbour formula
    ``1 - 2Xs(1 - cos(2π f P))``.
    """
    return detector_crosstalk_mtf(input_line_freq, 0 * u.cy / u.mm, xs, xd, pixel_pitch)


# ---------- CTE MTF (Charge Transfer Efficiency) ----------


def validate_cte_params(
    num_pixels: int, num_phases: int, cte: float, tdi_stages: int
) -> None:
    """Validate parameters for the CTE MTF model.

    Parameters
    ----------
    num_pixels : int
        Number of pixels between the target pixel and the readout register.
    num_phases : int
        Number of clock phases per pixel (typically 3 for 3-phase CCDs).
    cte : float
        Charge transfer efficiency per transfer (dimensionless, 0 to 1).
    tdi_stages : int
        Number of on-chip analog TDI stages (0 for non-TDI or digital TDI).
    """
    if num_pixels < 0:
        raise ValueError(f"num_pixels must be >= 0, got {num_pixels}.")
    if num_phases < 1:
        raise ValueError(f"num_phases must be >= 1, got {num_phases}.")
    if tdi_stages < 0:
        raise ValueError(f"tdi_stages must be >= 0, got {tdi_stages}.")
    if not (0.0 <= cte <= 1.0):
        raise ValueError(f"cte must be in [0, 1], got {cte}.")


def detector_cte_mtf_1d(
    input_line_freq: Quantity,
    num_pixels: int,
    num_phases: int,
    cte: float,
    pixel_pitch: Quantity,
    tdi_stages: int = 0,
) -> float | NDArray[np.float64]:
    """1D CTE MTF for a CCD detector (Janesick 2001).

    Models charge transfer inefficiency in the parallel or serial register.
    The MTF degrades with the number of charge transfers and the per-transfer
    inefficiency (1 - CTE):

        MTF(f) = exp(-N * (1 - cte) * (1 - cos(π f / f_nyq)))

    where N = (num_pixels + tdi_stages) * num_phases and
    f_nyq = 1 / (2 * pixel_pitch).

    Parameters
    ----------
    input_line_freq : Quantity
        Spatial frequency (e.g. in cy/mm).
    num_pixels : int
        Number of pixels between the target pixel and the readout register.
    num_phases : int
        Number of clock phases per pixel (typically 3 for 3-phase CCDs).
    cte : float
        Charge transfer efficiency per transfer (dimensionless, 0 to 1).
    pixel_pitch : Quantity["length"]
        Pixel pitch.
    tdi_stages : int, optional
        Number of on-chip analog TDI stages. Default is 0. Pass 0 for
        digital ("shift-and-add") TDI and for the ACT direction.

    Returns
    -------
    float | NDArray[np.float64]
        CTE MTF value(s).
    """
    num_transfers = (num_pixels + tdi_stages) * num_phases
    f_nyq = u.cy / (2 * pixel_pitch)
    cos_term = np.cos(np.pi * (input_line_freq / f_nyq).decompose())
    arg = num_transfers * (1.0 - cte) * (1.0 - cos_term)
    mtf = np.exp(-arg)
    _check_mtf_range(mtf, f"cte (N={num_transfers}, cte={cte:.6g})")

    scalar_input = np.ndim(input_line_freq) == 0
    return float(mtf) if scalar_input else np.asarray(mtf, dtype=float)
