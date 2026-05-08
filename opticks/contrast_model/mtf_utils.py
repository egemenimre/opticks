# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Shared sanity-check helpers for the MTF model code.

Two patterns are consolidated here:

- ``reject_unused_params`` — guards against caller-supplied keyword
  parameters that the chosen model / kernel / preset does not actually
  use. Used by ``validate_*_params`` functions and preset overrides.
- ``check_mtf_below_unity`` / ``check_mtf_nonneg`` — atomic output
  sanity checks composed at call sites. Passive (detector / optics)
  MTF stages call both; processing-stage MTF (resampling kernels with
  edge-boost behaviour like Bicubic / Lanczos) calls only the
  non-negativity primitive.

Each primitive treats non-finite values (NaN, ±Inf) as failures, so
silent NaN propagation is not possible.
"""

from collections.abc import Iterable, Mapping

import numpy as np
from numpy.typing import NDArray


def reject_unused_params(
    supplied: Mapping[str, object],
    required: Iterable[str],
    *,
    model_name: str,
) -> None:
    """Raise ``ValueError`` if any supplied parameter is not in ``required``.

    Parameters supplied as ``None`` are ignored (treated as "not supplied").

    Parameters
    ----------
    supplied : Mapping[str, object]
        Mapping from parameter name to caller-supplied value.
    required : Iterable[str]
        Names of parameters the chosen model / kernel / preset accepts.
    model_name : str
        Human-readable identifier of the model / kernel / preset, used in
        the error message (e.g. ``f"model {model}"``).
    """
    required_set = required if isinstance(required, (set, frozenset)) else set(required)
    for name, value in supplied.items():
        if value is not None and name not in required_set:
            raise ValueError(
                f"Parameter '{name}' is not used by {model_name}. "
                f"Remove it or choose a different model."
            )


def check_mtf_below_unity(
    values: float | NDArray[np.float64],
    model_name: str,
) -> None:
    """Raise ``ValueError`` if any MTF value exceeds 1 (or is non-finite).

    Use for passive MTF stages where ``MTF > 1`` is unphysical (detector,
    optics). For processing-stage MTF that legitimately exceeds 1 (e.g.
    Bicubic / Lanczos edge-boost), call only :func:`check_mtf_nonneg`.

    Parameters
    ----------
    values : float | NDArray[np.float64]
        MTF value(s) to check (scalar or array).
    model_name : str
        Identifier of the MTF model, used in the error message.
    """
    arr = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"MTF contains non-finite values for '{model_name}'.")
    if np.any(arr > 1.0):
        raise ValueError(
            f"MTF exceeds 1 for '{model_name}': max={float(np.max(arr)):.6g}. "
            f"Check physical parameters."
        )


def check_mtf_nonneg(
    values: float | NDArray[np.float64],
    model_name: str,
) -> None:
    """Raise ``ValueError`` if any MTF value is negative (or is non-finite).

    Parameters
    ----------
    values : float | NDArray[np.float64]
        MTF value(s) to check (scalar or array).
    model_name : str
        Identifier of the MTF model, used in the error message.
    """
    arr = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"MTF contains non-finite values for '{model_name}'.")
    if np.any(arr < 0.0):
        raise ValueError(
            f"MTF has negative values for '{model_name}': "
            f"min={float(np.min(arr)):.6g}. Check physical parameters."
        )
