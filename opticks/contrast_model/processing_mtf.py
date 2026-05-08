# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Image-processing-stage MTF models (resampling, sharpening, etc.).

Currently includes:

- Resampling MTF for orthorectification / geometric correction with
  configurable interpolation kernels (nearest-neighbor, bilinear,
  bicubic, Lanczos, ideal sinc). Models both up- and downsampling
  via the ``p_eff = max(input_pitch, output_pitch)`` rule
  (Schowengerdt; Holst, *CMOS/CCD Sensors and Camera Systems*,
  2nd ed., Ch. 10).
"""

from enum import Enum

import numpy as np
from astropy.units import Quantity
from numpy.typing import NDArray
from scipy.integrate import simpson

from opticks import u
from opticks.contrast_model.mtf_utils import check_mtf_nonneg, reject_unused_params


class ResamplingKernel(Enum):
    """Interpolation kernel selection for the resampling MTF.

    Each member carries its ``required_params`` - the set of optional
    parameter names the kernel accepts.  Any optional parameter *not*
    in this set is forbidden (triggers ``ValueError`` if supplied).
    """

    NEAREST_NEIGHBOR = ("nearest_neighbor", set())
    BILINEAR = ("bilinear", set())
    BICUBIC = ("bicubic", {"bicubic_a"})
    LANCZOS = ("lanczos", {"lanczos_n"})
    SINC = ("sinc", set())

    def __init__(self, label: str, required_params: set[str]):
        self.label = label
        self.required_params = required_params

    def __str__(self):
        return self.label


# ---------- validation helpers ----------


def validate_resampling_params(
    kernel: ResamplingKernel,
    bicubic_a: float | None,
    lanczos_n: int | None,
) -> None:
    """Validate that supplied/missing parameters match the selected kernel.

    Parameters
    ----------
    kernel : ResamplingKernel
        Selected resampling kernel.
    bicubic_a : float or None
        Keys cubic shape parameter. Allowed only for ``BICUBIC``;
        must be in [-1.0, 0.0].
    lanczos_n : int or None
        Lanczos lobe count. Allowed only for ``LANCZOS``;
        must be an integer >= 2.
    """
    supplied = {"bicubic_a": bicubic_a, "lanczos_n": lanczos_n}
    required = kernel.required_params

    reject_unused_params(supplied, required, model_name=f"kernel {kernel}")

    if bicubic_a is not None and not (-1.0 <= bicubic_a <= 0.0):
        raise ValueError(f"bicubic_a must be in [-1.0, 0.0], got {bicubic_a}.")

    if lanczos_n is not None:
        if not isinstance(lanczos_n, int) or isinstance(lanczos_n, bool):
            raise ValueError(
                f"lanczos_n must be an integer, got {type(lanczos_n).__name__}."
            )
        if lanczos_n < 2:
            raise ValueError(f"lanczos_n must be >= 2, got {lanczos_n}.")


# ---------- per-kernel MTF helpers ----------
# All take a normalised frequency nu = f * p_eff (dimensionless,
# in cycles per p_eff length) as a 1-D array and return |MTF(nu)|.
# Numerical-FT kernels (Bicubic, Lanczos) divide by H(0) so that the
# returned MTF is normalised to 1 at DC even when the truncated kernel
# does not have unit integral.


def _nearest_neighbor_mtf(nu: NDArray[np.float64]) -> NDArray[np.float64]:
    """Box-filter / pixel-replication MTF: ``|sinc(nu)|``."""
    return np.abs(np.sinc(nu))


def _bilinear_mtf(nu: NDArray[np.float64]) -> NDArray[np.float64]:
    """Triangle-filter / linear-interpolation MTF: ``sinc^2(nu)``."""
    return np.sinc(nu) ** 2


def _bicubic_mtf(nu: NDArray[np.float64], a: float) -> NDArray[np.float64]:
    """Keys piecewise-cubic MTF (numerical FT of the Keys kernel).

    Kernel (even, support [-2, 2])::

        h(x) = (a+2)|x|^3 - (a+3)|x|^2 + 1            for |x| <= 1
        h(x) = a|x|^3 - 5a|x|^2 + 8a|x| - 4a          for 1 < |x| <= 2
    """
    x = np.linspace(0.0, 2.0, 4001)
    h = np.where(
        x <= 1.0,
        (a + 2) * x**3 - (a + 3) * x**2 + 1,
        a * x**3 - 5 * a * x**2 + 8 * a * x - 4 * a,
    )
    cos_term = np.cos(2 * np.pi * nu[:, None] * x[None, :])
    h_at_nu = 2 * simpson(h[None, :] * cos_term, x=x, axis=1)
    h_at_zero = 2 * simpson(h, x=x)
    return np.abs(h_at_nu / h_at_zero)


def _lanczos_mtf(nu: NDArray[np.float64], n: int) -> NDArray[np.float64]:
    """Lanczos-N windowed-sinc MTF (numerical FT of the Lanczos kernel).

    Kernel (even, support [-n, n])::

        h(x) = sinc(x) * sinc(x/n)   for |x| < n
    """
    x = np.linspace(0.0, float(n), 2001 * n)
    h = np.sinc(x) * np.sinc(x / n)
    cos_term = np.cos(2 * np.pi * nu[:, None] * x[None, :])
    h_at_nu = 2 * simpson(h[None, :] * cos_term, x=x, axis=1)
    h_at_zero = 2 * simpson(h, x=x)
    return np.abs(h_at_nu / h_at_zero)


def _sinc_mtf(nu: NDArray[np.float64]) -> NDArray[np.float64]:
    """Ideal brick-wall (sinc-kernel) MTF: ``rect(nu)``.

    Returns 1 for ``|nu| < 0.5``, else 0.  This is the limit case;
    the underlying kernel ``sinc(x)`` is non-realisable (infinite
    support) and produces strong ringing in practice.
    """
    return np.where(np.abs(nu) < 0.5, 1.0, 0.0)


# ---------- dispatcher ----------


def resampling_mtf_1d(
    input_line_freq: Quantity,
    kernel: ResamplingKernel,
    input_pitch: Quantity,
    output_pitch: Quantity,
    bicubic_a: float | None = None,
    lanczos_n: int | None = None,
) -> float | NDArray[np.float64]:
    """1-D resampling MTF for a chosen interpolation kernel.

    Models the MTF of an interpolation step that resamples from a grid
    of spacing ``input_pitch`` (= local SSD on the ground in the
    geometric-correction case) onto a grid of spacing ``output_pitch``
    (= ortho grid).  The kernel is scaled to::

        p_eff = max(input_pitch, output_pitch)

    so the same expression covers both upsampling and downsampling
    (the kernel doubles as anti-alias prefilter on downsample)::

        MTF(f) = |H_kernel(f * p_eff)|

    The returned value is always non-negative; the absolute value is
    taken because Bicubic, Lanczos, and the bare sinc all have
    negative side-lobes that would be misleading inside an MTF chain
    product.

    Parameters
    ----------
    input_line_freq : Quantity
        Spatial frequency (e.g. in cy/m or cy/mm).
    kernel : ResamplingKernel
        Selected resampling kernel.
    input_pitch : Quantity["length"]
        Local input sample spacing (= SSD on ground for ortho).
    output_pitch : Quantity["length"]
        Output resampling grid pitch.
    bicubic_a : float, optional
        Keys cubic shape parameter; defaults to -0.5.  Only used
        when ``kernel == BICUBIC``.
    lanczos_n : int, optional
        Lanczos lobe count; defaults to 3.  Only used when
        ``kernel == LANCZOS``.

    Returns
    -------
    float | NDArray[np.float64]
        Resampling MTF value(s) in [0, 1].
    """
    p_eff = input_pitch if input_pitch >= output_pitch else output_pitch

    nu_qty = (input_line_freq * p_eff).decompose()  # type: ignore[union-attr]
    nu = np.atleast_1d(nu_qty.to_value(u.cy)).astype(float)

    if kernel == ResamplingKernel.NEAREST_NEIGHBOR:
        mtf = _nearest_neighbor_mtf(nu)
    elif kernel == ResamplingKernel.BILINEAR:
        mtf = _bilinear_mtf(nu)
    elif kernel == ResamplingKernel.BICUBIC:
        mtf = _bicubic_mtf(nu, bicubic_a if bicubic_a is not None else -0.5)
    elif kernel == ResamplingKernel.LANCZOS:
        mtf = _lanczos_mtf(nu, lanczos_n if lanczos_n is not None else 3)
    elif kernel == ResamplingKernel.SINC:
        mtf = _sinc_mtf(nu)
    else:
        raise ValueError(f"Unknown resampling kernel: {kernel}")

    check_mtf_nonneg(mtf, str(kernel))

    scalar_input = np.ndim(input_line_freq) == 0
    return float(mtf[0]) if scalar_input else np.asarray(mtf, dtype=float)
