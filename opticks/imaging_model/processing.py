# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Image processing component: resampling, sharpening, and related MTF models.
"""

from typing import Self

from astropy.units import Quantity
from pydantic import BaseModel, ConfigDict, model_validator

from opticks.contrast_model.mtf import MTF_Model_1D
from opticks.contrast_model.processing_mtf import (
    ResamplingKernel,
    validate_resampling_params,
)
from opticks.imaging_model.imager_component import ImagerComponent
from opticks.utils.parser_helpers import PositivePydanticQty


class ProcessingParams(BaseModel):
    """Image processing configuration parameters.

    Stores the fixed (pipeline-level) processing settings.  Runtime-varying
    inputs (local SSD) are passed as arguments to the ``get_*`` methods on
    ``Processing``.  The ``output_pitch`` is stored here as the default
    output grid pitch but can be overridden at call time.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    # Resampling kernel selection
    resampling_kernel: str | None = None
    # Default output grid pitch (can be overridden per call)
    output_pitch: PositivePydanticQty | None = None
    # Kernel-specific shape parameters
    bicubic_a: float | None = None
    lanczos_n: int | None = None

    @model_validator(mode="after")
    def _check_resampling_params(self) -> Self:
        if self.resampling_kernel is None:
            return self

        if self.output_pitch is None:
            raise ValueError("output_pitch is required when resampling_kernel is set.")

        kernel = ResamplingKernel[self.resampling_kernel.upper()]
        validate_resampling_params(kernel, self.bicubic_a, self.lanczos_n)

        return self


class Processing(ImagerComponent):
    """Image processing component (resampling, sharpening, etc.).

    Loaded from ``processing.yaml``. Holds the fixed processing configuration;
    the local SSD (``input_pitch``) is supplied at call time and can vary
    per image region.
    """

    name: str
    processing_params: ProcessingParams | None = None

    # ---------- begin modelling functions ----------

    def get_resampling_mtf_1d(
        self,
        input_pitch: Quantity,
        output_pitch: Quantity | None = None,
    ) -> MTF_Model_1D:
        """Return the resampling MTF model for a given local input pitch.

        Vary ``input_pitch`` (= local SSD on the ground) per image region to
        map the SSD-driven MTF variation across the frame.  ``output_pitch``
        defaults to the value in ``processing_params`` but can be overridden
        at call time (e.g. to compare different output resolutions).

        Parameters
        ----------
        input_pitch : Quantity["length"]
            Local input sample spacing (= SSD on the ground for ortho
            correction). Varies per image region.
        output_pitch : Quantity["length"], optional
            Output resampling grid pitch.  Defaults to
            ``processing_params.output_pitch``.

        Returns
        -------
        MTF_Model_1D
            Resampling MTF model for this (input_pitch, output_pitch) pair.

        Raises
        ------
        ValueError
            If ``processing_params`` is not configured, ``resampling_kernel``
            is not set, or no ``output_pitch`` is available.
        """
        p = self.processing_params
        if p is None:
            raise ValueError(
                "processing_params is not configured. "
                "Add a processing_params block to the YAML."
            )

        kernel = p.resampling_kernel
        if kernel is None:
            raise ValueError(
                "resampling_kernel is not configured in processing_params."
            )

        resolved_output_pitch = (
            output_pitch if output_pitch is not None else p.output_pitch
        )
        if resolved_output_pitch is None:
            raise ValueError(
                "output_pitch must be supplied either in processing_params "
                "or as an argument to get_resampling_mtf_1d()."
            )

        return MTF_Model_1D.resampling(
            kernel=kernel,
            input_pitch=input_pitch,
            output_pitch=resolved_output_pitch,
            bicubic_a=p.bicubic_a,
            lanczos_n=p.lanczos_n,
        )
