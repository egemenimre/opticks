# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pydantic import field_validator, model_validator

from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.parser_helpers import PydanticQty


class RWElectronics(ImagerComponent):
    """
    Class containing generic imager parameters.
    """

    name: str
    pixel_encoding: PydanticQty
    data_write_overhead: float
    compression_on: bool | None = None
    compression_ratio: float | None = None

    @field_validator("data_write_overhead")
    @classmethod
    def check_data_write_overhead(cls, v: float) -> float:
        if not 0.0 <= v <= 1.0:
            raise ValueError(
                f"data_write_overhead must be a fraction between 0.0 and 1.0, got {v}"
            )
        return v

    @model_validator(mode="after")
    def check_compression_ratio(self) -> "RWElectronics":
        if self.compression_on and self.compression_ratio is None:
            raise ValueError(
                "compression_ratio must be set when compression_on is True"
            )
        return self

    # ---------- begin modelling functions ----------
