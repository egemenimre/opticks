# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
from strictyaml import Bool, Float, Map, Optional, Str

from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.yaml_helpers import Qty

rw_electronics_schema = {
    "name": Str(),
    "pixel_encoding": Qty(),
    # TODO revalid: data write overhead is percentage
    "data_write_overhead": Float(),
    # TODO revalid: if compression ON, then we need the compression ratio
    Optional("compression_on"): Bool(),
    Optional("compression_ratio"): Float(),
}
"""Schema containing read-out / write electronics parameters."""


class RWElectronics(ImagerComponent):
    """
    Class containing generic imager parameters.
    """

    @classmethod
    def schema(cls) -> Map:
        return Map(rw_electronics_schema)

    @classmethod
    def _params_class_name(cls) -> str:
        return "RWElectronicsParams"

    # ---------- begin modelling functions ----------
