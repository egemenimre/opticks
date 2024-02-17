# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from strictyaml import Str, Optional, Float, Bool, Map

from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.yaml_helpers import Qty

rw_electronics_schema = {
    "name": Str(),
    "pixel encoding": Qty(),
    # TODO revalid: data write overhead is percentage
    "data write overhead": Float(),
    # TODO revalid: if compression ON, then we need the compression ratio
    Optional("compression on"): Bool(),
    Optional("compression ratio"): Float(),
}
"""Schema containing read-out / write electronics parameters."""


class RWElectronics(ImagerComponent):
    """
    Class containing generic imager parameters.
    """

    @classmethod
    def schema(cls) -> Map:
        return Map(rw_electronics_schema)

    # ---------- begin modelling functions ----------
