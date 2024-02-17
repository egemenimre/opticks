# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from strictyaml import Map, Str, Enum, Int, Optional

from opticks import u
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.yaml_helpers import Qty

detector_schema = {
    "name": Str(),
    "detector type": Enum(["pushbroom", "full frame"]),
    "pixel pitch": Qty(),
    "horizontal pixels": Int(),
    "vertical pixels": Int(),
    "horizontal pixels used": Int(),
    "vertical pixels used": Int(),
    Optional("binning", default=1): Int(),
    # TDI ignored for full frame detectors
    Optional("tdi stages", default=1): Int(),
    Optional("full well capacity"): Qty(),
    Optional("noise"): Map(
        {
            Optional("dark current"): Qty(),
            Optional("temporal dark noise"): Qty(),
        }
    ),
    "timings": Map(
        {
            # TODO revalid: required for frame-rate imager only
            Optional("frame rate", default=None): Qty(),
            "integration duration": Qty(),
            Optional("frame overhead duration", default=0 * u.ms): Qty(),
            Optional("frame overlap duration", default=0 * u.ms): Qty(),
        }
    ),
}
"""Schema containing detector parameters."""

# TODO schema: what to do with zero Optionals and None Optionals


class Detector(ImagerComponent):
    """
    Class containing generic Detector parameters.
    """

    @classmethod
    def schema(cls) -> Map:
        return Map(detector_schema)

    # ---------- begin modelling functions ----------
