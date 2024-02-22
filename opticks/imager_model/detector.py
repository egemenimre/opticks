# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
from pint import Quantity
from strictyaml import Map, Str, Enum, Int, Optional

from opticks import u
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.yaml_helpers import Qty


detector_schema = {
    "name": Str(),
    "detector_type": Enum(["pushbroom", "full frame"]),
    "pixel_pitch": Qty(),
    "horizontal_pixels": Int(),
    "vertical_pixels": Int(),
    "horizontal_pixels_used": Int(),
    "vertical_pixels_used": Int(),
    Optional("binning", default=1): Int(),
    # TDI ignored for full frame detectors
    Optional("tdi_stages", default=1): Int(),
    Optional("full_well_capacity"): Qty(),
    Optional("noise"): Map(
        {
            Optional("dark_current"): Qty(),
            Optional("temporal_dark_noise"): Qty(),
        }
    ),
    "timings": Map(
        {
            # TODO revalid: required for frame-rate imager only
            Optional("frame_rate", default=None): Qty(),
            "integration_duration": Qty(),
            Optional("frame_overhead_duration", default=0 * u.ms): Qty(),
            Optional("frame_overlap_duration", default=0 * u.ms): Qty(),
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

    def _params_class_name(cls) -> str:
        return "DetectorParams"

    # ---------- begin modelling functions ----------
