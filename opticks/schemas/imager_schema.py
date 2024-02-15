# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from strictyaml import Map, Str, Enum, Int, Optional, Float, Bool

from opticks import u
from opticks.utils.yaml_helpers import Qty

# An Imager is made up of two necessary and one optional part
# Optics, Detector and Read-out / Write Electronics (Optional)
# There are separate schemas for each.


optics_schema = {
    "name": Str(),
    "focal length": Qty(),
    "aperture diameter": Qty(),
    "image diameter on focal plane": Qty(),
}
"""Schema containing optical parameters."""

# TODO what to do with zero Optionals and None Optionals

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
