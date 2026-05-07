# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""Top-level imaging pipeline composed of an Imager and an optional Processing."""

from opticks.imaging_model.imager import Imager
from opticks.imaging_model.processing import Processing


class ImagingChain:
    """Top-level imaging pipeline: an Imager and an optional Processing.

    Holds an Imager (the physical imaging hardware: optics, detector,
    read/write electronics) and an optional Processing component
    (resampling, sharpening, etc.) that operates on data after the
    imager produces it.

    Parameters
    ----------
    imager : Imager
        Imaging hardware (optics + detector + optional read/write electronics).
    processing : Processing, optional
        Post-acquisition processing component.
    """

    def __init__(
        self,
        imager: Imager,
        processing: Processing | None = None,
    ):
        self.imager = imager
        self.processing = processing
