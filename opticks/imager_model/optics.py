# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import numpy as np
from pint import Quantity
from strictyaml import Str, Map

from opticks import u
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.yaml_helpers import Qty

optics_schema = {
    "name": Str(),
    "focal length": Qty(),
    "aperture diameter": Qty(),
    "image diameter on focal plane": Qty(),
}
"""Schema containing optical parameters."""


class Optics(ImagerComponent):
    """
    Class containing generic Optics parameters.
    """

    @classmethod
    def schema(cls) -> Map:
        return Map(optics_schema)

    # ---------- begin modelling functions ----------

    @property
    def f_number(self) -> float:
        return (
            (self.params["focal length"] / self.params["aperture diameter"])
            .to_reduced_units()
            .m
        )

    @property
    def full_optical_fov(self) -> Quantity:
        return 2 * np.arctan(
            (self.params["image diameter on focal plane"] / 2.0)
            / self.params["focal length"]
        ).to(u.deg)

    @property
    def aperture_area(self) -> Quantity:
        return np.pi * (self.params["aperture diameter"] / 2.0) ** 2

    @property
    def aperture_solid_angle(self) -> Quantity:
        # solid angle = 2pi h/r
        return (
            np.pi
            / (self.params["focal_length"] / (self.params["aperture diameter"] / 2.0))
            ** 2
            * u.steradian
        )

    def spatial_cutoff_freq(self, ref_wavelength) -> Quantity:
        u.define("lp = 1 * dimensionless = lp")
        # perfect incoherent optics
        return (1.0 * u.lp) / (ref_wavelength * self.f_number).to(u.mm)
