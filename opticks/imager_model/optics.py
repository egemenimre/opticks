# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import numpy as np
from pint import Quantity
from strictyaml import Map, Str

from opticks import u
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.yaml_helpers import Qty

optics_schema = {
    "name": Str(),
    "focal_length": Qty(),
    "aperture_diameter": Qty(),
    "image_diam_on_focal_plane": Qty(),
}
"""Schema containing optical parameters."""


class Optics(ImagerComponent):
    """
    Class containing generic Optics parameters.
    """

    @classmethod
    def schema(cls) -> Map:
        return Map(optics_schema)

    @classmethod
    def _params_class_name(cls) -> str:
        return "OpticsParams"

    # ---------- begin modelling functions ----------

    @property
    def f_number(self) -> float:
        r"""
        F-number.

        Computed as: $\frac{\text{focal length}}{\text{aperture diameter}}$

        """
        return (
            (self.params.focal_length / self.params.aperture_diameter)
            .to_reduced_units()
            .m
        )

    @property
    def full_optical_fov(self) -> Quantity:
        r"""
        Full optical Field of View.

        Actual FoV depends on the detector size, but cannot be wider than this value.

        """
        return 2 * np.arctan(
            (self.params.image_diam_on_focal_plane / 2.0) / self.params.focal_length
        ).to(u.deg)

    @property
    def aperture_area(self) -> Quantity:
        r"""
        Aperture area.

        Computed as: $\pi \times \left( \frac{ \text{aperture diameter}}{2} \right)^2$
        """
        return np.pi * (self.params.aperture_diameter / 2.0) ** 2

    @property
    def aperture_solid_angle(self) -> Quantity:
        r"""
        Aperture solid angle in steradians.
        """
        return (
            np.pi
            / (self.params.focal_length / (self.params.aperture_diameter / 2.0)) ** 2
            * u.steradian
        )

    @u.check(None, "[length]")
    def spatial_cutoff_freq(self, ref_wavelength: Quantity) -> Quantity:
        r"""
        Spatial cutoff frequency, assumes perfect incoherent optics.

        Determines the theoretical limit of the optical resolution, or the smallest
        object resolvable by the optical system.

        Computed as: `1/(wavelength F_nr)` in line pairs per mm

        Parameters
        ----------
        ref_wavelength : Quantity
            Reference wavelength

        Returns
        -------
        Quantity
            Spatial cutoff frequency (in lp/mm)
        """
        u.define("line_pair = 1 * cycle = lp = lp")
        # perfect incoherent optics
        return (1.0 * u.lp) / (ref_wavelength * self.f_number).to(u.mm)

    @u.check(None, None, "[length]")
    def ideal_mtf(self, input_line_freq: Quantity, wavelength: Quantity) -> float:
        """
        Ideal optical MTF for the given input line frequency.

        Assumes uniformly illuminated circular aperture, no significant aberrations.

        Returns the MTF value between 0 and 1.

        Parameters
        ----------
        wavelength : Quantity
            Wavelength at which MTF is computed
        input_line_freq: Quantity
            Input line frequency (in lp/mm)

        Returns
        -------
        Quantity
            MTF value
        """

        # normalised optical frequency
        f_ov_fc = input_line_freq / self.spatial_cutoff_freq(wavelength)

        # This is the alternative formulation
        # psi = np.arccos(f_ov_fc)
        # mtf_ideal_optical = 2/np.pi * (psi-np.cos(psi)*np.sin(psi))

        return 2 / np.pi * (np.arccos(f_ov_fc) - f_ov_fc * np.sqrt(1 - f_ov_fc**2)).m
