# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

import numpy as np
from pint import Quantity
from strictyaml import YAML

from opticks import u
from opticks.imager_model.detector import Detector
from opticks.imager_model.optics import Optics
from opticks.imager_model.rw_electronics import RWElectronics


class Imager:
    """
    Class containing generic imager parameters.

    An Imager is made up of two necessary and one optional part:
    Optics, Detector and Read-out / Write Electronics (Optional)

    Parameters
    ----------
    optics_data : Optics
        Optics data object
    detector_data : Detector
        Detector data object
    rw_electronics_data : RWElectronics
        Read-out/Write Electronics data object
    """

    def __init__(
        self,
        optics_data: Optics,
        detector_data: Detector,
        rw_electronics_data: RWElectronics = None,
    ):

        self.optics = optics_data
        self.detector = detector_data
        if rw_electronics_data:
            self.rw_electronics = rw_electronics_data
        else:
            self.rw_electronics = None

    @classmethod
    def from_yaml_data(
        cls, optics_data: YAML, detector_data: YAML, rw_electronics_data: YAML = None
    ) -> "Imager":
        """
        Initialise an Imager from components defined in YAML data.

        Parameters
        ----------
        optics_data : YAML
            YAML containing optics design data
        detector_data : YAML
            YAML containing detector design data
        rw_electronics_data : YAML
            YAML containing read-write electronics design data

        Returns
        -------
        Imager
            Imager object from the input data
        """

        # extract data structures from YAML data
        optics = Optics(optics_data)
        detector = Detector(detector_data)
        rw_electronics = None
        if rw_electronics_data:
            rw_electronics = RWElectronics(rw_electronics_data)

        # return a new instance of Imager
        return cls(optics, detector, rw_electronics)

    @classmethod
    def from_yaml_file(
        cls, optics_data: Path, detector_data: Path, rw_electronics_data: Path = None
    ) -> "Imager":
        """
        Initialise an Imager from components in YAML files.

        Parameters
        ----------
        optics_data : Path
            Filepath containing optics design data file (YAML)
        detector_data : Path
            Filepath containing detector design data file (YAML)
        rw_electronics_data : Path
            Filepath containing read-write electronics design data file (YAML)

        Returns
        -------
        Imager
            Imager object from the input data
        """
        optics = Optics.from_yaml_file(optics_data)
        detector = Detector.from_yaml_file(detector_data)
        rw_electronics = None
        if rw_electronics_data:
            rw_electronics = RWElectronics.from_yaml_file(rw_electronics_data)

        # return a new instance of Imager
        return cls(optics, detector, rw_electronics)

    @classmethod
    def from_yaml_text(
        cls, optics_data: str, detector_data: str, rw_electronics_data: str = None
    ) -> "Imager":
        """
        Initialise an Imager from components in YAML text.

        Parameters
        ----------
        optics_data : str
            Text containing optics design data (YAML)
        detector_data : str
            Text containing detector design data (YAML)
        rw_electronics_data : str
            Text containing read-write electronics design data (YAML)

        Returns
        -------
        Imager
            Imager object from the input data
        """
        optics = Optics.from_yaml_text(optics_data)
        detector = Detector.from_yaml_text(detector_data)
        rw_electronics = None
        if rw_electronics_data:
            rw_electronics = RWElectronics.from_yaml_text(rw_electronics_data)

        # return a new instance of Imager
        return cls(optics, detector, rw_electronics)

    # ---------- begin modelling functions ----------

    def ifov(self, with_binning: bool = True) -> Quantity:
        """
        Computes the Instantaneous field of view (works in vertical and horizontal).

        Assumes constant IFOV per pixel.

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not

        Returns
        -------
        Quantity
            IFOV angle

        """
        return 2 * np.arctan(
            (self.detector.pix_pitch(with_binning) / 2.0)
            / self.optics.params.focal_length
        ).to(u.mdeg)

    def pix_solid_angle(self, with_binning=True) -> Quantity:
        """
        Pixel solid angle (of a pyramid).

        Parameters
        ----------
        with_binning : bool
            Return the value with binning or not

        Returns
        -------
        Quantity
            Pixel solid angle in steradians
        """
        #

        pix_solid_angle = 4 * np.arcsin(
            np.sin(self.ifov(with_binning) / 2.0)
            * np.sin(self.ifov(with_binning) / 2.0)
        )

        # correct the unit from rad to sr (or rad**2)
        return (pix_solid_angle * u.rad).to(u.steradian)

    @property
    def horizontal_fov(self) -> Quantity:
        """
        Computes the full field of view in the horizontal direction.

        Assumes constant IFOV per pixel.

        Parameters
        ----------

        Returns
        -------
        Quantity
            Horizontal FOV angle

        """
        return 2 * np.tan(
            self.ifov(False) * self.detector.params.horizontal_pixels_used / 2.0
        ).to(u.deg)

    @property
    def vertical_fov(self) -> Quantity:
        """
        Computes the full field of view in the vertical direction.

        Assumes constant IFOV per pixel.

        Parameters
        ----------

        Returns
        -------
        Quantity
            Vertical FOV angle

        """
        return 2 * np.tan(
            self.ifov(False) * self.detector.params.vertical_pixels_used / 2.0
        ).to(u.deg)
