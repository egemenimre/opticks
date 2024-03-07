# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from pathlib import Path

from strictyaml import YAML

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
