# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

from collections import namedtuple
from pathlib import Path
from typing import Iterable

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
        Initialises an Imager from components defined in YAML data.

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
        Initialises an Imager from components in YAML files.

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
        Initialises an Imager from components in YAML text.

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

    @u.check(None, "[length]", None, None)
    def q_factor(
        self,
        wavelength: Quantity | np.ndarray[Quantity],
        band_id: str,
        with_binning=True,
    ):
        """
        Computes the Q Factor.

        Q = wavelength x f_number / pixel pitch

        Parameters
        ----------
        wavelength : Quantity | ArrayLike[Quantity]
            Wavelength at which Q is computed
        band_id : str
            band ID (to select the correct band or channel)
        with_binning : bool, optional
            Return the value with binning or not

        Returns
        -------
        float
            Q factor value
        """

        # select the correct channel
        channel = self.detector.params.channels.all[band_id]

        pixel_pitch = channel.pixel_pitch(with_binning)

        q = (wavelength * self.optics.f_number / pixel_pitch).to_reduced_units()

        return q

    def horizontal_fov(self, band_id: str) -> Quantity:
        """
        Computes the full field of view in the horizontal direction.

        Used pixels only.

        Parameters
        ----------
        band_id : str
            band ID (to select the correct band or channel)

        Returns
        -------
        Quantity
            Horizontal FOV angle

        """

        # select the correct channel
        channel = self.detector.get_channel(band_id)

        return 2 * np.arctan(
            (
                (
                    channel.pixel_pitch(with_binning=False)
                    * channel.horizontal_pixels
                    / 2.0
                )
                / self.optics.params.focal_length
            )
        ).to(u.deg)

    def vertical_fov(self, band_id: str) -> Quantity:
        """
        Computes the full field of view in the vertical direction.

         Used pixels only.

        Parameters
        ----------
        band_id : str
            band ID (to select the correct band or channel)

        Returns
        -------
        Quantity
            Vertical FOV angle

        """
        # select the correct channel
        channel = self.detector.get_channel(band_id)

        return 2 * np.arctan(
            (
                (
                    channel.pixel_pitch(with_binning=False)
                    * channel.vertical_pixels
                    / 2.0
                )
                / self.optics.params.focal_length
            )
        ).to(u.deg)

    def ifov(self, band_id: str, with_binning=True) -> Quantity:
        """
        Computes the average instantaneous field of view
        (works in vertical and horizontal).

        The IFoV is computed simply as the FoV divided by the number of
        (binned or unbinned) pixels.

        Parameters
        ----------
        band_id : str
            band ID (to select the correct band or channel)
        with_binning : bool, optional
            Return the value with binning or not

        Returns
        -------
        Quantity
            IFOV angle

        """

        # select the correct channel
        channel = self.detector.get_channel(band_id)

        # horizontal pixels with binning included
        horiz_pixels = channel.pixel_count_line(with_binning).to(u.pixel).m

        return (self.horizontal_fov(band_id) / horiz_pixels).to(u.mdeg)

    def pix_solid_angle(self, band_id: str, with_binning=True) -> Quantity:
        """
        Pixel solid angle (of a pyramid).

        Parameters
        ----------
        band_id : str
            band ID (to select the correct band or channel)
        with_binning : bool, optional
            Return the value with binning or not

        Returns
        -------
        Quantity
            Pixel solid angle in steradians
        """

        pix_solid_angle = 4 * np.arcsin(
            np.sin(self.ifov(band_id, with_binning) / 2.0)
            * np.sin(self.ifov(band_id, with_binning) / 2.0)
        )

        # correct the unit from rad to sr (or rad**2)
        return (pix_solid_angle * u.rad).to(u.steradian)

    def data_write_rate(
        self,
        band_id: str | Iterable[str],
        with_binning: bool = True,
        with_compression: bool = True,
    ) -> Quantity:
        """
        Data write rate with or without compression.

        For a pushbroom only a single line is computed. For a full-frame,
        the entire frame is computed.

        If a list of channels is given, then  the returned result is the sum of all
        the requested channels.

        Note that the unused pixels are also read, this assumes that the
        detector does not have ROI functionality.

        Parameters
        ----------
        band_id : str or Iterable[str]
            band ID (to select the correct band or channel)
        with_binning : bool, optional
            Return the value with binning or not
        with_compression : bool
            Return the value with compression or not

        Returns
        -------
        Quantity
            Pixel read rate with or without binning (Mbit/s)
        """
        # TDI data is processed but not written unless raw data is needed
        with_tdi = False

        # data rate after encoding
        enc_data_rate = (
            self.detector.pix_read_rate(band_id, with_binning, with_tdi)
            * self.rw_electronics.params.pixel_encoding
        )

        # data rate after compression and other processing
        if with_compression:
            process_output_data_rate = (
                enc_data_rate / self.rw_electronics.params.compression_ratio
            )
        else:
            process_output_data_rate = enc_data_rate

        # data rate after overheads
        write_data_rate = process_output_data_rate * (
            1 + self.rw_electronics.params.data_write_overhead
        )

        return write_data_rate.to("Mbit/s")

    @u.check(None, "[length]", None)
    def projected_horiz_img_extent(
        self,
        distance: Quantity | np.ndarray[Quantity],
        band_id: str,
    ) -> Quantity:
        """
        Computes the projected horizontal image extent at some distance.

        The image extent is computed on a "flat surface" at the 'distance'.

        Parameters
        ----------
        distance : Quantity | np.ndarray[Quantity]
            _description_
        band_id : str
            band ID (to select the correct band or channel)

        Returns
        -------
        Quantity
            Projected horizontal image extent
        """

        return 2 * distance * np.tan(self.horizontal_fov(band_id) / 2.0)

    @u.check(None, "[length]", None)
    def projected_vert_img_extent(
        self,
        distance: Quantity | np.ndarray[Quantity],
        band_id: str,
    ) -> Quantity:
        """
        Computes the projected vertical image extent at some distance.

        The image extent is computed on a "flat surface" at the 'distance'.

        Parameters
        ----------
        distance : Quantity | np.ndarray[Quantity]
            _description_
        band_id : str
            band ID (to select the correct band or channel)

        Returns
        -------
        Quantity
            Projected horizontal image extent
        """

        return 2 * distance * np.tan(self.vertical_fov(band_id) / 2.0)

    @u.check(None, "[length]", None, None, None)
    def spatial_sample_distance(
        self,
        distance: Quantity | np.ndarray[Quantity],
        band_id: str,
        with_binning=True,
        location="centre",
    ) -> tuple[Quantity, Quantity]:
        """
        Computes the Spatial Sampling Distance (SSD) at some distance.

        The SSD is computed on a "flat surface" at the 'distance'.
        The location can be 'centre', 'centre left', 'centre right',
        'centre top', 'centre bottom' or 'corner'.

        - 'centre' corresponds to the boresight LoS vector corresponding
        to the channel pixels.
        - 'centre left', 'centre right' are equivalent and they correspond
        to the horizontal centre points or 3 o'clock and 9 o'clock
        positions of the channel pixels.
        - 'centre top', 'centre bottom' are equivalent and they correspond
        to the vertical centre points or 12 o'clock and 6 o'clock
        positions of the channel pixels
        - 'corner' corresponds to any corner of the channel pixels.

        As the ifov is constant, the horizontal and vertical SSD are not
        necessarily equal. This is particularly evident with large pixel
        sizes and large FoVs.

        The result is a 'namedtuple' and the horizontal and vertical
        SSD values can be queried with 'horiz' and 'vert', respectively.

        Parameters
        ----------
        distance : Quantity | np.ndarray[Quantity]
            _description_
        band_id : str
            band ID (to select the correct band or channel)
        with_binning : bool, optional
            Return the value with binning or not
        location : str, optional
            Location on the detector

        Returns
        -------
        (Quantity, Quantity)
            Spatial Sampling Distance with or without binning (horizontal and vertical)
        """

        # select the binned/unbinned ifov values
        ifov = self.ifov(band_id, with_binning)

        d = distance

        # compute the SSD, depending on the location
        if location == "centre":
            ssd_h = 2 * d * np.tan(ifov / 2.0)
            ssd_v = ssd_h

        elif location == "centre left" or location == "centre right":

            s_h = self.projected_horiz_img_extent(distance, band_id)
            fov_h = self.horizontal_fov(band_id)

            ssd_h = s_h / 2.0 - d * np.tan(fov_h / 2.0 - ifov)

            a_v = np.sqrt(d**2 + (s_h / 2) ** 2)
            ssd_v = 2 * a_v * np.tan(ifov / 2)

        elif location == "centre top" or location == "centre bottom":

            s_v = self.projected_vert_img_extent(distance, band_id)
            fov_v = self.vertical_fov(band_id)

            ssd_v = s_v / 2.0 - d * np.tan(fov_v / 2.0 - ifov)

            a_h = np.sqrt(d**2 + (s_v / 2) ** 2)
            ssd_h = 2 * a_h * np.tan(ifov / 2)

        elif location == "corner":

            s_v = self.projected_vert_img_extent(distance, band_id)
            fov_v = self.vertical_fov(band_id)

            s_h = self.projected_horiz_img_extent(distance, band_id)
            fov_h = self.horizontal_fov(band_id)

            a_h = np.sqrt(d**2 + (s_v / 2) ** 2)
            gamma_h = np.atan((s_h / 2) / a_h)
            ssd_h = s_h / 2.0 - a_h * np.tan(gamma_h - ifov)

            a_v = np.sqrt(d**2 + (s_h / 2) ** 2)
            gamma_v = np.atan((s_v / 2) / a_v)
            ssd_v = s_v / 2.0 - a_v * np.tan(gamma_v - ifov)

        else:
            raise ValueError(f"Invalid location flag: {location}")

        SSD = namedtuple("SSD", "horiz, vert")

        ssd = SSD(ssd_h.to_reduced_units().to(u.m), ssd_v.to_reduced_units().to(u.m))

        return ssd
