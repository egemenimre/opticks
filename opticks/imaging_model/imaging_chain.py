# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""Top-level imaging pipeline composed of an Imager and an optional Processing."""

from collections import namedtuple

import numpy as np
from astropy.units import Quantity

from opticks import u
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

    @u.quantity_input(distance="length")
    def projected_horiz_img_extent(
        self,
        distance: Quantity,
        band_id: str,
    ) -> Quantity:
        """
        Computes the projected horizontal image extent at some distance.

        The image extent is computed on a "flat surface" at the `distance`.

        Parameters
        ----------
        distance : Quantity
            distance to the target or plane
        band_id : str
            band ID (to select the correct band or channel)

        Returns
        -------
        Quantity
            Projected horizontal image extent
        """

        return 2 * distance * np.tan(self.imager.horizontal_fov(band_id) / 2.0)

    @u.quantity_input(distance="length")
    def projected_vert_img_extent(
        self,
        distance: Quantity,
        band_id: str,
    ) -> Quantity:
        """
        Computes the projected vertical image extent at some distance.

        The image extent is computed on a "flat surface" at the `distance`.

        Parameters
        ----------
        distance : Quantity
            distance to the target or plane
        band_id : str
            band ID (to select the correct band or channel)

        Returns
        -------
        Quantity
            Projected vertical image extent
        """

        return 2 * distance * np.tan(self.imager.vertical_fov(band_id) / 2.0)

    @u.quantity_input(distance="length")
    def spatial_sample_distance(
        self,
        distance: Quantity,
        band_id: str,
        with_binning=True,
        location="centre",
    ) -> tuple[Quantity, Quantity]:
        """
        Computes the Spatial Sampling Distance (SSD) at some distance.

        The SSD is computed on a "flat surface" at the `distance`.
        The location can be `centre`, `centre left`, `centre right`,
        `centre top`, `centre bottom` or `corner`.

        - `centre` corresponds to the boresight LoS vector corresponding
          to the channel pixels.
        - `centre left`, `centre right` are equivalent and they correspond
          to the horizontal centre points or 3 o'clock and 9 o'clock
          positions of the channel pixels.
        - `centre top`, `centre bottom` are equivalent and they correspond
          to the vertical centre points or 12 o'clock and 6 o'clock
          positions of the channel pixels
        - `corner` corresponds to any corner of the channel pixels.

        As the ifov is constant, the horizontal and vertical SSD are not
        necessarily equal. This is particularly evident with large pixel
        sizes and large FoVs.

        The result is a `namedtuple` and the horizontal and vertical
        SSD values can be queried with `horiz` and `vert`, respectively.

        Parameters
        ----------
        distance : Quantity
            distance between the imager and the target
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
        ifov = self.imager.ifov(band_id, with_binning)

        d = distance

        # compute the SSD, depending on the location
        if location == "centre":
            ssd_h = 2 * d * np.tan(ifov / 2.0)
            ssd_v = ssd_h

        elif location == "centre left" or location == "centre right":
            s_h = self.projected_horiz_img_extent(distance, band_id)
            fov_h = self.imager.horizontal_fov(band_id)

            ssd_h = s_h / 2.0 - d * np.tan(fov_h / 2.0 - ifov)

            a_v = np.sqrt(d**2 + (s_h / 2) ** 2)
            ssd_v = 2 * a_v * np.tan(ifov / 2)

        elif location == "centre top" or location == "centre bottom":
            s_v = self.projected_vert_img_extent(distance, band_id)
            fov_v = self.imager.vertical_fov(band_id)

            ssd_v = s_v / 2.0 - d * np.tan(fov_v / 2.0 - ifov)

            a_h = np.sqrt(d**2 + (s_v / 2) ** 2)
            ssd_h = 2 * a_h * np.tan(ifov / 2)

        elif location == "corner":
            s_v = self.projected_vert_img_extent(distance, band_id)
            fov_v = self.imager.vertical_fov(band_id)

            s_h = self.projected_horiz_img_extent(distance, band_id)
            fov_h = self.imager.horizontal_fov(band_id)

            a_h = np.sqrt(d**2 + (s_v / 2) ** 2)
            gamma_h = np.arctan((s_h / 2) / a_h)
            ssd_h = s_h / 2.0 - a_h * np.tan(gamma_h - ifov)

            a_v = np.sqrt(d**2 + (s_h / 2) ** 2)
            gamma_v = np.arctan((s_v / 2) / a_v)
            ssd_v = s_v / 2.0 - a_v * np.tan(gamma_v - ifov)

        else:
            raise ValueError(f"Invalid location flag: {location}")

        SSD = namedtuple("SSD", "horiz, vert")

        ssd = SSD(ssd_h.decompose().to(u.m), ssd_v.decompose().to(u.m))  # type: ignore[union-attr]

        return ssd
