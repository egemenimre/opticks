# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.

import numpy as np
from astropy.units import Quantity
from numpy import ndarray
from prysm.fttools import pad2d
from prysm.geometry import circle

from opticks import u
from opticks.utils.prysm_utils import Grid


class Aperture:
    grid: Grid
    """Grid object associated with the aperture"""

    def __init__(self, data: ndarray, grid: Grid) -> None:
        """Aperture class.

        Holds the Aperture data as defined by `prysm`. The aperture data is an
        ndarray mask of `bool` or 0 and 1 (or anything in between).

        Parameters
        ----------
        data : ndarray
            aperture data
        grid : Grid
            Grid object associated with the aperture
        """
        self.data = data
        self.grid = grid

    @classmethod
    def circle_aperture(
        cls,
        aperture_diam: Quantity,
        samples,
        with_units=True,
    ) -> "Aperture":
        """Circle aperture model.

        This is valid for many if not most refractive optics.

        This is a thin wrapper around the aperture building procedure
        from `prysm`.

        The output is an Aperture object containing `ndarray` of bools,
        therefore the unit input is not important for the aperture generation.
        But the backing Grid object does use units.

        Parameters
        ----------
        aperture_diam : Quantity
            aperture diameter in mm
        samples : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions.  Order is numpy axis convention, (row, col)
        with_units : bool, optional
            Grid returned with units (in mm and radians)

        Returns
        -------
        aperture: Aperture
            Aperture mask object
        """

        # aperture radius
        aperture_radius = aperture_diam / 2.0

        # generate the square grid in polar coords
        grid = Grid.from_size(samples, aperture_diam)
        r, t = grid.polar()

        # generate aperture (circle mask applied to the square grid)
        aperture = circle(aperture_radius, r)

        # strip units if needed
        if not with_units:
            grid = grid.strip_units(u.mm)

        return Aperture(aperture, grid)

    @classmethod
    def circle_aperture_with_obscuration(
        cls,
        aperture_diam: Quantity,
        obscuration_ratio: float,
        samples,
        with_units=True,
    ) -> "Aperture":
        """Circle aperture model with circular centre obscuration.

        This is valid for many reflective telescopes. The obscuration
        is usually the secondary mirror.

        This is a thin wrapper around the aperture building procedure from `prysm`.

        The output is an Aperture object containing `ndarray` of bools, therefore
        the unit input is not important for the aperture generation.

        Parameters
        ----------
        aperture_diam : Quantity
            aperture diameter in mm
        obscuration_ratio: float
            obscuration ratio (between 0 and 1)
        samples : int or tuple of int
            number of samples per dimension.  If a scalar value, broadcast to
            both dimensions.  Order is numpy axis convention, (row, col)
        with_units : bool, optional
            Grid returned with units (in mm and radians)

        Returns
        -------
        aperture: Aperture
            Aperture mask object
        """

        # aperture radius
        aperture_radius = aperture_diam / 2.0

        # generate the square grid in polar coords
        grid = Grid.from_size(samples, aperture_diam)
        r, t = grid.polar()

        # generate full aperture (circle mask applied to the square grid)
        full_aperture = circle(aperture_radius, r)

        # generate circular centre obscuration
        obscuration = circle(aperture_radius * obscuration_ratio, r)

        # apply full aperture minus obscuration as the mask
        aperture = full_aperture ^ obscuration  # or full_aperture & ~obscuration

        # can also use floats of 0 and 1 instead of bool
        # this also enables varying the amount of light through the aperture
        # aperture = aperture.astype(float)

        # strip units if needed
        if not with_units:
            grid = grid.strip_units(u.mm)

        return Aperture(aperture, grid)

    def scale_for_norm_sum_psf(self) -> "Aperture":
        """Generates a new, scaled Aperture that results in a
        Point Spread Function (PSF) that has a sum of 1.0.

        The aperture data (`data`) is divided by
        `sqrt(data.sum())`.

        """

        ap_data = self.data

        # scale for PSF sum = 1.0
        aperture_unity_sum = ap_data / np.sqrt(ap_data.sum())

        return Aperture(aperture_unity_sum, self.grid)

    def scale_for_norm_peak_psf(self, Q: float, Q_pad: float = 1) -> "Aperture":
        """Generates a new, scaled Aperture that results in a
        Point Spread Function (PSF) that has a peak of 1.0.

        `Q_pad` is used to pad the aperture if needed.
        It is passed on to the underlying `pad2d`.

        The aperture data (`data`) is multiplied by
        `Q x Q_pad x sqrt(data.size) / data.sum()`.

        Parameters
        ----------
        Q : float
            Q factor used in scaling pupil samples to psf samples
        Q_pad : float, optional
            padding factor, by default 1

        Returns
        -------
        Aperture
            New Aperture object with scaled data
        """

        ap_data = self.data

        # scale for PSF sum = 1.0
        aperture_padded = pad2d(ap_data, Q=Q_pad)  # type: ignore[arg-type]
        aperture_unity_peak = aperture_padded * (
            Q_pad * Q * np.sqrt(ap_data.size) / ap_data.sum()
        )

        return Aperture(aperture_unity_peak, self.grid)
