# opticks: Sizing Tool for Optical Systems
#
# Copyright (C) 2024 Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
from pint import Quantity
from strictyaml import Bool, Float, Map, Optional, Str

from opticks.imager_model.detector import Detector
from opticks.imager_model.imager_component import ImagerComponent
from opticks.utils.yaml_helpers import Qty

rw_electronics_schema = {
    "name": Str(),
    "pixel_encoding": Qty(),
    # TODO revalid: data write overhead is percentage
    "data_write_overhead": Float(),
    # TODO revalid: if compression ON, then we need the compression ratio
    Optional("compression_on"): Bool(),
    Optional("compression_ratio"): Float(),
}
"""Schema containing read-out / write electronics parameters."""


class RWElectronics(ImagerComponent):
    """
    Class containing generic imager parameters.
    """

    @classmethod
    def schema(cls) -> Map:
        return Map(rw_electronics_schema)

    @classmethod
    def _params_class_name(cls) -> str:
        return "RWElectronicsParams"

    # ---------- begin modelling functions ----------

    def data_write_rate(
        self,
        detector: Detector,
        with_binning: bool = True,
        with_compression: bool = True,
    ) -> Quantity:
        """
        Data write rate with or without compression.

        For a pushbroom only a single line is computed. For a full-frame,
        the entire frame is computed.

        Note that the unused pixels are also read, this assumes that the
        detector does not have ROI functionality.

        Parameters
        ----------
        detector : Detector
            Detector properties
        with_binning : bool
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
            detector.pix_read_rate(Channel, with_binning, with_tdi)
            * self.params.pixel_encoding
        )

        # data rate after compression and other processing
        if with_compression:
            process_output_data_rate = enc_data_rate / self.params.compression_ratio
        else:
            process_output_data_rate = enc_data_rate

        # data rate after overheads
        write_data_rate = process_output_data_rate * (
            1 + self.params.data_write_overhead
        )

        return write_data_rate.to("Mbit/s")
