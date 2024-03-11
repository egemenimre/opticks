import marimo

__generated_with = "0.2.13"
app = marimo.App()


@app.cell
def __():
    import marimo as mo

    mo.md("# Basic Monochrome MWIR Imager on a Drone")
    return mo,


@app.cell
def __(mo):
    mo.md(
        f"""
        ## Setting the Scene

        Computing the basic parameters of a monochrome imager on a drone. Valid for a single band.
        """
    )
    return


@app.cell
def __():
    import numpy as np

    from opticks import u

    # positional params
    # ---------------------
    distance = 10.0 * u.km
    return distance, np, u


@app.cell
def __(distance, mo):
    mo.md(
        f"""
        Reference parameters for the camera position and motion with respect to the target are given below:

        ```
        target distance : {distance:~P}   
        ```
        """
    )
    return


@app.cell
def __(mo):
    # init the ref wavelength slider in a separate cell
    slider = mo.ui.slider(3000, 5000, label="reference wavelength")

    mo.md(
        f"""
        Certain parameters vary with the light characteristics and therefore wavelength dependent.
        """
    )
    return slider,


@app.cell
def __(mo, slider, u):
    ref_wavelength = slider.value * u.nm

    mo.md(
        f"""
        Set the reference wavelength with the slider (limited by the detector sensitivity range):

        {slider} 
        ```
        reference wavelength : {ref_wavelength:~P}    
        ```

        Note the widely used bands:

        - blue: 450-485 nm
        - green: 500-565 nm
        - red: 625-750 nm
        - nir: 750-1400 nm
        - swir: 1400-3000 nm
        - mwir: 3000-8000 nm
        """
    )
    return ref_wavelength,


@app.cell
def __(mo):
    mo.md(
        f"""
        ## Loading the Imager Parameters

        An imager is made up of three parts: Optics, Detector and (optionally) Read-out/Write Electronics.
    """
    )
    return


@app.cell
def __(mo):
    optics_file = mo.ui.file(
        filetypes=[".yml", ".yaml", "txt"],
        multiple=False,
        kind="button",
    )

    detector_file = mo.ui.file(
        filetypes=[".yml", ".yaml"],
        multiple=False,
        kind="button",
    )

    rw_electronics_file = mo.ui.file(
        filetypes=[".yml", ".yaml"],
        multiple=False,
        kind="button",
    )
    return detector_file, optics_file, rw_electronics_file


@app.cell
def __(detector_file, mo, optics_file, rw_electronics_file):
    mo.md(
        f""" 
        The design parameters of each part is kept in a separate YAML file. Therefore, the imager needs to be initialised with at least the optics and detector and, if needed, Read-out/Write YAML data.

        First, the yaml files containing the optics and detector data are defined and then the `Imager` object is constructed from the yaml files (in this example the Read-out/Write Electronics are omitted).

        Please upload your optics and detector YAML files (loaded as UTF-8). Otherwise the default files will be used:

        {optics_file} optics file : {optics_file.name(0)}


        {detector_file} detector file : {detector_file.name(0)}

        {rw_electronics_file} rw electronics file : {rw_electronics_file.name(0)}

    """
    )
    return


@app.cell
def __(detector_file, optics_file, rw_electronics_file):
    from pathlib import Path

    from opticks.imager_model.imager import Imager

    if not optics_file.name(0) and not detector_file.name(0):

        # Going with the defaults
        file_directory = Path("docs", "examples", "sample_drone_imager")
        optics_file_path = file_directory.joinpath("optics.yaml")
        detector_file_path = file_directory.joinpath("detector.yaml")
        rw_electronics_file_path = file_directory.joinpath("rw_electronics.yaml")

        imager = Imager.from_yaml_file(
            optics_file_path, detector_file_path, rw_electronics_file_path
        )

    else:

        # extract yaml data from the user-supplied files (loaded by marimo as binary)
        optics_yaml = optics_file.contents(0).decode("utf-8")
        detector_yaml = detector_file.contents(0).decode("utf-8")

        # RW Electronics is optional
        rw_electronics_yaml = None
        if rw_electronics_file.name(0):
            rw_electronics_yaml = rw_electronics_file.contents(0).decode("utf-8")

        # Init imager object
        imager = Imager.from_yaml_text(
            optics_yaml, detector_yaml, rw_electronics_yaml
        )

    # shorthands
    optics = imager.optics
    detector = imager.detector
    rw_electronics = imager.rw_electronics

    # select the first channel
    channel = next(iter(detector.params.channels.all.values()))

    # binning status
    binning_on = False if channel.binning == 1 else True
    return (
        Imager,
        Path,
        binning_on,
        channel,
        detector,
        detector_file_path,
        detector_yaml,
        file_directory,
        imager,
        optics,
        optics_file_path,
        optics_yaml,
        rw_electronics,
        rw_electronics_file_path,
        rw_electronics_yaml,
    )


@app.cell
def __(mo, optics):
    mo.md(
        f"""
        ## Extracting the Imager Parameters

        ### Optical Parameters

        Basic optical parameters are given below for the optics with the name "_{optics.params.name}_":

        ```
        focal length : {optics.params.focal_length:~}
        aperture diameter : {optics.params.aperture_diameter:~}
        image diameter on focal plane : {optics.params.image_diam_on_focal_plane:~}
        ```
    """
    )
    return


@app.cell
def __(mo, optics, ref_wavelength):
    mo.md(
        f"""
        The derived optical parameters are given below for the optics with the name "_{optics.params.name}_":

        ```
        f-number : {optics.f_number:.4}
        full optical fov : {optics.full_optical_fov:~P.4}
        aperture area : {optics.aperture_area.to("cm**2"):~P.6}
        spatial cut-off freq  : {optics.spatial_cutoff_freq(ref_wavelength):~P.5} (at {ref_wavelength:~P})
        ```
    """
    )
    return


@app.cell
def __(binning_on, channel, detector, mo):
    mo.md(
        f"""
        ### Detector Parameters

        Basic design parameters are given below for the detector with the name "_{detector.params.name}_" and channel "_{channel.name}_":

        ```
        detector type : {detector.params.detector_type}
        horizontal x vertical pixels (detector) : {detector.params.horizontal_pixels} x {detector.params.vertical_pixels}
        horizontal x vertical pixels (channel) : {channel.horizontal_pixels} x {channel.vertical_pixels}
        binning : {channel.binning if binning_on else 'None'}
        pixel pitch : {channel.pixel_pitch(False):~P} {f"({channel.pixel_pitch(True):~P} binned)" if binning_on else ''}
        TDI stages : {'None' if channel.tdi_stages == 1 else channel.tdi_stages}

        Timing parameters:
        integration duration : {detector.params.timings.integration_duration:~}
        frame overhead duration : {detector.params.timings.frame_overhead_duration:~}
        frame overlap duration : {detector.params.timings.frame_overlap_duration:~} 

        ```
    """
    )
    return


@app.cell
def __(binning_on, channel, detector, mo):
    mo.md(
        f"""
        The derived detector parameters are given below for the detector with the name "_{detector.params.name}_" and channel "_{channel.name}_":
        
        ```
        Nyquist freq : {channel.nyquist_freq(False):~P.4} {f"({channel.nyquist_freq(True):~P.4} binned)" if binning_on else ''}
        number of pixels (full frame) : {detector.pixel_count:~P.4}
        number of pixels (full frame, channel) : {channel.pixel_count_frame(False):~P.4} {f"({channel.pixel_count_frame(True):~P.4} binned)" if binning_on else ''}
        ```
    """
    )
    return


@app.cell
def __(binning_on, channel, mo, optics):
    mo.md(
        f"""
        ### Imager Geometry Parameters
        
        Basic derived parameters are given below for the imager:
        
        ```
        ifov : {channel.ifov(optics, False):~P.4} {f"({channel.ifov(optics, True):~P.4} binned)" if binning_on else ''}
        pixel solid angle : {channel.pix_solid_angle(optics, False):~P.4} {f"({channel.pix_solid_angle(optics, True):~P.4} binned)" if binning_on else ''}
        horizontal full fov: {channel.horizontal_fov(optics):~P.4} 
        vertical full fov: {channel.vertical_fov(optics):~P.4}
        ```
    """
    )
    return


@app.cell
def __(mo):
    mo.md(
        f"""
        ### Geometric Projection Parameters

    """
    )
    return


@app.cell
def __(channel, distance, np, optics, u):
    # Ground sample distance at distance
    spatial_sample_distance_native = (
        (
            distance
            * (
                channel.pixel_pitch(False)
                * channel.binning
                / optics.params.focal_length
            )
        )
        .to_reduced_units()
        .to(u.m)
    )
    spatial_sample_distance = (
        (
            distance
            * (
                channel.pixel_pitch(True)
                * channel.binning
                / optics.params.focal_length
            )
        )
        .to_reduced_units()
        .to(u.m)
    )

    # image width and height assuming flat plate and constant Instantaneous FoV
    image_width = (
        2
        * np.tan(channel.ifov(optics, False) * channel.horizontal_pixels / 2.0)
        * distance
    )

    image_height = (
        2
        * np.tan(channel.ifov(optics, False) * channel.vertical_pixels / 2.0)
        * distance
    )
    return (
        image_height,
        image_width,
        spatial_sample_distance,
        spatial_sample_distance_native,
    )


@app.cell
def __(
    binning_on,
    distance,
    image_height,
    image_width,
    mo,
    spatial_sample_distance,
    spatial_sample_distance_native,
):
    mo.md(
        f"""
        ```
        target distance : {distance:~P}   
        spatial sample distance : {spatial_sample_distance_native:~P.4} {f"({spatial_sample_distance:~P.4} binned)" if binning_on else ''} (at target distance)
        image width : {image_width:~.4P} (at target distance)
        image height : {image_height:~.4P} (at target distance)
        ```
    """
    )
    return


@app.cell
def __(detector, mo):
    timings = detector.params.timings

    mo.md(
        f"""
    ### Timings

        ```
        frame rate :  {timings.frame_rate:~P.6} ({timings.frame_duration:~P.4})
        max integration duration : {timings.max_integration_duration:~P.4}
        actual integration duration : {timings.integration_duration:~P.4} 
        ```
    """
    )
    return timings,


@app.cell
def __(binning_on, channel, detector, mo, rw_electronics, timings):
    if rw_electronics:
        rw_electronics_text = f"""    
        ```
        pixel read rate : {channel.pix_read_rate(timings.frame_rate, False, False):~P.6} {f'({channel.pix_read_rate(timings.frame_rate, True, False):~P.6} binned)' if binning_on else ''} ({channel.name} channel)

        data write rate (uncompressed, incl. overheads) : {rw_electronics.data_write_rate(channel, detector, False, False):~P.6} {f'({rw_electronics.data_write_rate(channel, detector, True, False):~P.6} binned)' if binning_on else ''} ({channel.name} channel)
        data write rate (compressed, incl. overheads) : {rw_electronics.data_write_rate(channel, detector, False, True):~P.6} {f'({rw_electronics.data_write_rate(channel, detector, True, True):~P.6} binned)' if binning_on else ''} ({channel.name} channel)

        ```
        """
    else:
        rw_electronics_text = "No Read-out/Write Data loaded."

    mo.md(
        f"""
    ### Read-out/Write Electronics

        {rw_electronics_text}

    """
    )
    return rw_electronics_text,


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
