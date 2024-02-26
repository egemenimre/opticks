import marimo

__generated_with = "0.2.7"
app = marimo.App()


@app.cell
def __():
    import marimo as mo

    mo.md("# Basic Pushbroom Imager on a Satellite")
    return (mo,)


@app.cell
def __(mo):
    mo.md(
        f"""
        ## Setting the Scene

        Computing the basic parameters of a satellite pushbroom imager. 
        """
    )
    return


@app.cell
def __():
    from opticks import u
    import numpy as np

    # constants
    # ---------
    r_earth = 6378.137 * u.km
    mu = 398600.5 * u.km**3 / u.s**2


    # sat positional params
    # ---------------------
    sat_altitude = 540.0 * u.km

    n = np.sqrt(mu / (r_earth + sat_altitude) ** 3)
    ground_vel = (n * r_earth).to("m/s")

    # converted to generic params
    distance = sat_altitude
    target_rel_velocity = ground_vel
    return (
        distance,
        ground_vel,
        mu,
        n,
        np,
        r_earth,
        sat_altitude,
        target_rel_velocity,
        u,
    )


@app.cell
def __(ground_vel, mo, sat_altitude):
    mo.md(
        f"""
        Reference parameters for the camera position and motion with respect to the target are given below:

        ```
        sat altitude : {sat_altitude:~P}    
        ground velocity : {ground_vel:~P}
        ```
        """
    )
    return


@app.cell
def __(mo):
    # init the ref wavelength slider in a separate cell
    slider = mo.ui.slider(450, 2500, label="reference wavelength")

    mo.md(
        f"""
        Certain parameters vary with the light characteristics and therefore wavelength dependent.
        """
    )
    return (slider,)


@app.cell
def __(mo, slider, u):
    ref_wavelength = slider.value * u.nm

    mo.md(
        f"""
        Set the reference wavelength with the slider:

        {slider} 
        ```
        reference wavelength : {ref_wavelength:~P}    
        ```

        Note the widely used bands:

        - blue: 450-485 nm
        - green: 500-565 nm
        - red: 625-750 nm
        - nir: 780-2500 nm
        """
    )
    return (ref_wavelength,)


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
        file_directory = Path("docs", "examples", "sample_sat_pushbroom")
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
        rw_electronics_yaml = None

        if rw_electronics_file.name(0):
            rw_electronics_yaml = rw_electronics_file.contents(0).decode("utf-8")

        # Init imager object
        imager = Imager.from_yaml_text(optics_yaml, detector_yaml, rw_electronics_yaml)

    # shorthands
    optics = imager.optics
    detector = imager.detector
    rw_electronics = imager.rw_electronics

    # binning status
    binning_on = False if detector.params.binning == 1 else True
    return (
        Imager,
        Path,
        binning_on,
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
        aperture area : {optics.aperture_area:~P.6}
        spatial cutoff freq  : {optics.spatial_cutoff_freq(ref_wavelength):~P.5} (at {ref_wavelength:~P})
        ```
    """
    )
    return


@app.cell
def __(binning_on, detector, mo):
    mo.md(
        f"""
        ### Detector Parameters

        Basic design parameters are given below for the detector with the name "_{detector.params.name}_":

        ```
        detector type : {detector.params.detector_type}
        horizontal x vertical pixels (total) : {detector.params.horizontal_pixels} x {detector.params.vertical_pixels}
        horizontal x vertical pixels (used) : {detector.params.horizontal_pixels_used} x {detector.params.vertical_pixels_used}
        binning : {detector.params.binning if binning_on else 'None'}
        pixel pitch : {detector.pix_pitch(False):~P.4} {f"({detector.pix_pitch(True):~P.4} binned)" if binning_on else ''}
        TDI stages : {'None' if detector.params.tdi_stages == 1 else detector.params.tdi_stages}

        Timing parameters:
        integration duration : {detector.params.timings.integration_duration:~}  
        frame overhead duration : {detector.params.timings.frame_overhead_duration:~}  
        frame overlap duration : {detector.params.timings.frame_overlap_duration:~}   
        ```
    """
    )
    return


@app.cell
def __(binning_on, detector, mo):
    mo.md(
        f"""
        The derived detector parameters are given below for the detector with the name "_{detector.params.name}_":

        ```
        Nyquist freq : {detector.nyquist_freq(False):~P.4} {f"({detector.nyquist_freq(True):~P.4} binned)" if binning_on else ''}
        number of pixels (full frame) : {detector.pixel_count_frame(False):~P.4} {f"({detector.pixel_count_frame(True):~P.4} binned)" if binning_on else ''}
        number of pixels (full frame, used) : {detector.pixel_count_frame(False,True):~P.4} {f"({detector.pixel_count_frame(True,True):~P.4} binned)" if binning_on else ''}
        ```
    """
    )
    return


@app.cell
def __(binning_on, imager, mo):
    mo.md(
        f"""
        ### Imager Geometry Parameters

        Basic derived parameters are given below for the imager:

        ```
        ifov : {imager.ifov(False):~P.4} {f"({imager.ifov(True):~P.4} binned)" if binning_on else ''}
        pixel solid angle : {imager.pix_solid_angle(False):~P.4} {f"({imager.pix_solid_angle(True):~P.4} binned)" if binning_on else ''}
        horizontal full fov: {imager.horizontal_fov:~P.4} (used pixels only)
        vertical full fov: {imager.vertical_fov:~P.4} (used pixels only)
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
def __(distance, imager, np, u):
    # Ground sample distance at nadir
    spatial_sample_distance_native = (
        (
            distance
            * (
                imager.detector.pix_pitch(False)
                * imager.detector.params.binning
                / imager.optics.params.focal_length
            )
        )
        .to_reduced_units()
        .to(u.m)
    )
    spatial_sample_distance = (
        (
            distance
            * (
                imager.detector.pix_pitch(True)
                * imager.detector.params.binning
                / imager.optics.params.focal_length
            )
        )
        .to_reduced_units()
        .to(u.m)
    )

    # swath assuming flat plate and constant Instantaneous FoV
    swath = (
        2
        * np.tan(
            imager.ifov(False) * imager.detector.params.horizontal_pixels_used / 2.0
        )
        * distance
    )
    return spatial_sample_distance, spatial_sample_distance_native, swath


@app.cell
def __(
    binning_on,
    mo,
    spatial_sample_distance,
    spatial_sample_distance_native,
    swath,
):
    mo.md(
        f"""
        ```
        spatial sample distance : {spatial_sample_distance_native:~P.4} {f"({spatial_sample_distance:~P.4} binned)" if binning_on else ''}
        swath : {swath:~.4P}  (disregarding Earth curvature)
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
        line duration : {timings.frame_duration:~P.4} 
        line rate :  {timings.frame_rate:~P.6} 
        max integration duration : {timings.max_integration_duration:~P.4}
        actual integration duration : {timings.integration_duration:~P.4} 
        total TDI column duration : {timings.total_tdi_col_duration:~P.4} ({detector.params.tdi_stages}x stages)   
        ```
    """
    )
    return (timings,)


@app.cell
def __(binning_on, detector, mo, rw_electronics):
    if rw_electronics:
        rw_electronics_text = f"""    
        ```
        pixel read rate (without TDI) : {detector.pix_read_rate(False, False):~P.4} {f'({detector.pix_read_rate(True, False):~P.4} binned)' if binning_on else ''}
        pixel read rate (with TDI) : {detector.pix_read_rate(False, True):~P.4} {f'({detector.pix_read_rate(True, True):~P.4} binned)' if binning_on else ''}

        data write rate (uncompressed, incl. overheads) : {rw_electronics.data_write_rate(detector, False, False):~P.4} {f'({detector.pix_read_rate(True, False):~P.4} binned)' if binning_on else ''}
        data write rate (compressed, incl. overheads) : {rw_electronics.data_write_rate(detector, False, True):~P.4} {f'({detector.pix_read_rate(True, True):~P.4} binned)' if binning_on else ''}
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
    return rw_electronics_text


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
