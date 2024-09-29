# MTF Scenarios

## Picking the Correct Scope of MTF Models

It is important to understand the problem and the physical sources that add blurring to the image. For each physical phenomenon, we have to quantify the extent, use a model and quantify it as a separate MTF item. All items relevant for the problem are then multiplied to compute the final MTF. Also, some of the MTF sources are directional and, particularly for line scanners, the along-track and across-track MTF should be computed separately.

While it does depend on the problem, the usual scopes for the MTF computation can be Imager, Platform and End-to-End System.

Imager MTF (as tested in the lab):

- Optics (Aberrated)
- Detector Sampling
- Detector Diffusion
- Cryocoolers, shutters and similar vibration sources on the Imager

This is usually the requirement (or its breakdown) towards the Imager Supplier. Each item can be directly measured or at least derived via lab measurements. The total Imager MTF is the combination of these items.

Platform MTF

- Imager MTF (see above)
- Target-to-Imager relative motion during exposure (motion blur)
- Target-to-Imager relative rotation during exposure (for line scanners only)
- Platform Stability during exposure (drift/smear due to attitude) - above and around the attitude control frequency
- Platform Vibration during exposure (drift/smear and jitter) - below the attitude control frequency (to be computed for each vibration source separately)
- Multiple image acquisition impacts (if applicable)

This is usually the requirement (or its breakdown) towards the Platform Supplier (for example a satellite or aircraft system integrator). It is usually not possible to measure each item separately -or even in combination- in the lab. However, they can be modelled and quantified to keep track of the performance. They can be measured once the platform is undergoing tests in the operational conditions (a satellite in space or an aircraft in the air), using dedicated MTF targets.

While there are numerous sources for drift/smear and jitter, they are not necessarily handled with a large number of MTF models that are multiplied. Rather, their amplitude can be summed properly, taking into account the statistical properties and a single drift/smear magnitude (in along-track and across-track) as well as a single jitter magnitude (in along-track and across-track, though there usually the same) are computed. Then only two MTF contributions are needed per axis: one for the total drift/smear and one for the total jitter.

End-to-End System MTF

- Platform MTF (see above)
- Atmospheric aerosol
- Atmospheric turbulence (where applicable)

This is usually the requirement (or its breakdown) towards the System Supplier, usually defined (even if via proxy parameters) by the Customer or the Product Team. It is usually not possible to measure each item separately -or even in combination- in the lab. However, they can be modelled and quantified to keep track of the performance. They can be measured once the platform is undergoing tests in the operational conditions (a satellite in space or an aircraft in the air), using dedicated MTF targets.

The MTF computation is tailored for the problem at hand. For example the Hubble Space Telescope, only Platform MTF applies as there is no atmosphere for the images. An Earth Observation satellite may have practically all of them. An imager on a drone will also have practically all of them, but the imager is likely not a line scanner and atmospheric turbulence may not feature.

## Practical Example: High-Res Earth observation Satellite

As an example, for a high-resolution Earth Observation satellite operating in VNIR, employing with a TDI-capable line scanner, the following blurring sources will require a separate MTF item in the overall MTF computation (or the MTF budget).

- Imager MTF
  - Optics (Aberrated)
  - Detector Sampling
  - Detector Diffusion
  - No drift/smear or jitter due to cryocoolers or shutter mechanisms assumed
  
- Platform MTF
  - Imager MTF (see above)
  - Drift/smear due to ground (or satellite) motion during exposure (motion blur)
  - Drift/smear due to yaw steering error (or imperfect yaw steering control) during exposure and during a full TDI column
  - Drift/smear due to roll, pitch or yaw rate (or imperfect attitude rate control) during exposure and during a full TDI column
  - Drift/smear due to residual yaw steering error (on the detector edge pixels) during exposure and during a full TDI column
  - Jitter due to high frequency platform vibrations (e.g., reaction wheels) during exposure and during a full TDI column
  - Drift/smear due to low frequency platform vibrations (e.g., flexible body dynamics or propellant sloshing) during exposure and during a full TDI column

- End-to-End System MTF
  - Platform MTF (see above)
  - Atmospheric aerosol
  - Atmospheric turbulence

  It should be noted that, many of the Platform MTF contributors for a single exposure (or a single image acquisition of the TDI column) are small enough to be ignored for a high resolution satellite. A 70 cm GSD would correspond to about 0.1 ms line rate (or maximum theoretical exposure duration), during which many of the dynamics are too slow to be of significant concern. Nevertheless, a quick check is recommended for all the items, either in terms of quantity (e.g., drift rate times the exposure duration in terms of pixels) or frequency, where the frequency of the disturbance is significantly smaller than the line rate frequency and it acts not as a drift but as a bias, without introducing a blurring. In this example with the 70 cm GSD and 10 stage-TDI, a single acquisition is at 10000 Hz (0.1 ms) and the total TDI column is at 1000 Hz (1 ms). This is way too fast for most "low frequency-high amplitude" disturbances to cause an appreciable blurring. The higher frequency effects (primarily high acceleration jerks) may cause blurring, and only the "very high frequency-very low amplitude" jitter will feature in the images.

  A full example is given [here](docs/examples/sat_mtf_budget.ipynb).
