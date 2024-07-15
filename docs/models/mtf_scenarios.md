# MTF Scenarios

## Picking the Correct Scope of MTF Models

It is important to understand the problem and the physical sources that add blurring to the image. For each physical phenomenon, we have to quantify the extent, use a model and quantify it as a separate MTF item. All items relevant for the problem are then multiplied to compute the final MTF. Also, some of the MTF sources are directional and, particularly for line scanners, the along-track and across-track MTF should be computed separately.

Imager MTF (as tested on-ground):

- Optics (Aberrated)
- Detector Sampling
- Detector Diffusion
- Cryocoolers, shutters and similar vibration sources on the Imager

Platform MTF

- Imager MTF (see above)
- Target-to-Imager relative motion during exposure (motion blur)
- Target-to-Imager relative rotation during exposure (for line scanners only)
- Platform Stability during exposure (drift/smear due to attitude) - above and around the attitude control frequency
- Platform Vibration during exposure (drift/smear and jitter) - below the attitude control frequency (to be computed for each vibration source separately)
- Multiple image acquisition impacts (if applicable)

difficult to test on ground. Some items possible.

System MTF

- Platform MTF (see above)
- Atmospheric aerosol
- Atmospheric turbulence (where applicable)

The MTF computation is tailored for the problem at hand. For example the Hubble Space Telescope, only Platform MTF applies as there is no atmosphere for the images. An Earth Observation satellite may have practically all of them. An imager on a drone will also have practically all of them, but the imager is likely not a line scanner and atmospheric turbulence may not feature.

## Practical Example: High-Res Earth observation Satellite

As an example, for a high-resolution Earth Observation satellite operating in VisNIR, employing with a TDI-capable line scanner, the following blurring sources will require a separate MTF item in the overall MTF computation (or the MTF budget).

- Imager MTF
  - Optics (Aberrated)
  - Detector Sampling
  - Detector Diffusion

- Platform MTF
  - Imager MTF (see above)
  - Drift/smear due to ground (or satellite) motion during exposure
  - Drift/smear due to yaw steering error (or imperfect yaw steering control)
  - Drift/smear due to residual yaw steering error (on the pixel edges)
