# Computing PSF and MTF with `opticks`

*For a primer on Modulation Transfer Function (MTF), its contributors and how they are combined, you may want to have a look at [the discussion on sharpness](sharpness_pt1.md) and [the correct contributors for your scenario](mtf_scenarios.md)*

## Generating Basic MTF Models

The Modulation Transfer Function (MTF) is useful to define the sharpness of the system. In most cases, multiple MTF contributors are computed and then combined to yield the "System MTF". A full example to build a System MTF is [given here](../examples/sat_mtf_budget.ipynb).

In `opticks`, each MTF contributor is modelled by an {py:class}`.MTF_Model_1D` object. It should be noted that, in reality the MTF is a 2D surface defined for each wavelength (monochromatic) or a combination of wavelengths (polychromatic) and also for each point on the image plane. This {py:class}`.MTF_Model_1D` definition is just a 1D slice in x or y directions, defined for a single point on the image plane. Therefore the optical MTF will be different for the centre of the optical axis or at the edge.

To put {py:class}`.MTF_Model_1D` into use, the following simply generates an Ideal Optical System MTF with a given {py:class}`.Optics` model and a wavelength (as the MTF is wavelength dependent):

```python
from opticks import u
from opticks.imaging_model.optics import Optics

optics = Optics.from_yaml_file("optics_file.yaml")
wvl = 650 * u.nm

mtf_model = MTF_Model_1D.ideal_optics(wvl, optics)
```

The following MTF types are supported:

- Ideal Optics ({py:meth}`.MTF_Model_1D.ideal_optics`): Circular aperture, no aberrations
- External Data ({py:meth}`.MTF_Model_1D.external_data`): External discrete data with interpolation
- Aberrated Optics with Empirical Model ({py:meth}`.MTF_Model_1D.emp_model_aberrated_optics`): Aberrated optics with an empirical model for aberrations
- Slicing an External 2D MTF ({py:meth}`.MTF_Model_1D.from_mtf_2d`): This is in a sense similar to External Data, but acts as the bridge between more detailed (aberrated) optical models from `prysm`
- Detector Sampling ({py:meth}`.MTF_Model_1D.detector_sampling`): Detector sampling MTF model
- Detector Diffusion ({py:meth}`.MTF_Model_1D.detector_diffusion`): Carrier diffusion MTF model (Crowell & Labuda 1969), with five geometry/boundary variants
- Detector Diffusion Preset ({py:meth}`.MTF_Model_1D.detector_diffusion_preset`): Detector diffusion MTF from a predefined detector category (BSI CCD, BSI sCMOS, BSI IR FPA, FSI CMOS, FSI CCD)
- Detector Crosstalk ({py:meth}`.MTF_Model_1D.detector_crosstalk`): Nearest-neighbour crosstalk MTF model with separate side and diagonal coupling coefficients
- Detector CTE ({py:meth}`.MTF_Model_1D.detector_cte`): Charge transfer efficiency MTF for CCDs (direction-dependent — instantiate one model per direction)
- Motion Blur ({py:meth}`.MTF_Model_1D.motion_blur`): Motion blur MTF model
- Drift/Smear ({py:meth}`.MTF_Model_1D.smear`): Drift/Smear MTF model, the more generalised form of motion blur
- Jitter ({py:meth}`.MTF_Model_1D.jitter`): Jitter MTF model
- Resampling ({py:meth}`.MTF_Model_1D.resampling`): Post-acquisition resampling MTF (geometric correction / orthorectification), with five reconstruction kernels (nearest-neighbour, bilinear, bicubic, Lanczos, sinc)
- Combined ({py:meth}`.MTF_Model_1D.combined`): The MTF model that combines any of the items above, used to group or sum them together
- Fixed ({py:meth}`.MTF_Model_1D.fixed`): MTF model for a fixed MTF value for the entire spatial domain, used for crude modelling of disturbances and imperfections

The full [System MTF example](../examples/sat_mtf_budget.ipynb) illustrates the use of almost all of them. "Slicing an External 2D MTF" is explained in the next section, as it enables detailed Wavefront models to be incorporated.

## Detector Diffusion, Crosstalk, and CTE MTF

### Detector Diffusion MTF

The detector diffusion MTF models lateral carrier diffusion in the semiconductor substrate. The underlying physics and model derivation are explained in [Detector MTF](sharpness_pt2_det.md). Five model variants are available, corresponding to different detector geometries and boundary conditions: BSI-1 (scientific CCD), BSI-2 (sCMOS), BSI-3 (IR FPA), FSI-1 (consumer CMOS), and FSI-2 (bulk CCD). These map to the `DetectorDiffusionModel` enum.

For full control over model parameters, use {py:meth}`.MTF_Model_1D.detector_diffusion` directly:

```python
from opticks import u
from opticks.contrast_model.mtf import MTF_Model_1D
from opticks.contrast_model.detector_mtf import DetectorDiffusionModel

# Silicon absorption coefficient at 850 nm (~NIR)
alpha = 0.065 / u.um

# BSI-2 (sCMOS): field-free depth L_a and depletion depth L_b required
mtf_diffusion = MTF_Model_1D.detector_diffusion(
    model=DetectorDiffusionModel.BSI_2,
    absorption_coeff=alpha,
    diffusion_length=1.5 * u.um,
    field_free_depth=5.0 * u.um,
    depletion_depth=1.5 * u.um,
)
```

For common detector categories, use {py:meth}`.MTF_Model_1D.detector_diffusion_preset`, which looks up default parameters and allows overrides:

```python
from opticks.contrast_model.detector_mtf import DetectorDiffusionPreset

# FSI CCD preset with a custom diffusion length override
mtf_diffusion = MTF_Model_1D.detector_diffusion_preset(
    preset=DetectorDiffusionPreset.FSI_CCD,
    absorption_coeff=0.21 / u.um,      # red channel (~650 nm)
    diffusion_length=3.0 * u.um,       # override preset default
)
```

The absorption coefficient is wavelength-dependent (shorter wavelengths are absorbed nearer the surface, longer wavelengths penetrate deeper). The diffusion MTF should therefore be computed separately for each spectral channel. A [model derivation notebook](misc/detector_diffusion_derivation.ipynb) and a [numerical examples notebook](misc/detector_diffusion_mtf.ipynb) with real detector parameters are available.

### Detector Crosstalk MTF

The crosstalk MTF models nearest-neighbour charge sharing using a discrete kernel with side ($X_s$) and diagonal ($X_d$) coupling coefficients. The 1D MTF slice is:

$$\text{MTF}_\text{xtalk}(f) = 1 - 2(X_s + 2X_d)(1 - \cos(2\pi f p))$$

```python
# 3% side crosstalk, 0.5% diagonal crosstalk, 5 µm pitch
mtf_crosstalk = MTF_Model_1D.detector_crosstalk(
    pixel_pitch=5 * u.um,
    crosstalk_xs=0.03,
    crosstalk_xd=0.005,
)
```

When only side crosstalk is relevant (the common case), omit `crosstalk_xd` or pass 0.

```{note}
The diffusion MTF and the electrical component of crosstalk describe the same underlying physics. Do not multiply them in an MTF budget unless the crosstalk coefficient refers to optical crosstalk only. See [Detector MTF](sharpness_pt2_det.md) for guidance.
```

### Detector CTE MTF

The CTE MTF models charge transfer inefficiency in CCD detectors. During readout, charge is clocked through neighbouring pixels toward the output amplifier; each transfer leaves a small fraction behind in silicon traps. The effect is cumulative and **strictly directional** — it produces an asymmetric tail on every point source pointing away from the readout register. Unlike diffusion (isotropic) or aperture (symmetric), CTE blur applies only along the clocking axis.

$$\text{MTF}_\text{CTE}(f) = \exp\!\left(-N(1-\varepsilon)\!\left(1 - \cos\!\left(\frac{\pi f}{f_\text{nyq}}\right)\right)\right)$$

where $N = (\text{num\_pixels} + \text{tdi\_stages}) \times \text{num\_phases}$, $\varepsilon$ = CTE per transfer, and $f_\text{nyq} = 1/(2p)$.

**CTE MTF is not constant across the image** — it degrades roughly linearly from pixels nearest the readout (sharp) to those at the far end (blurry). This is why `num_pixels` is a runtime argument, not a fixed sensor property: instantiate one model per pixel position of interest.

```python
from opticks import u
from opticks.contrast_model.mtf import MTF_Model_1D

pitch = 13 * u.um

# Non-TDI 3-phase CCD: ALT (parallel register), worst-case pixel
mtf_cte_alt = MTF_Model_1D.detector_cte(
    pixel_pitch=pitch,
    num_pixels=4096,   # rows from pixel to readout edge
    num_phases=3,      # 3-phase CCD
    cte=0.99999,       # high-quality CTE (not CTI — pass 0.99999, not 0.00001)
)

# TDI CCD with 32 on-chip analog stages: ALT only — each stage adds transfers
mtf_cte_alt_tdi = MTF_Model_1D.detector_cte(
    pixel_pitch=pitch,
    num_pixels=4096,
    num_phases=3,
    cte=0.99999,
    tdi_stages=32,     # analog TDI stages in the parallel register
)

# ACT (across-track) MTF — serial register only, no TDI
mtf_cte_act = MTF_Model_1D.detector_cte(
    pixel_pitch=pitch,
    num_pixels=6000,   # columns from pixel to readout amplifier
    num_phases=3,
    cte=0.99999,
)
```

```{note}
`tdi_stages` applies only to **on-chip analog TDI** in CCDs, where each stage is a real charge transfer in the parallel register. For **digital ("shift-and-add") TDI** — CMOS sensors that read out each row independently and accumulate in firmware — set `tdi_stages = 0`, since no charge transfer occurs between stages. TDI is also parallel-register only, so always pass `tdi_stages = 0` for the ACT direction.
```

**Radiation aging.** Space radiation creates traps in the silicon lattice that capture and re-emit electrons, increasing CTI ($1 - \varepsilon$) over the mission lifetime. EO satellite MTF budgets should use the predicted End-of-Life (EOL) CTE, not the launch (BOL) value.

| Condition | CTE ($\varepsilon$) | Inefficiency $1-\varepsilon$ | MTF at Nyquist ($N=2000$) |
| --- | --- | --- | --- |
| High quality (BOL) | $0.999999$ | $10^{-6}$ | $\approx 0.996$ (negligible) |
| Degraded (EOL, post-radiation) | $0.99995$ | $5 \times 10^{-5}$ | $\approx 0.818$ (severe) |

## Point Spread Functions and Optics MTF

The MTF associated with real optical systems are more complex and are usually defined via the Point Spread Function (PSF).

For this, [`prysm`](https://github.com/brandondube/prysm/) is used as the Wavefront calculation backend and some wrappers provided by `opticks`. The methodology is as follows:

1. We define an aperture (conveniently through the {py:class}`.Aperture`), though `prysm` can also be used. See the [tutorial here](../tutorials/aperture.ipynb) for an example on defining a complex aperture.
2. We then attach the aperture to an {py:class}`.Optics` object via {py:meth}`.Optics.set_aperture_model`.
3. For each wavelength, we form one Pupil Function (essentially a Wavefront with an Amplitude and a Phase Error) and attach it to the {py:class}`.Optics` object via {py:meth}`.Optics.add_mono_pupil_function` or {py:meth}`.Optics.add_poly_pupil_function`.
4. We call {py:meth}`.Optics.compute_psf` to compute the PSF. The result is a `prysm` `RichData` object. It is possible to plot this 2D PSF using the `RichData.plot2d` method.
5. We call {py:meth}`.Optics.mtf` to retrieve the 2D MTF (lazily computed from the cached PSF). The result is also a `prysm` `RichData` object.
6. Finally we call {py:meth}`.Optics.to_MTF_Model_1D` (or {py:meth}`.MTF_Model_1D.from_mtf_2d`) with the direction of the slice, for example in "x" or "y" direction. This yields the 1D MTF Model.

This process is illustrated in an [example notebook](../examples/optics_aperture_psf.ipynb).

## Off-Axis / Field-Dependent MTF

The PSF and MTF computed above apply to a single point on the image plane. Off-axis pixels see different — and generally worse — wavefront errors due to field-dependent aberrations (coma, astigmatism, field curvature). The [theory and modelling tiers](sharpness_pt2_opt.md#field-dependence-of-the-optical-mtf) are discussed in the optics MTF documentation.

`opticks` supports three paths for computing field-dependent MTF:

**Tier A — Measurements / ray-trace data.** Ingest per-field MTF curves via {py:meth}`.MTF_Model_1D.external_data`, or per-field Zernike coefficients via {py:meth}`.OptPathDiff.from_zernike` and {py:meth}`.Optics.add_mono_pupil_function`. The user registers one pupil function per field point using a naming convention like `f"field_{hx:.2f}_{hy:.2f}"` — `Optics.pupil_functs` already supports multiple named entries.

**Tier B — Analytical Seidel + prysm PSF.** Use {py:meth}`.FieldAberrationModel.to_zernikes` to get the Zernike vector at any field point, then follow the standard `OptPathDiff.from_zernike` → `add_mono_pupil_function` → `compute_psf` pipeline.

**Tier C — Analytical Seidel + empirical ATF (fastest).** Use {py:meth}`.Optics.field_mtf_model_1d` for a one-call path that computes $W_{\mathrm{rms}}(H)$ from the Seidel model and feeds it into the ATF:

```python
from opticks import u
from opticks.contrast_model.optics_mtf import FieldAberrationModel

# define Seidel coefficients (in waves, at the edge of the field)
field_model = FieldAberrationModel.from_waves(
    w040=0.02, w131_edge=0.08, w222_edge=0.05, w220_edge=0.03,
    reference_wavelength=650 * u.nm,
)

wvl = 650 * u.nm

# MTF at on-axis (h=0), mid-field (h=0.5), and edge (h=1.0)
mtf_on_axis = optics.field_mtf_model_1d(field_model, h=0.0, wavelength=wvl)
mtf_mid     = optics.field_mtf_model_1d(field_model, h=0.5, wavelength=wvl)
mtf_edge    = optics.field_mtf_model_1d(field_model, h=1.0, wavelength=wvl)
```

## Resampling MTF

The resampling MTF models the post-acquisition reconstruction step (typically orthorectification onto a fixed ground grid). The underlying physics, the $p_\text{eff} = \max(p_\text{in}, p_\text{out})$ rule, and per-kernel formulas are in [Processing MTF](sharpness_pt5.md). Five kernels are available via the `ResamplingKernel` enum: `NEAREST_NEIGHBOR`, `BILINEAR`, `BICUBIC`, `LANCZOS`, `SINC`.

For a one-shot model, use {py:meth}`.MTF_Model_1D.resampling`:

```python
from opticks import u
from opticks.contrast_model.mtf import MTF_Model_1D
from opticks.contrast_model.processing_mtf import ResamplingKernel

# Bilinear, downsample (input SSD 0.5 m onto a 1 m ortho grid).
# p_eff = max(0.5 m, 1.0 m) = 1.0 m  →  the kernel acts as the AA prefilter.
mtf_resampling = MTF_Model_1D.resampling(
    kernel=ResamplingKernel.BILINEAR,
    input_pitch=0.5 * u.m,
    output_pitch=1.0 * u.m,
)

# Bicubic with the standard Keys parameter (a = -0.5)
mtf_bicubic = MTF_Model_1D.resampling(
    kernel="bicubic",          # string also accepted
    input_pitch=2.0 * u.m,     # upsample: SSD wider than ortho grid
    output_pitch=1.0 * u.m,
    bicubic_a=-0.5,
)

# Lanczos-3
mtf_lanczos = MTF_Model_1D.resampling(
    kernel=ResamplingKernel.LANCZOS,
    input_pitch=1.0 * u.m,
    output_pitch=1.0 * u.m,
    lanczos_n=3,
)
```

For pipeline use, the {py:class}`.Processing` component carries the kernel choice and the (default) ortho grid pitch in YAML:

```yaml
# processing.yaml
name: Sample Processing
processing_params:
  resampling_kernel: bilinear
  output_pitch: 1 m
```

```python
from pathlib import Path
from opticks import u
from opticks.imaging_model.processing import Processing

processing = Processing.from_yaml_file(Path("processing.yaml"))

# Sweep input_pitch across the swath to map SSD-driven MTF variation:
for ssd in (0.5, 0.7, 1.0, 1.5, 2.0) * u.m:
    mtf = processing.get_resampling_mtf_1d(input_pitch=ssd)
    # … combine with optics/detector/motion via MTF_Model_1D.combined(...)

# output_pitch can be overridden per call (e.g. to compare grid choices)
mtf_2m_grid = processing.get_resampling_mtf_1d(
    input_pitch=0.5 * u.m,
    output_pitch=2.0 * u.m,
)
```

The full imaging pipeline — instrument plus processing — is held by an {py:class}`.ImagingChain`:

```python
from opticks.imaging_model.imager import Imager
from opticks.imaging_model.imaging_chain import ImagingChain

imager = Imager.from_yaml_file(
    "optics.yaml", "pan_detector.yaml", "rw_electronics.yaml"
)
processing = Processing.from_yaml_file("processing.yaml")

chain = ImagingChain(imager, processing)

# chain.imager has the hardware (FOV, GSD, data rates, hardware MTFs);
# chain.processing has the post-acquisition MTFs.
```

`ImagingChain` is a pure container — there is no `from_yaml_file` factory. Build the components separately and compose. This keeps `Imager` focused on the hardware physics and lets future post-acquisition stages (sharpening, denoising) slot in next to `Processing` without changing the imager.

```{note}
Bicubic and Lanczos can produce **MTF $> 1$** in the mid-passband (real edge-boost behaviour, not a numerical artefact). For ML / detection inputs prefer Bilinear, which has no overshoot. See the [ringing caveat](sharpness_pt5.md#ringing-and-edge-boost) for the trade-off.
```

## Combining MTF Models

Usually more than one MTF contributor is present, and we sum them properly to either group them (for example Imager MTF or Dynamic MTF) or eventually to combine all of them (End-to-End or System MTF).

The process is detailed for [a Satellite Imager here](../examples/sat_mtf_budget.ipynb). Essentially, all MTF sources are initialised via various flavours of the {py:class}`.MTF_Model_1D` object classmethods, and they are then combined via the {py:meth}`.MTF_Model_1D.combined` method.

For example the following is the final step of the MTF computation in the along-track and across-track directions, combining imager, single image dynamic and TDI MTF contributions.

```python
# compute the total system MTF in the ACT and ALT directions
mtf_model_act_sys = MTF_Model_1D.combined(
    mtf_model_imager_mid, mtf_model_single_act, mtf_model_tdi_act
)
mtf_model_alt_sys = MTF_Model_1D.combined(
    mtf_model_imager_mid, mtf_model_single_alt, mtf_model_tdi_alt
)
```

## Plotting MTF Curves

The {py:class}`.MTF_Model_1D` provides some plotting functionality for the classic MTF curves. The following combines multiple MTF curves.

```python
# Select midrange wavelength
ref_wavelength = channel.centre_wavelength
cutoff_freq = optics.spatial_cutoff_freq(ref_wavelength)
input_line_freq = np.linspace(0.0, cutoff_freq, 100)

# plot the MTF data
# -----------------
imager_mtf_plot = MTF_Plot_1D()

mtf_plot_list = {
    "Detector": mtf_det_sampling_model,  # detector MTF
    f"Ideal Optical {ref_wavelength_mid:~P}": mtf_ideal_opt_mid,  # Ideal Optical MTF
    f"Aberrated Optical {ref_wavelength_mid:~P}": mtf_aberr_opt_mid,  # Aberrated Optical MTF
    f"Imager {ref_wavelength_mid:~P}": mtf_model_imager_mid,  # total imager MTF
}

# populate the mtf plot with the models
imager_mtf_plot.populate_plot(input_line_freq, mtf_plot_list, nyq_limit=nyq_freq)

# set MTF plot style tweaks
imager_mtf_plot.set_plot_style(x_max=input_line_freq.max(), title="Imager MTF")
```

We begin by initialising the {py:class}`.MTF_Plot_1D`. We then add each item to a dict, where the key is the label of the MTF curve and the value is the {py:class}`.MTF_Model_1D` model. To populate the plot, we also need the discrete x axis values (`input_line_freq` here), limited by the spatial cut-off frequency. We can also plot a vertical "Nyquist Frequency" (`nyq_limit` key). The final step is to set the plot style with decorators like labels and titles, as well as limits and sizes. The keys are passed on to the underlying `matplotlib` implementation.

The end result looks like this:

![Static MTF](images/static_mtf.png "Sample MTF plot")
