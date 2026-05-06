# Processing Contributors to the MTF

The topic of sharpness and contrast performance is [continued from here](sharpness_pt4.md).

## Overview

The MTF contributors discussed in [Detector MTF](sharpness_pt2_det.md), [Dynamic MTF](sharpness_pt3.md), and [Atmospheric MTF](sharpness_pt4.md) all act on the light path *up to* the moment the detector finishes integrating charge. Once the imager has produced raw samples on its physical pixel grid, additional processing is applied — typically geometric correction (orthorectification onto a fixed ground grid), and possibly sharpening or denoising. The *resampling* step alone introduces a real, measurable MTF loss that belongs in the system MTF budget.

This page covers the resampling MTF. Sharpening and other processing-stage contributors are not yet implemented in `opticks`.

## Position in the System MTF Chain

The resampling MTF is the **final multiplier** in the system chain:

$$\text{MTF}_\text{system}(f) = \text{MTF}_\text{optics}(f) \cdot \text{MTF}_\text{detector}(f) \cdot \text{MTF}_\text{motion}(f) \cdot \text{MTF}_\text{atmos}(f) \cdot \text{MTF}_\text{resampling}(f)$$

It is applied on the *output* of the imager — meaning it operates on samples already band-limited and shaped by everything before it. It is not "an optical effect" in the usual sense; it is the transfer function of the discrete reconstruction kernel used to interpolate the input samples onto the output grid.

## Resampling MTF

### Geometric-Correction Use Case

The driving use case in `opticks` is **orthorectification** of pushbroom or framing imagery onto a fixed ground grid. Two pitches matter:

- **`input_pitch`** = the local **Spatial Sample Distance** (SSD) on the ground at the pixel of interest. The detector pixel pitch is constant, but the projected SSD varies across the image due to viewing geometry (off-nadir angle, terrain, platform attitude). For example, a satellite imaging at 30° off-nadir produces SSDs that change by tens of percent from one side of the swath to the other.

- **`output_pitch`** = the **fixed ortho grid pitch** chosen by the processing chain (e.g. 1 m for a Level-1B product). It is constant for the entire delivered image.

Resampling onto the ortho grid is therefore **upsample** in regions where SSD > output pitch (the kernel must bridge wide input gaps) and **downsample** where SSD < output pitch (the kernel must act as an anti-alias prefilter before decimation). Both regimes appear in the same image and must be modelled with the same kernel.

### The $p_\text{eff} = \max(p_\text{in}, p_\text{out})$ Rule

The standard treatment (Holst[^1] Ch. 10; Schowengerdt[^2]) is to scale the kernel by the **larger** of the two pitches:

$$p_\text{eff} = \max(p_\text{in}, p_\text{out})$$

The resampling MTF is then the magnitude of the Fourier transform of the kernel evaluated at the dimensionless frequency $\nu = f \cdot p_\text{eff}$:

$$\text{MTF}_\text{resampling}(f) = \left| H_\text{kernel}(\nu) \right|, \qquad \nu = f \cdot p_\text{eff}$$

This single rule covers both regimes:

- **Upsample** ($p_\text{in} > p_\text{out}$): $p_\text{eff} = p_\text{in}$. The kernel is scaled to the (wider) input grid so it bridges the gaps between input samples — this is the "reconstruction" role.
- **Downsample** ($p_\text{in} < p_\text{out}$): $p_\text{eff} = p_\text{out}$. The kernel is scaled to the (wider) output grid so its first zero falls near the output Nyquist — this is the "anti-alias prefilter" role.
- **Same pitch** ($p_\text{in} = p_\text{out}$): the kernel just smooths.

A consequence of the `max` rule is that **the kernel is the only resampling-stage MTF term** — there is no separate explicit anti-alias prefilter. The kernel scaled to $p_\text{out}$ on downsample already plays that role. Adding a second prefilter MTF would double-count the same effect.

### Per-Kernel Formulas

`opticks` implements five reconstruction kernels via the `ResamplingKernel` enum. Let $\text{sinc}(x) = \sin(\pi x) / (\pi x)$ (the [normalised sinc function](https://en.wikipedia.org/wiki/Sinc_function)). All formulas are normalised so that $\text{MTF}(0) = 1$.

#### Nearest-Neighbour (Pixel Replication)

$$h(x) = \text{rect}(x), \qquad \text{MTF}(\nu) = \left|\text{sinc}(\nu)\right|$$

Zero-order hold. Cheapest kernel; produces the most high-frequency residual (the sinc roll-off is the slowest of the family). Use only for pre-visualisation, never for analytical workflows.

#### Bilinear (Linear Interpolation)

$$h(x) = \text{tri}(x), \qquad \text{MTF}(\nu) = \text{sinc}^2(\nu)$$

Triangular kernel; the squared sinc is just the convolution of two box functions. Always non-negative — bilinear cannot produce contrast reversal. The most aggressive roll-off below Nyquist of the practical kernels.

#### Bicubic (Keys 1981[^3] — parameter $a$, default $-0.5$)

$$h(x) = \begin{cases} (a+2)|x|^3 - (a+3)|x|^2 + 1, & |x| \leq 1 \\ a|x|^3 - 5a|x|^2 + 8a|x| - 4a, & 1 < |x| \leq 2 \\ 0, & \text{otherwise} \end{cases}$$

The MTF is the (numerical) Fourier transform of $h(x)$, normalised by $H(0)$. The $a$ parameter shapes the boost/roll-off trade-off; $a = -0.5$ (the Keys cubic) is the standard photographic choice. Allowed range: $a \in [-1.0, 0.0]$.

#### Lanczos-N (Windowed Sinc — parameter $N$, default $3$)

$$h(x) = \begin{cases} \text{sinc}(x) \cdot \text{sinc}(x/N), & |x| < N \\ 0, & \text{otherwise} \end{cases}$$

Windowed sinc with $2N$ taps. The MTF is again the (numerical) Fourier transform normalised by $H(0)$. $N \geq 2$ (integer); $N = 3$ is the common default. Lanczos has the flattest passband of the practical kernels — the price is mild edge-boost behaviour (see caveat below).

#### Sinc (Ideal / Brick-wall)

$$h(x) = \text{sinc}(x), \qquad \text{MTF}(\nu) = \text{rect}(\nu)$$

The ideal reconstruction filter: $\text{MTF} = 1$ for $|\nu| < \tfrac{1}{2}$ and $\text{MTF} = 0$ for $|\nu| > \tfrac{1}{2}$. Mathematically optimal but not realisable — the kernel has infinite support and produces severe ringing. Useful as a theoretical reference, not as an actual processing choice.

### Ringing and Edge-Boost

Two of the five kernels can produce **MTF $> 1$** in the mid-passband:

- **Bicubic** with negative $a$ (the standard case) lifts the response above unity around $\nu \approx 0.3$–$0.4$.
- **Lanczos** likewise overshoots in the mid-band, often more strongly than Bicubic.

This is **not a numerical artefact** — it reflects the real edge-enhancement behaviour of these kernels. In the spatial domain it manifests as halos around sharp edges (overshoot and undershoot at step transitions). In an MTF chain it means the resampling stage can *boost* contrast at certain frequencies, partially compensating the roll-off from earlier stages.

```{note}
For ML inference pipelines, target detection workflows, or any application sensitive to edge artefacts, prefer Bilinear (no overshoot) over Bicubic or Lanczos — the higher mid-passband MTF of the latter comes with halos that downstream models may interpret as features.

For Sinc and (to a lesser extent) Lanczos, the underlying kernel produces ringing/halos near and beyond the $p_\text{eff}$-Nyquist frequency in the spatial domain. `opticks` returns $|H(\nu)|$ in the MTF chain (so the multiplied product remains well-defined), but the spatial-domain artefacts are real.
```

The shared sanity helper for processing-stage MTFs (`_check_processing_mtf`) therefore enforces only **finite** and **non-negative** values — there is no upper bound at 1, in contrast to detector-stage MTFs which never exceed unity. This is forward-compatible with future sharpening/boost filters that intentionally lift the MTF above unity.

### Summary

| Kernel             | Cost | Passband shape          | Overshoot? | Typical use                          |
| ------------------ | ---- | ----------------------- | ---------- | ------------------------------------ |
| Nearest-Neighbour  | Low  | Slow sinc roll-off      | No         | Quick-look / pre-visualisation       |
| Bilinear           | Low  | Aggressive sinc² fall   | No         | Default for ML / detection inputs    |
| Bicubic ($a=-0.5$) | Med  | Mid-band boost          | Yes        | Standard photographic resampling     |
| Lanczos-3          | Med  | Flattest passband       | Yes        | High-quality interpolation           |
| Sinc               | High | Brick-wall (rectangular)| Severe     | Theoretical reference only           |

The code-side walkthrough — `MTF_Model_1D.resampling`, the `Processing` component, and `ImagingChain` end-to-end — is in [`opticks` MTF](opticks_mtf.md#resampling-mtf).

---

[^1]: Electro-Optical Imaging System Performance, 6th Ed.; G. C. Holst; SPIE Press, 2017 — Ch. 10 (Image Reconstruction).

[^2]: Remote Sensing: Models and Methods for Image Processing, 3rd Ed.; R. A. Schowengerdt; Academic Press, 2007.

[^3]: Cubic Convolution Interpolation for Digital Image Processing; R. G. Keys; IEEE Transactions on Acoustics, Speech, and Signal Processing, Vol. 29, No. 6, pp. 1153–1160, 1981.
