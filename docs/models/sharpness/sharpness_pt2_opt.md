# Static Contributors to the MTF: Optics

The topic of sharpness and contrast performance is [continued from here](sharpness_pt1.md).

## Optical MTF

### Ideal Optical MTF

The *ideal* optical MTF for a clear circular diffraction-limited aperture with monochromatic illumination is given as[^1]:

$$\text{MTF}_\text{ideal opt}(f) = \frac{2}{\pi} \left[ \arccos \left( \frac{f}{f_c} \right) - \frac{f}{f_c}  \sqrt{1- \left( \frac{f}{f_c} \right)^2} \right]$$

where $f$ is the input line frequency and $f_c$ is the [spatial cut-off frequency](imager_geom.md#spatial-cut-off-frequency), which is a inversely proportional to the wavelength.

This can also be written as:

$$\text{MTF}_\text{ideal opt}(f) = \frac{2}{\pi} \left[ \arccos(\nu) - \nu \sqrt{1- (\nu)^1} \right]$$

where $\nu =\left( \frac{f}{f_c} \right) $.

An alternative formulation is:

$$\text{MTF}_\text{ideal opt}(f) = \frac{2}{\pi} \left[ \psi - \cos(\psi) \sin(\psi) \right]$$

where $\psi$ is equal to $\arccos(\nu)$.

The MTF value (or the curve) is wavelength dependent. The optical MTF is usually expressed either as multiple curves for each wavelength within the limits of the imager (for example red, green and blue for a colour imager) or as a single curve with the weighted average of multiple curves.

If there is an obstruction in front of the aperture (such as the secondary mirrors in most reflecting telescopes), the MTF curve is affected. It lowers the MTF in lower frequencies and increases it in the higher frequencies.

The real optical MTF will be lower than this value, due to real world design limitations, materials, manufacturing, integration as well as mechanical and thermal loads. This is usually simulated in a software like Zemax and eventually measured in the lab.

### Aberration Transfer Factor (ATF) and Aberrated Optical MTF

The Aberration Transfer Factor[^2] (ATF) is an empirical model that combines all sources of optical aberrations into a single total wavefront error:

$$\text{ATF}(f) = 1- \left( \frac{W_{RMS}}{0.18} \right)^2 \left[ 1 - 4 ( \nu -0.5 )^2 \right] $$

where $W_{RMS}$ is the RMS of the total wavefront error, or how much the actual wavefront deviates from the ideal wavefront. The unit of this deviation is the multiple wavelengths. For $W_{RMS} = 0$, ATF is equal to 1 for all input frequencies, corresponding to no optical aberrations.

All aberration sources are (RMS) summed and then the resulting wavelength error can be inserted in the ATF to compute the total ATF of the optical system.

$$\text{MTF}_\text{aberr opt}(f) = \text{MTF}_\text{ideal opt}(f) \times \text{ATF}(f)$$

Multiplying the ATF value with the ideal optical MTF, we can reach a more realistic MTF with the aberrations. As the $W_{RMS}$ value increases, the ATF value decreases and the resulting MTF also decreases, corresponding to a degradation in image quality.

A $W_{RMS} = 1/14$ is likely a good starting point to model manufacturing errors [^3]. Some [sample fabrication tolerances are given here](https://www.telescope-optics.net/fabrication.htm). For example, surface roughness for Commercial Optics can be a single wavelength (Peak-to-Valley), whereas for Precision Optics it could be about quarter of a wavelength and for High Precision Optics it could be as low as 5% of a wavelength. Satellite imagers would also be as high as 5% of a wavelength.

### Field Dependence of the Optical MTF

#### Why a single PSF is not enough

The ideal and aberrated optical MTF models above describe the contrast transfer for **one point on the image plane**. When we build a pupil function and propagate it to the focus via the DFT, the result is the PSF for whatever wavefront (OPD) was loaded into the pupil — there is no notion of *where* in the field of view that PSF applies.

In reality, the wavefront error is a function of the field position. Coma, astigmatism and field curvature all grow with the distance from the optical axis, so an off-axis pixel sees a different — and generally worse — PSF and MTF than the on-axis pixel.

A pure output-grid shift (such as `prysm`'s `focus_fixed_sampling(shift=...)`) does **not** produce a field-dependent PSF. It only re-centres the sampling window over an already-computed diffraction pattern. A genuinely field-dependent PSF must come from a **field-dependent OPD**.

#### Third-Order Wave-Aberration Expansion (Seidel / Hopkins)

The classical third-order (Seidel) wave-aberration function in the Hopkins notation is[^5]:

$$W(H, \rho, \phi) = W_{040}\,\rho^{4} + W_{131}\,H\,\rho^{3}\cos\phi + W_{222}\,H^{2}\,\rho^{2}\cos^{2}\phi + W_{220}\,H^{2}\,\rho^{2} + W_{311}\,H^{3}\,\rho\cos\phi$$

where $H$ is the normalised radial field coordinate ($0$ at centre, $1$ at the edge of the field), $\rho$ is the normalised pupil radius, and $\phi$ is the azimuthal angle in the pupil. The five terms are:

| Term | Name | Field dependence |
| --- | --- | --- |
| $W_{040}\,\rho^{4}$ | Spherical aberration | None (field-independent) |
| $W_{131}\,H\,\rho^{3}\cos\phi$ | Coma | $\propto H$ |
| $W_{222}\,H^{2}\,\rho^{2}\cos^{2}\phi$ | Astigmatism | $\propto H^{2}$ |
| $W_{220}\,H^{2}\,\rho^{2}$ | Field curvature (Petzval) | $\propto H^{2}$ |
| $W_{311}\,H^{3}\,\rho\cos\phi$ | Distortion | $\propto H^{3}$, does not affect MTF |

Each coefficient (e.g. $W_{131}$) represents the *peak* wavefront deviation at the edge of the field ($H = 1$). The total RMS wavefront error at field position $H$ is[^4]:

$$W_{\mathrm{rms}}(H) = \sqrt{\sum_i \sigma_i^{2}(H)}$$

where each per-term RMS contribution $\sigma_i$ is obtained via the standard Mahajan variance coefficients[^4] (which account for the balanced Zernike representation of each Seidel term).

#### Plugging into the ATF

Once $W_{\mathrm{rms}}(H)$ is known, the ATF model from the preceding section becomes field-aware:

$$\text{MTF}_\text{aberr opt}(f;\, H) = \text{MTF}_\text{ideal opt}(f) \times \text{ATF}\!\left(f,\; W_{\mathrm{rms}}(H)\right)$$

This does not introduce a new MTF formula — it simply makes the existing ATF framework sensitive to the field position by supplying a field-dependent $W_{\mathrm{rms}}$.

#### Modelling Tiers

Three practical paths are available, in decreasing order of fidelity:

**Tier A — Measurements or ray-trace data (highest fidelity).** When lab MTF measurements or Zemax-exported data exist for multiple field points, they can be ingested directly:

- **Per-field MTF curves** via {py:meth}`.MTF_Model_1D.external_data` — the user supplies spatial-frequency and MTF arrays for each field point.
- **Per-field Zernike coefficients** from Zemax or CODE V — each coefficient set is converted to an OPD via {py:meth}`.OptPathDiff.from_zernike`, attached to the optics via {py:meth}`.Optics.add_mono_pupil_function`, and the PSF/MTF computed via {py:meth}`.Optics.compute_psf` / {py:meth}`.Optics.mtf`. `OptPathDiff.from_zernike` accepts three index orderings via the `ordering` argument — **ANSI/OSA** (default, 0-based), **Noll** (1-based), and **Fringe** (1-based). Both Zemax OpticStudio and CODE V default to Fringe ordering, so coefficients exported from either tool can be ingested directly by setting `ordering="fringe"`.

No new API is needed for Tier A — all required methods already exist. The user manages a collection of per-field models keyed by the field position.

**Tier B — Analytical Seidel model + prysm PSF.** {py:meth}`.FieldAberrationModel.to_zernikes` returns an ANSI-ordered Zernike coefficient vector at any field point $(H_x, H_y)$. This vector drives the full `prysm` PSF pipeline through {py:meth}`.OptPathDiff.from_zernike` and {py:meth}`.Optics.add_mono_pupil_function`, yielding an actual field-dependent PSF image and 2D MTF at the cost of one DFT propagation per field point. This is lower fidelity than real measurements but accounts for diffraction and aberration coupling properly.

**Tier C — Analytical Seidel model + empirical ATF (fastest).** The {py:class}`.FieldAberrationModel` class stores the four Seidel coefficients at the edge of the field and computes $W_{\mathrm{rms}}(H)$ at any field position. The result feeds directly into {py:meth}`.MTF_Model_1D.emp_model_aberrated_optics` (via {py:meth}`.Optics.field_mtf_model_1d`) to produce a 1D MTF curve — no `prysm` propagation needed. This is the fastest path for MTF-budget work during early system sizing, but the empirical ATF does not model aberration-specific PSF shapes (e.g. comatic tails).

#### Typical Seidel Coefficient Values

The following table gives rough $W_{\mathrm{rms}}$ ranges at the edge of the field ($H = 1$) for different classes of optics, as a guideline for early sizing[^2][^4][^6]:

| Optics class | $W_{\mathrm{rms}}$ at edge | Notes |
| --- | --- | --- |
| Diffraction-limited space telescope | $\leq \lambda/14$ | Hubble-class, well-corrected |
| High-quality satellite imager | $\lambda/14$ to $\lambda/10$ | Multi-element corrected design |
| Precision commercial lens | $\lambda/10$ to $\lambda/4$ | Machine-vision, microscopy |
| Consumer camera lens (wide-angle) | $\lambda/4$ to $\lambda/2$ | Significant field curvature |
| Uncorrected singlet | $\geq \lambda$ | Dominated by spherical + coma |

---

The topic of sharpness and contrast performance is [continued here](sharpness_pt2_det.md).

[^1]: The Infrared & Electro-Optical Systems Handbook; J. S. Accetta, David L. Shumaker (Ed.); Infrared Information Analysis Center, 1993.

[^2]: The Art and Science of Optical Design; R. R. Shannon; Cambridge University Press; 1997.

[^3]: CMOS/CCD sensors and camera systems; G. C. Holst, T.S. Lomheim; JCD publishing, SPIE, 2nd ed., 2011.

[^4]: Aberration Theory Made Simple, 2nd Ed.; V. N. Mahajan; SPIE Press, 2011.

[^5]: Wave Theory of Aberrations; H. H. Hopkins; Oxford University Press, 1950.

[^6]: Basic Wave-front Aberration Theory for Optical Metrology; J. C. Wyant, K. Creath; Chapter 1 of Applied Optics and Optical Engineering, Vol. XI, Academic Press, 1992.
