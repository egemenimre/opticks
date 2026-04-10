# Static Contributors to the MTF

The topic of sharpness and contrast performance is continued from [here](sharpness_pt1.md).

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

A $W_{RMS} = 1/14$ is likely a good starting point to model manufacturing errors [^11]. Some [sample fabrication tolerances are given here](https://www.telescope-optics.net/fabrication.htm). For example, surface roughness for Commercial Optics can be a single wavelength (Peak-to-Valley), whereas for Precision Optics it could be about quarter of a wavelength and for High Precision Optics it could be as low as 5% of a wavelength. Satellite imagers would also be as high as 5% of a wavelength.

## Detector MTF Models

### Detector Sampling MTF

Each detector pixel performs spatial averaging of the incident irradiance over its active area. More precisely, the scene irradiance is integrated with the detector's spatial responsivity function (which, for a uniform pixel, is a rectangle function of width equal to the pixel pitch). In the spatial domain this is a convolution, and in the frequency domain it becomes a multiplication with the Fourier transform of the rectangle — the sinc function[^3] (See also [here](https://spie.org/publications/spie-publication-resources/optipedia-free-optics-information/tt52_21_detector_footprint_mtf)):

$$\text{MTF}_\text{det sampling}(f) = \frac{\sin(\pi p f)}{\pi p f} = \text{sinc}(p f)$$

where $p$ is the pixel pitch and $f$ is the spatial frequency. Note the [normalised sinc function notation](https://en.wikipedia.org/wiki/Sinc_function) used.

This is the most fundamental detector MTF contributor — it applies universally to all detector types (CCD, CMOS, infrared) and is determined entirely by the pixel geometry. The key features of this MTF contributor are:

- At zero frequency ($f = 0$): MTF = 1 (the full signal is captured).
- At Nyquist ($f = f_{ny} = 1/2p$): $\text{MTF} = \text{sinc}(0.5) = 2/\pi \approx 0.637$. This means that even for a perfectly sampled signal at the Nyquist limit, the detector averages away about 36% of the modulation.
- At $f = 1/p$: MTF = 0. The spatial period of the input equals the pixel pitch, so the signal is entirely averaged out within a single pixel.
- Beyond $f = 1/p$: the sinc function goes negative, corresponding to contrast reversal — bright and dark regions swap. This is related to aliasing but is distinct from it; aliasing is a sampling effect while the sinc rolloff is an averaging (integration) effect.

For detectors with a fill factor less than unity (i.e., the photosensitive area is smaller than the pixel pitch), the effective averaging width is reduced and the MTF rolls off more slowly. In this case, $p$ in the formula should be replaced by the active detector width $d$ rather than the pixel pitch. Most modern sensors have fill factors close to 1 (aided by microlenses in CMOS sensors), so $p$ is commonly used directly.

When pixel binning is applied (combining $n \times n$ pixels), the effective pixel pitch increases to $n \times p$, and the sampling MTF degrades correspondingly — the sinc function becomes narrower in the frequency domain.

### Detector Diffusion MTF

When photons are absorbed in the detector substrate, they generate charge carriers (electron-hole pairs). These carriers do not stay exactly where they are generated — they diffuse laterally through the semiconductor material before being collected by the pixel's depletion region. This lateral spread means that a photon absorbed under one pixel can contribute signal to a neighbouring pixel, blurring the image and reducing the MTF.

The extent of this effect depends on the detector material and architecture:

- **CCD detectors**: Moderate diffusion effect. In front-illuminated CCDs, carriers generated deep in the substrate (from longer wavelength photons, particularly in the red and near-IR) have to diffuse further to reach the depletion region, leading to a wavelength-dependent MTF loss[^5]. Back-illuminated CCDs can have worse diffusion MTF unless thinned, as carriers generated near the back surface must traverse the full substrate thickness.

- **CMOS detectors**: Generally smaller diffusion effect than CCDs, and modern designs incorporate several architectural features that further reduce it[^5]:
  - **Thin epitaxial layers** (often <10 µm) limit the lateral diffusion distance.
  - **Deep Trench Isolation (DTI)**: Physical SiO2 barriers etched between pixels (typically 3–4 µm deep) that block lateral carrier diffusion at the pixel boundary. DTI effectively eliminates electrical crosstalk in well-isolated designs and is standard in modern small-pitch sensors.
  - **Pinned Photodiodes (PPD)**: The depletion structure extends deeper into the substrate and can be over-depleted, increasing carrier drift velocity and reducing the time (and distance) carriers diffuse before collection.
  - **Back-Side Illumination (BSI)**: The substrate is thinned (to ~3–5 µm) so that photons enter from the back and are absorbed close to the depletion region. This shortens the diffusion path but reverses the geometry compared to front-side illumination.

  The combination of these features means that diffusion MTF is often a secondary concern for visible-wavelength CMOS imagers. However, at longer wavelengths (NIR, >800 nm), photons penetrate deeper and the diffusion effect increases — CMOS MTF at 850 nm can be significantly lower than at visible wavelengths.

- **Infrared detectors (HgCdTe, InSb)**: Diffusion is a **dominant** MTF contributor. These detectors typically have thick absorber layers (to achieve sufficient quantum efficiency at long wavelengths) and relatively large diffusion lengths. In HgCdTe detectors, the minority carrier diffusion length can be tens of micrometres — comparable to or larger than the pixel pitch — making lateral diffusion the primary MTF-limiting mechanism in many IR focal plane arrays[^6]. This gets more important as the pixel pitch gets smaller. InSb detectors (covering 1–5.5 µm) face similar challenges. The effect worsens at longer wavelengths where thicker absorbers are required.

#### Diffusion MTF Model

##### General Model — Crowell & Labuda (1969)

The most general analytical diffusion MTF model was derived by Crowell & Labuda[^10] for the silicon diode array camera tube. Their Eq. (5) gives the quantum efficiency as a function of spatial frequency $k$ (here written as $f$ in cycles/length), accounting for surface recombination velocity $S$, back-surface reflectivity $R$, and a finite substrate.

In the C&L geometry, light enters the **field-free region** first (from the illuminated surface at $x = 0$), passes through the field-free region of thickness $L_a$, and is collected at the **depletion region** boundary at $x = L_a$. The substrate back surface is at $x = L_b$, so the depletion region width is $L_b - L_a$.

The frequency-dependent quantum efficiency $\eta(f)$ is:

$$\eta(f) = \frac{\alpha_{abs}}{\alpha_{abs}^2 - 1/L_f^2} \left[ \frac{\beta^{-}\, e^{-\alpha_{abs}\, L_a} - \beta^{+}\, e^{\alpha_{abs}\, L_a}}{\beta^{-} - \beta^{+}\, e^{2 L_a / L_f}} + \alpha_{abs}\, L_f \right] + \frac{e^{-\alpha_{abs}\, L_a} - e^{-\alpha_{abs}\, L_b}}{1 + \alpha_{abs}\, L_f} - \frac{(1 - R)\, e^{-\alpha_{abs}\, L_b}}{1 + \alpha_{abs}\, L_f}$$

where:

$$\beta^{\pm} = \left(1 \pm \frac{S\, L_f}{D}\right) e^{\pm L_a / L_f}$$

and the frequency-dependent auxiliary length is:

$$L_f = \frac{L_{Diff}}{\sqrt{1 + (2\pi\, L_{Diff}\, f)^2}}$$

The diffusion MTF is the normalised form:

$$\text{MTF}_\text{diff}(f) = \frac{\eta(f)}{\eta(0)}$$

At $f = 0$, $L_f = L_{Diff}$, so $\eta(0)$ is obtained by substituting $L_f = L_{Diff}$ in the expression above.

The parameters are:

- $f$ : spatial frequency
- $\alpha_{abs}(\lambda)$ : photon absorption coefficient of the substrate material (wavelength-dependent)
- $L_a$ : field-free region thickness (distance from illuminated surface to the depletion region boundary)
- $L_b$ : total substrate thickness (distance from illuminated surface to the back surface); depletion width is $L_b - L_a$
- $L_{Diff}$ : minority carrier diffusion length ($L_{Diff} = \sqrt{D\tau}$, where $D$ is the diffusion coefficient and $\tau$ is the carrier lifetime)
- $S$ : surface recombination velocity at the illuminated surface (cm/s)
- $D$ : minority carrier diffusion coefficient (cm²/s); note $S/D$ has units of 1/length
- $R$ : reflectivity of the back surface (dimensionless, $0 \leq R \leq 1$)

The $\beta^{\pm}$ terms encode the boundary condition at the illuminated surface: $S$ controls how quickly carriers recombine there. The $(1-R)$ term accounts for carriers that reach the back surface and are either reflected back into the substrate (fraction $R$) or lost (fraction $1-R$).

##### Simplified Models

The general C&L equation is simplified for five practically important cases by applying boundary condition limits and geometry choices. Full derivations are in the [derivation notebook](misc/detector_diffusion_derivation.ipynb); numerical examples with real detector parameters are in the [numerical examples notebook](misc/detector_diffusion_mtf.ipynb). In all cases $R = 0$ (no back-surface reflection).

| Model | Geometry | Detector category | Conditions | Physical meaning | Example detector |
| --- | --- | --- | --- | --- | --- |
| **BSI-1** | BSI | Scientific CCD | $S \to \infty$, finite $L_b$ | Dead back surface; all carriers reaching the illuminated surface recombine | Teledyne e2v CCD47-10 |
| **BSI-2** | BSI | Scientific sCMOS | $S = 0$, finite $L_b$ | Perfectly passivated back surface; carriers reflect back toward junction | Gpixel GSENSE400BSI |
| **BSI-3** | BSI | IR FPA | $L_b \to \infty$, general $S$ | Thick substrate; finite surface recombination | Teledyne H2RG (MCT) |
| **FSI-1** | FSI | Consumer CMOS | Dead far surface, finite $L_a$ | Light enters depletion first; finite field-free bulk | Sony IMX174 |
| **FSI-2** | FSI | Bulk CCD | Dead far surface, $L_a \to \infty$ | Semi-infinite bulk; equivalent to the Seib[^4] / Fiete[^7] model | Generic FSI CCD |

**BSI geometry** ($L_a$ = field-free depth, $L_b$ = depletion depth): light enters the field-free region first (from the back surface), diffuses to the depletion edge at $x = L_a$, and the depletion layer extends from $L_a$ to $L_b$.

**FSI geometry** (reversed): light enters the depletion region first ($L_D$ = depletion width), and the field-free bulk of thickness $L_a$ is behind it. The BSI-1 through BSI-3 equations follow directly from C&L Eq. (5); the FSI cases require an independent derivation from the diffusion equation with reversed boundary conditions.

The strong wavelength dependence of $\alpha_{abs}(\lambda)$ is physically significant in all cases: longer wavelengths have smaller $\alpha_{abs}$, so photons are absorbed deeper, carriers travel further before collection, and the diffusion MTF degrades accordingly.

##### Applicability

| Detector type | Model | $L_a$ | $L_b$ / $L_D$ | Example |
| --- | --- | --- | --- | --- |
| Sci-CCD (BSI, thinned) | BSI-1 | 15–20 µm | 10–20 µm | e2v CCD47-10 |
| sCMOS (BSI, thin epi) | BSI-2 | 3–8 µm | 1–2 µm | Gpixel GSENSE400BSI |
| MCT IR array (BSI) | BSI-3 | 5–15 µm | $L_b \to \infty$ | Teledyne H2RG |
| Consumer CMOS (FSI) | FSI-1 | 3–10 µm | 1–3 µm | Sony IMX174 |
| Bulk CCD (FSI) | FSI-2 | $L_a \to \infty$ | 1–5 µm | Generic FSI CCD |

The model assumes an uninterrupted substrate with no lateral barriers to diffusion. For modern CMOS sensors with Deep Trench Isolation (DTI), the physical trenches between pixels block lateral carrier movement, and the actual diffusion MTF will be better than the model predicts. In such cases the model provides a conservative (pessimistic) estimate. For CCDs without pixel-level isolation barriers, the model is directly applicable.

Note that the model is isotropic — the diffusion MTF is the same in the ALT and ACT directions for a given spatial frequency.

#### Typical Diffusion Length Values

| Detector Type | Typical $L_{Diff}$ | Notes |
| --- | --- | --- |
| Front-illuminated CCD (visible) | 1–5 µm | Wavelength dependent; worse in red/NIR |
| BSI CMOS (visible) | <2 µm | Thin substrates limit diffusion |
| HgCdTe (MWIR, 3–5 µm) | 10–30 µm | Can exceed pixel pitch |
| HgCdTe (LWIR, 8–12 µm) | 15–50 µm | Thick absorbers, long diffusion lengths |
| InSb (MWIR) | 10–40 µm | Similar magnitude to HgCdTe MWIR |

### Detector Crosstalk MTF

Crosstalk occurs when signal generated by a photon in one pixel leaks into adjacent pixels. There are two distinct mechanisms[^5][^7]:

- **Electrical crosstalk**: Photogenerated charge carriers diffuse laterally in the substrate and are collected by a neighbouring pixel's depletion region. This is closely related to the diffusion effect discussed above, but treated here as an inter-pixel coupling rather than a continuous spread. It is the dominant crosstalk mechanism in most detectors.

- **Optical crosstalk**: Photons scatter, reflect, or diffract within the detector structure (e.g., off metal interconnect layers, microlens edges, or colour filter boundaries) and are absorbed in a neighbouring pixel. This is more significant in small-pitch sensors with tall pixel structures.

The net effect is that a point source illuminating a single pixel produces a response in its neighbours, reducing modulation at high spatial frequencies where the pattern alternates between pixels.

#### Nearest-Neighbour Crosstalk Model

The crosstalk model uses the **Discrete Impulse Response** method[^7]. A unit signal (1) is generated in the target pixel, and a fraction of that charge leaks to its available neighbours. We define separate coefficients for the two types of neighbour:

- $X_s$: crosstalk to a **side** (edge-adjacent) neighbour
- $X_d$: crosstalk to a **diagonal** (corner-adjacent) neighbour

For a centre pixel with 8 neighbours (4 side, 4 diagonal), the discrete kernel is:

$$
\begin{bmatrix}
X_d & X_s & X_d \\
X_s & (1 - 4X_s - 4X_d) & X_s \\
X_d & X_s & X_d
\end{bmatrix}
$$

The kernel sums to unity (signal is conserved). Its 2D Discrete Fourier Transform gives the transfer function:

$$H(f_x, f_y) = (1 - 4X_s - 4X_d) + 2X_s \cos(2\pi f_x p) + 2X_s \cos(2\pi f_y p) + 4X_d \cos(2\pi f_x p)\cos(2\pi f_y p)$$

where $p$ is the pixel pitch, $f_x$ and $f_y$ are the spatial frequencies in the two detector axes, and $X_s$ and $X_d$ are dimensionless fractions. Because the kernel is symmetric, $H$ is real-valued and equals the MTF directly.

Taking a 1D slice along one axis ($f_y = 0$) yields:

$$\text{MTF}_\text{xtalk}(f) = 1 - 2(X_s + 2X_d)\left(1 - \cos(2\pi f p)\right)$$

When diagonal crosstalk is negligible ($X_d = 0$), this reduces to the classical formula $\text{MTF}_\text{xtalk}(f) = 1 - 2X_s(1 - \cos(2\pi f p))$.

This model has the property that:

- At zero frequency ($f = 0$): MTF = 1 (signal is conserved, just redistributed)
- At Nyquist ($f = f_{ny} = 1/2p$): $\text{MTF} = 1 - 4(X_s + 2X_d)$ (maximum degradation, where the pattern alternates every pixel)

The constraint $4X_s + 4X_d < 1$ is required for the kernel centre weight to remain positive. In practice, typical values for a high-performance detector are $X_s$ = 1%–4% and $X_d$ = 0.2%–1%.

#### Crosstalk Across Detector Types

| Detector Type | Typical $X_s$ | Primary Mechanism | Notes |
| --- | --- | --- | --- |
| CCD (visible) | 0.01–0.03 | Electrical | Well-controlled in modern designs |
| BSI CMOS (visible) | 0.01–0.05 | Optical + Electrical | Optical crosstalk increases at small pitch (<2 µm)[^5] |
| HgCdTe (IR) | 0.02–0.10 | Electrical | Lateral diffusion in thick absorber; **can be dominant** MTF limiter[^6] |
| InSb (IR) | 0.02–0.08 | Electrical | Similar to HgCdTe; depends on pixel pitch vs. diffusion length |

For HgCdTe detectors, crosstalk can reduce the total system MTF by over 30% at the 50% MTF frequency[^6]. Ion implantation guard-ring structures around each pixel can suppress electrical crosstalk, and this is a common design mitigation in modern IR focal plane arrays.

#### Relationship Between Diffusion and Crosstalk

The diffusion MTF and the electrical component of crosstalk MTF describe the **same underlying physics** — lateral movement of charge carriers in the substrate — but model it at different levels of abstraction[^5][^6]:

- **Diffusion MTF** is a continuous, physics-based model. It uses the material diffusion length $L_d$ as its parameter and does not reference pixel boundaries. It models the full spatial spread of the charge cloud, including long-range tails that can extend over multiple pixels.

- **Crosstalk MTF** is a discrete, measurement-based model. It uses the inter-pixel coupling coefficients $X_s$ and $X_d$ as its parameters, typically obtained from detector characterisation data. It captures only the nearest-neighbour signal sharing and does not model the shape of the underlying spread.

Because of this shared physical origin, **the two should not be blindly multiplied** in an MTF budget when the crosstalk is predominantly electrical. Doing so would double-count the same effect. The recommended approach is:

- If material properties are known (diffusion length), use **diffusion MTF alone** for the electrical component.
- If only measured inter-pixel coupling data is available, use **crosstalk MTF alone**.
- If both optical and electrical crosstalk are significant (e.g., small-pitch BSI CMOS), the diffusion MTF can be used for the electrical component and the crosstalk MTF for the **optical component only**, with $X_s$ representing only the optical coupling fraction.

### CTE (Charge Transfer Efficiency) MTF

This section applies **only to CCD detectors** and is not relevant for CMOS or infrared focal plane arrays, which read out pixels individually rather than shifting charge through a register.

In a CCD, the accumulated charge in each pixel is read out by sequentially transferring it through neighbouring pixels towards the output amplifier. Each transfer is not perfectly efficient — a small fraction of the charge is left behind in traps within the silicon (and released later, into the wrong pixel).

This *charge transfer inefficiency* (CTI = $1 - \varepsilon$, where $\varepsilon$ is the CTE) causes a directional smearing of the image along the transfer direction, which is similar to motion blur in the transfer direction. The effect is cumulative: pixels further from the readout register undergo more transfers and suffer more degradation. For a CCD with $n$ columns, the last pixel to be read out experiences $n$ transfers (or $n \times n_\text{phases}$ for a multi-phase CCD). The standard model[^8] for the resulting MTF is:

$$\text{MTF}_\text{CTE}(f) = \exp\left(-N(1-\varepsilon)\left(1 - \cos\left(2 \pi fp\right)\right)\right)$$

where:

- $N$ is the number of charge transfers (depends on pixel position and CCD architecture)
- $\varepsilon$ is the CTE per transfer (dimensionless)
- $p$ is pixel pitch
- $f$ is the spatial frequency

It is also possible to convert the equation to Nyquist frequency, where $f_{ny}= 1 / 2p$.

Note that, Boreman [^5] uses a reversed terminology, where $\varepsilon$ is the fractional charge left behind at each transfer, or CTI as defined here.

At zero frequency ($f = 0$), the cosine term equals 1 and the MTF is unity — charge is conserved overall. The worst degradation occurs at the Nyquist frequency, where the cosine term equals $-1$, giving $\text{MTF} = \exp(-2N(1-\varepsilon))$.

#### CTE Values and Impact

Modern scientific CCDs achieve very high CTE values[^8]:

| Application | Typical CTE | CTI ($1-\varepsilon$) |
| --- | --- | --- |
| Consumer / industrial CCDs | 0.99999 | $10^{-5}$ |
| Scientific CCDs (new) | 0.999999 | $10^{-6}$ |
| Space-qualified CCDs (new) | 0.9999995 | $5 \times 10^{-7}$ |
| Radiation-damaged CCDs (space, after years) | 0.99999–0.9999 | $10^{-5}$–$10^{-4}$ |

Three-phase CCDs have significantly higher CTE than two-phase CCDs.

While CTI values appear tiny, the cumulative effect over thousands of transfers is significant. For example, a 4096-pixel CCD with 3-phase architecture and $\varepsilon = 0.99999$ undergoes $N = 4096 \times 3 = 12288$ transfers for the furthest pixel. At the Nyquist frequency: $\text{MTF} = \exp(-2 \times 12288 \times 10^{-5}) = \exp(-0.246) \approx 0.78$. Also note that the CTI values are different in horizontal and vertical directions.

Radiation damage in space progressively degrades CTE by creating charge traps in the silicon lattice[^9]. This is a major concern for long-duration space missions (e.g., the [Hubble Space Telescope ACS/WFC](https://hst-docs.stsci.edu/wfc3ihb/chapter-5-wfc3-detector-characteristics-and-performance/5-4-wfc3-ccd-characteristics-and-performance/#id-5.4WFC3CCDCharacteristicsandPerformance-cte5.4.11ChargeTransferEfficiency) suffered from degradation and required dedicated CTE correction algorithms after years in orbit).

CMOS detectors do not suffer from MTF degradation due to CTE because each pixel has its own amplifier and is read out individually — there is no charge transfer between pixels.

### Directionality: ALT vs ACT for Pushbroom Detectors

For pushbroom (line scanner) detectors, the along-track (ALT) and across-track (ACT) directions have fundamentally different characteristics, and the detector MTF contributors are **not necessarily the same in both directions**:

- **Sampling MTF**: In the ACT direction, the sampling MTF is determined by the physical pixel pitch as usual. In the ALT direction, the effective "pixel size" is determined by the integration time multiplied by the image velocity on the focal plane — this is often different from the physical pixel pitch. The sinc MTF must therefore be computed separately for each direction, using the appropriate effective pitch.

- **Diffusion MTF**: The lateral diffusion of charge carriers is an **isotropic** process in the semiconductor — carriers spread equally in all directions. The diffusion MTF is therefore the same in both ALT and ACT for a given spatial frequency. However, because the Nyquist frequencies may differ between ALT and ACT (due to different effective pixel sizes), the relative impact of diffusion at each direction's Nyquist will differ.

- **Crosstalk MTF**: Like diffusion, the nearest-neighbour coupling is generally **isotropic** for square pixels. However, for rectangular pixels or detectors with asymmetric structures (e.g., transfer gates or bus lines running in one direction), the crosstalk coefficients $X_s$ and $X_d$ can differ between ALT and ACT.

- **CTE MTF**: This is **inherently directional**. Charge is transferred along one specific axis — typically along the column (ALT direction for a pushbroom) for the parallel register, and along the row (ACT direction) for the serial register. Each direction has a different number of transfers and potentially different CTE values. For a TDI (Time Delay Integration) CCD, the parallel CTE applies in the ALT direction with $N$ equal to the number of TDI stages times the number of phases.

As a consequence, for pushbroom and TDI systems, the full detector MTF budget should be computed separately for the ALT and ACT directions.

### Summary: Detector MTF Sources by Detector Type

The following table summarises the applicability and typical impact of each detector-level MTF contributor:

| MTF Source | CCD (Visible) | CMOS (Visible) | HgCdTe / InSb (IR) |
| --- | --- | --- | --- |
| **Sampling** (sinc) | Always applies | Always applies | Always applies |
| **Diffusion** | Moderate (1–5 µm $L_d$) | Small (<2 µm $L_d$, thin substrates) | **Large / dominant** (10–50 µm $L_d$) |
| **Crosstalk** | Small (well-controlled) | Small–Moderate (optical at small pitch) | **Large** (electrical, thick absorber) |
| **CTE** | Applies (degrades with radiation) | **Not applicable** | **Not applicable** |

For visible-wavelength systems, the detector sampling MTF typically dominates, with diffusion and crosstalk as secondary contributors. For infrared systems, diffusion and crosstalk can be the **primary** MTF-limiting mechanisms, often exceeding the degradation from optical aberrations.

---

The topic of sharpness and contrast performance is continued [here](sharpness_pt3.md).

[^1]: The Infrared & Electro-Optical Systems Handbook; J. S. Accetta, David L. Shumaker (Ed.); Infrared Information Analysis Center, 1993.

[^2]: The Art and Science of Optical Design; R. R. Shannon; Cambridge University Press; 1997.

[^3]: A System Engineering Approach to Imaging; N. S. Kopeika; SPIE Press, 1998.

[^4]: Carrier diffusion degradation of modulation transfer function in charge coupled imagers; D. H. Seib; IEEE Transactions on Electron Devices, Vol. 21, No. 3, 1974.

[^5]: Modulation Transfer Function in Optical and Electro-Optical Systems, 2nd Ed.; G. D. Boreman; SPIE Tutorial Texts, 2001.

[^6]: Electro-Optical Imaging System Performance, 6th Ed.; G. C. Holst; SPIE Press, 2017.

[^7]: Modeling the Imaging Chain of Digital Cameras; R. D. Fiete; SPIE Tutorial Texts, 2010.

[^8]: Scientific Charge-Coupled Devices; J. R. Janesick; SPIE Press, 2001.

[^9]: Charge Transfer Efficiency in Proton Damaged CCDs; T. Hardy, R. Murowinski, M. J. Deen; IEEE Transactions on Nuclear Science, Vol. 45, No. 2, 1998.

[^10]: The Silicon Diode Array Camera Tube; M. H. Crowell, E. F. Labuda; The Bell System Technical Journal, Vol. 48, No. 5, pp. 1481–1528, 1969.

[^11]: CMOS/CCD sensors and camera systems; G. C. Holst, T.S. Lomheim; JCD publishing, SPIE, 2nd ed., 2011.
