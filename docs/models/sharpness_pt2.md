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

$$\text{MTF}_\text{aberr opt}(f) = \text{MTF}_\text{ideal opt}(f) \times \text{ATF}(f) $$  

Multiplying the ATF value with the ideal optical MTF, we can reach a more realistic MTF with the aberrations. As the $W_{RMS}$ value increases, the ATF value decreases and the resulting MTF also decreases, corresponding to a degradation in image quality.

Some sample fabrication tolerances are given [here](https://www.telescope-optics.net/fabrication.htm). For example, surface roughness for Commercial Optics can be a single wavelength (Peak-to-Valley), whereas for Precision Optics it could be about quarter of a wavelength and for High Precision Optics it could be as low as 5% of a wavelength. Satellite imagers would also be as high as 5% of a wavelength.

## Detector Sampling MTF

Each detector (or a single pixel of the detector) performs spatial averaging of the irradiance, or more precisely, we integrate the irradiance with the detector responsivity, over the detector area (See [here](https://spie.org/publications/spie-publication-resources/optipedia-free-optics-information/tt52_21_detector_footprint_mtf) for more information). In the frequency domain this corresponds to:

$$\text{MTF}_\text{detector}(f) = \frac{\sin(\pi p f)}{\pi p f} = \text{sinc}(p f)$$

where $p$ is equal to pixel pitch. Note the [normalised sinc function notation](https://en.wikipedia.org/wiki/Sinc_function) used.

As the frequency increases, and the wave periods become comparable to the detector pixel size, the pixels cannot represent the sine wave properly and there is a reduction in modulation. At a line frequency corresponding to the inverse of the pixel pitch, modulation goes down to zero, as the input sine wave is completely inside the pixel pitch. Even higher input frequencies will then be completely undersampled. This results in contrast reversal and MTF values will be negative.

## Diffusion MTF

This is similar to a crosstalk in the detector due to electronics, see chap 13 pg 404 of a Sys eng approach [^3].

In the absence of any better information, it can be taken as a 90% constant reduction to the MTF.

TBW

---

The topic of sharpness and contrast performance is continued [here](sharpness_pt3.md).

[^1]: The Infrared & Electro-Optical Systems Handbook; J. S. Accetta, David L. Shumaker (Ed.); Infrared Information Analysis Center, 1993.

[^2]: The Art and Science of Optical Design; R. R. Shannon; Cambridge University Press; 1997.

[^3]: A System Engineering Approach to Imaging; N. S. Kopeika; SPIE Press, 1998.
