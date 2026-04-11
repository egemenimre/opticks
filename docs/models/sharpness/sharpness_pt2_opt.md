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

A $W_{RMS} = 1/14$ is likely a good starting point to model manufacturing errors [^11]. Some [sample fabrication tolerances are given here](https://www.telescope-optics.net/fabrication.htm). For example, surface roughness for Commercial Optics can be a single wavelength (Peak-to-Valley), whereas for Precision Optics it could be about quarter of a wavelength and for High Precision Optics it could be as low as 5% of a wavelength. Satellite imagers would also be as high as 5% of a wavelength.

---

The topic of sharpness and contrast performance is [continued here](sharpness_pt2_det.md).

[^1]: The Infrared & Electro-Optical Systems Handbook; J. S. Accetta, David L. Shumaker (Ed.); Infrared Information Analysis Center, 1993.

[^2]: The Art and Science of Optical Design; R. R. Shannon; Cambridge University Press; 1997.

[^11]: CMOS/CCD sensors and camera systems; G. C. Holst, T.S. Lomheim; JCD publishing, SPIE, 2nd ed., 2011.
