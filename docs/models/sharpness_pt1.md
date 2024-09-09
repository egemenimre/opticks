# Modelling Contrast and Sharpness Performance

## Contrast and Sharpness for an Imaging System

An imager "sees" a scene and reproduces it in an image. How faithfully contrast and sharpness information is reproduced, or conversely how much it is degraded by the overall system is a significant performance characteristic.

For example, for a (theoretical) high contrast scene with a sharp "black to white" transition, the imager will likely generate an image that takes a few or several pixels to transition from pure black to a progressively lighter grey to pure white. A good imager will need less pixels for this transition, reproducing the sharp contrast more correctly.

As a minimum, both the optics and the detector will play a role in determining the sharpness performance of the imager (though there are also other concerns like aliasing). The theoretical limit of the resolution is given by the [spatial cut-off frequency](imager_geom.md#spatial-cut-off-frequency) for the optics and [Nyquist Frequency or Limit](imager_geom.md#nyquist-frequency-or-nyquist-limit) for the detector. In reality, the quality of the optics as well as defocus (due to design, manufacturing, deformations under mechanical loads or thermo-elastic distortions) will result in a degraded sharpness for the optics.

Furthermore, the "shaking" of the imager (for example the moving platform that the imager is mounted on, cryocoolers or other vibration sources) will degrade the sharpness. The imaging duration (or integration time) may cause motion blur if there is a relative motion between the scene and the imager. If the "end-to-end" sharpness needs to be estimated (for example for target identification or recognition tasks), atmospheric effects need to be considered as well. Therefore, all characteristics of the entire imaging system as well as the complete path of the light should be considered when modelling the overall image sharpness.

The concept of *frequency* in this context means how "quickly" consecutive pixels/regions change in contrast. For example, for a scene with black and white strips, "low frequency" means large black and white strips, whereas "high frequency" means thin black and white strips. Similar to this sharp "square wave", the frequency can also be used for a smooth and continuous "sine wave" change of the scene brightness in the spatial domain. Note that, sine wave frequencies are measured in "cycles per unit width on the image plane", usually "cycles/mm". Square wave frequencies are measured in "lines or line pairs per unit width on the image plane", usually "line pairs/mm", where each line pair consists one black and one white line.

The discussion above is related to the concept of resolution - the system may resolve lower frequencies, but as the frequency gets higher (and as the strips get thinner), the system will not be able to reproduce white and black - only light grey and dark grey, and for even higher frequencies, only a constant middle grey.

In practical terms, this means that a good optical system will resolve most of the detail it can theoretically resolve, whereas a bad one (for example with low quality optics or too much vibration) will fall well short of its theoretical limits and will generate a more "muddy" image.

## Modulation Transfer Function (MTF) and Modulation Contrast Function (MCF)

Modulation Transfer Function (MTF) is essentially how well the optical system converts the spatial modulation of the target object into the spatial modulation of the "image object". In practice, it is a measure of how well the input "resolution" or spatial detail information is transferred through an element of the imager system. In other words, it is a measure of how much an input data of a given frequency (or level of detail) is degraded.

Modulation Contrast Function (MCF) is similar, as it is essentially how well the optical system converts the target object contrast into image object contrast.

Both can be measured with sine wave or square wave targets. For a sine wave target (or sine wave response), MCF is equal to MTF. However, for a square wave target (or square wave response), MTF is close to, but not equal to MCF. For square wave response, the MCF is usually higher than MTF. From a Fourier analysis perspective, a square wave comprises an infinite number of frequencies and a sine wave is just a first order approximation with a single frequency.[^3]

This means that, for a sine wave input, the image contrast can be easily computed multiplying the MTF and the object contrast in the scene, or rather the Modulation Contrast of the Object (MCO). MCO depends on the object irradiance with respect to the background. Depending on the wavelength, the primary source of object irradiance can be reflectance or emittance.

While the equations will be given in the following sections, a couple of general remarks regarding the MTF of an optical system are due here. First and foremost, for an ideal optical system, a larger aperture will result in a better MTF. A smaller aperture will diffract the light more and will result in a lower image quality. This also means that, for an optical system with larger aperture, the PSF will be narrower. The second point is that, the optical MTF depends on the wavelength. As wavelength increases, the MTF for the same (ideal) optical system decreases and usually a larger aperture area is needed to compensate. For example, an ideal optical system for SWIR or LWIR will require a larger aperture (and overall a larger optical system) than an equivalent visible range optical system to get a similar level of resolution.

While MTF can be defined for just the optics and the detector of an imager in a narrow sense, other factors such as vibrations (usually called jitter) can be modelled as an MTF contributor. If considering the "end-to-end MTF performance", the drop in MTF due to atmosphere should be taken into account as well. In the end, depending on the definition or scope of the "system MTF", all MTF contributors are (as long as they are independent) simply multiplied to generate the system MTF, representing how well the input frequency be reproduced in the final image.

Even though MTF can be evaluated as a single value for a single input line frequency, it is more useful to evaluate it as the plot of possible line frequencies, starting from low frequencies (usually corresponding to high MTF values) to higher frequencies (usually with decreasing MTF values), all the way to the Nyquist limit of the detector, which sets the practical limit of the resolution. The following plot shows the Optics, Detector and the resulting Imager MTF decrease with the increasing input line frequency, as well as the Nyquist limit.

![Static MTF](images/static_mtf.png "Sample MTF plot")

A rule of thumb is that the resolution limit of the system is at the frequency for which the MTF is equal to 0.1[^1].

As can be imagined, the ratio of the optics limits to the resolution compared to the detector limits is a useful metric as to how good the output images will *theoretically* be. The equation is given simply as the ratio of the spatial cut-off frequency to the Nyquist limit, or for the given wavelength, the ratio of the F-number to the pixel pitch.

$$ Q = \frac {\lambda F_\#}{ \text{pix pitch}}  $$

It is desirable to have the Q value between 1 and 2. Below 1, the images are undersampled. Beyond 2, the image may become too blurry. The emphasis on *theoretically* should be noted, as the equation does not take into account optics defects like defocussing or surface imperfections. Nor does it take into account the real-world effects that cause blurring, such as vibrations.

It should also be noted that, MTF is different for every single point on the image plane (for example focus is usually less towards the edges). Furthermore, its value can be different for different directions, particularly when the lens is not axially symmetric. For example, optical imperfections may vary in tangential and sagittal directions or motion blur will result in lower MTF in the direction of the motion. The methodology discussed here is valid for a single point on the image plane and takes into account "cuts" of this 3D MTF surface in the vertical and horizontal directions.

---

The topic of sharpness and contrast performance is continued [here](sharpness_pt2.md).

[^1]: A Tutorial on Electro-Optical/Infrared (EO/IR) Theory and Systems; G. M. Koretsky, J. F. Nicoll, M. S. Taylor; Institute for Defense Analyses, IDA Document D-4642, 2013.

[^3]: The Art and Science of Optical Design; R. R. Shannon; Cambridge University Press; 1997.
