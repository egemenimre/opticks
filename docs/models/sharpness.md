# Modelling Contrast and Sharpness Performance

## Contrast and Sharpness for an Imaging System

An imager "sees" a scene and reproduces it in an image. How faithfully contrast and sharpness information is reproduced, or conversely how much it is degraded by the overall system is a significant performance characteristic.

For example, for a (theoretical) high contrast scene with a sharp "black to white" transition, the imager will likely generate an image that takes a few or several pixels to transition from pure black to a progressively lighter grey to pure white. A good imager will need less pixels for this transition, reproducing the sharp contrast more correctly.

As a minimum, both the optics and the detector will play a role in determining the sharpness performance of the imager (though there are also other concerns like aliasing). The theoretical limit of the resolution is given by the [spatial cut-off frequency](imager_geom#spatial-cutoff-frequency) for the optics and [Nyquist Frequency / Limit](imager_geom#nyquist-frequency-Limit) for the detector. In reality, the quality of the optics as well as defocus (due to design, manufacturing, deformations under mechanical loads or thermo-elastic distortions) will result in a degraded sharpness for the optics.

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

## Static Contributors to the MTF

### Optical MTF

#### Ideal Optical MTF

The *ideal* optical MTF for a clear circular diffraction-limited aperture with monochromatic illumination is given as[^2]:

$$\text{MTF}_\text{ideal opt}(f) = \frac{2}{\pi} \left[ \arccos \left( \frac{f}{f_c} \right) - \frac{f}{f_c}  \sqrt{1- \left( \frac{f}{f_c} \right)^2} \right]$$

where $f$ is the input line frequency and $f_c$ is the [spatial cut-off frequency](imager_geom/#spatial-cut-off-frequency), which is a inversely proportional to the wavelength.

This can also be written as:

$$\text{MTF}_\text{ideal opt}(f) = \frac{2}{\pi} \left[ \arccos(\nu) - \nu \sqrt{1- (\nu)^2} \right]$$

where $\nu =\left( \frac{f}{f_c} \right) $.

An alternative formulation is:

$$\text{MTF}_\text{ideal opt}(f) = \frac{2}{\pi} \left[ \psi - \cos(\psi) \sin(\psi) \right]$$

where $\psi$ is equal to $\arccos(\nu)$.

The MTF value (or the curve) is wavelength dependent. The optical MTF is usually expressed either as multiple curves for each wavelength within the limits of the imager (for example red, green and blue for a colour imager) or as a single curve with the weighted average of multiple curves.

If there is an obstruction in front of the aperture (such as the secondary mirrors in most reflecting telescopes), the MTF curve is affected. It lowers the MTF in lower frequencies and increases it in the higher frequencies.

The real optical MTF will be lower than this value, due to real world design limitations, materials, manufacturing, integration as well as mechanical and thermal loads. This is usually simulated in a software like Zemax and eventually measured in the lab.

#### Aberration Transfer Factor (ATF) and Aberrated Optical MTF

The Aberration Transfer Factor[^3] (ATF) is an empirical model that combines all sources of optical aberrations into a single total wavefront error:

$$\text{ATF}(f) = 1- \left( \frac{W_{RMS}}{0.18} \right)^2 \left[ 1 - 4 ( \nu -0.5 )^2 \right] $$

where $W_{RMS}$ is the RMS of the total wavefront error, or how much the actual wavefront deviates from the ideal wavefront. The unit of this deviation is the multiple wavelengths. For $W_{RMS} = 0$, ATF is equal to 1 for all input frequencies, corresponding to no optical aberrations.

All aberration sources are (RMS) summed and then the resulting wavelength error can be inserted in the ATF to compute the total ATF of the optical system.

$$\text{MTF}_\text{aberr opt}(f) = \text{MTF}_\text{ideal opt}(f) \times \text{ATF}(f) $$  

Multiplying the ATF value with the ideal optical MTF, we can reach a more realistic MTF with the aberrations. As the $W_{RMS}$ value increases, the ATF value decreases and the resulting MTF also decreases, corresponding to a degradation in image quality.

Some sample fabrication tolerances are given [here](https://www.telescope-optics.net/fabrication.htm). For example, surface roughness for Commercial Optics can be a single wavelength (Peak-to-Valley), whereas for Precision Optics it could be about quarter of a wavelength and for High Precision Optics it could be as low as 5% of a wavelength. Satellite imagers would also be as high as 5% of a wavelength.

### Detector Sampling MTF

Each detector (or a single pixel of the detector) performs spatial averaging of the irradiance, or more precisely, we integrate the irradiance with the detector responsivity, over the detector area (See [here](https://spie.org/publications/spie-publication-resources/optipedia-free-optics-information/tt52_21_detector_footprint_mtf) for more information). In the frequency domain this corresponds to:

$$\text{MTF}_\text{detector}(f) = \frac{\sin(\pi p f)}{\pi p f} = \text{sinc}(p f)$$

where $p$ is equal to pixel pitch. Note the [normalised sinc function notation](https://en.wikipedia.org/wiki/Sinc_function) used.

As the frequency increases, and the wave periods become comparable to the detector pixel size, the pixels cannot represent the sine wave properly and there is a reduction in modulation. At a line frequency corresponding to the inverse of the pixel pitch, modulation goes down to zero, as the input sine wave is completely inside the pixel pitch. Even higher input frequencies will then be completely undersampled. This results in contrast reversal and MTF values will be negative.

### Diffusion MTF

Crosstalk in the detector due to electronics, see chap 13 pg 404 of a Sys eng approach

## Dynamic Contributors to the MTF

### Overview

Dynamic MTF is the broad name given to multiple sources of image sharpness loss, commonly due to a "shaking" of the imager or a relative motion between the imager and the platform. This results in the loss of sharpness; sometimes manifesting itself as blurring roughly equal in all directions and sometimes as a directional smear in the image.

Depending on the imaging setting, one or more of the following sources can be present:

- Jitter: High frequency random imager shake *during the exposure duration*, causing a blur on the image.
- Motion blur: Relative motion between the imager and the scene *during the exposure duration*, causing a smear on the image.
- Drift/smear: Directional imager drift *during the exposure duration*, causing a smear on the image.

Jitter is prominent in all imagers attached to equipment that generates high frequency vibration sources, such as engines in aircraft or reaction wheels in satellites. Cryocoolers with pistons also introduce vibrations close to the detector. Motion blur is prominent in satellite imagers at Low Earth Orbit (as opposed to Geosynchronous Orbit) and aircraft, both of which fly over the target. Drift/smear is prominent in imagers that are not physically fixed to a stable platform and/or with relatively long exposures where the slow shaking of the imaging platform becomes visible.

The common theme in all these sources is the exposure duration. Generally speaking, vibration sources with a period much shorter than (or a frequency higher than) the exposure duration are called jitter and they cause random artefacts such as blurring. The detector pixel captures information from around neighbouring areas in all directions, rather than the area corresponding to what the pixel would normally "see". Vibration sources with a frequency close to exposure duration cannot complete multiple cycles (or even a single cycle) during the exposure duration, therefore they cause a drift/smear in the direction they apply. The vibration sources with periods much longer than the exposure duration have no significant impact on the sharpness (or MTF).

It follows that the same vibration sources may cause different type of artefacts for different exposure durations. As the exposure duration becomes longer, even low frequency vibration sources start causing blur (rather than smear), and very slow motion sources start causing smear. For a 10 ms exposure (or "1/100" s in photography terms and 100 Hz in frequency terms), a 1000 Hz (1 ms) vibration source would cause jitter. A 50 Hz (20 ms) vibration source would cause a drift. A 1 Hz (1 s) vibration source would be too slow to impact the imaging quality (though the image shake would be visible in a video). If the exposure is increased to 100 ms (10 Hz), both 1000 Hz and 50 Hz vibration sources would cause jitter and the 1 Hz vibration source could cause drift/smear.

To summarise:

- For the vibration frequencies much lower than the exposure frequency, they are treated as bias, as the vibration is too slow to cause any blurring in the image.
- For the vibration frequencies higher than the exposure frequency, they are treated as random blurring, as the vibration completes at least one cycle and there are usually many more frequencies and phases.
- For the vibration frequencies lower than the exposure frequency (but still in the same order), the sine wave is completed only partially but it is very difficult to know the actual phase. Therefore, assuming a linear drift/smear is actually a worst-case assumption, with the drift rate approximating the "climb" part of the sine wave.

While it may then sound tempting to reduce the exposure duration to maximise the sharpness, it also controls the number of photons received, and therefore a lot of critical parameters, chiefly Signal-to-Noise Ratio (SNR): the longer is the exposure, usually the better is the SNR (up to saturation level). Therefore, the exposure duration needs to be optimised between the image sharpness and SNR needs.

The discussion above is valid for most usecases, where single images are generated that are independent from each other. However, if the images of the same location need to be combined to create another image (e.g., Time Delay Integration (TDI) applications), then the images should match exactly and the same pixel in the consecutive images should "see" the same area on the ground. In reality, this is not possible, and the jitter as well as drift/smear within the *total combined imaging duration* (e.g., the total TDI column duration) will also impact the total imaging sharpness.

Continuing from the example above with the 10 ms exposure duration, if we introduce a 5 stage TDI, then the imager should stay stable not only for just the 10 ms exposure duration, but also within the 50 ms total TDI column duration. The jitter and drift/smear within this total TDI column duration will also introduce a loss in sharpness.

When computing the *total* MTF for the system, all the relevant contributors should be combined and must be then multiplied with the Static or Imager MTF. Therefore all the contributors should be thought together and trade-offs will have to be considered for the overall Imaging System, including the Imager and the Platform that the Imager is mounted on. As an example, the following are some practical considerations when designing an Imaging System:

- Using an excellent imager on an unstabilised platform may result in mediocre image quality.
- Using a detector with better Quantum Efficiency (or optics with better transmission) may result in shorter exposure durations and/or TDI stages, resulting in better image sharpness.
- As the jitter and drift is measured in the percentage of the detector pixel size, binning the pixels will increase the effective pixel size and reduce the relative error, yielding sharper images, albeit at reduced image resolution.  

### Jitter

As explained above, jitter is due to the high frequency "shaking" of the imager (or more specifically, the line of sight (LoS) vectors), at a frequency that is higher than that of the imaging. It is usually caused by the shaking of the imaging platform (e.g. the aircraft or the ground vehicle that carries the imager) or the mechanical cooling system of the detector where present. The shaking usually involves multiple frequencies and multiple phases, therefore the combined effect can be modelled as a Gaussian variation of all LoS vectors at the same time, and in both horizontal and vertical directions - though the vibration profiles in the horizontal and vertical directions may differ. As this introduces irradiance into a pixel from neighbouring areas, this effect manifests itself as a blur (or a reduction in MTF) on the image. As mentioned above, the shorter the imaging or exposure duration, the higher the frequency that causes blurring and therefore the smaller the effects of jitter.

The jitter MTF is computed as a Gaussian of the form (from [^2] Vol 4 pg 69 eqn 2.6):

$$\text{MTF}_\text{jitter}(f) = exp(-2 (\pi  j_{LoS} f)^2)$$

where $j_{LoS}$ is the standard deviation (or 1 sigma) of the jitter amplitude (in the angular sense)) and $f$ is the input line frequency.

This equation can be written in terms of jitter defined as distance on the image plane (in pixels):

$$\text{MTF}_\text{jitter}(f) = exp(-2 (\pi  (j_{pix} p) f)^2)$$

where $j_{pix}$ is the standard deviation (or 1 sigma) of the jitter amplitude (as percentage of the pixel), $p$ is the pixel pitch and $f$ is the input line frequency.s

When multiple images are combined (as in TDI), then not only the blurring in the individual images that make up the final image is important, but also the motion during the combination of multiple images. Even if we have extremely short duration (and therefore sharp) images of the same scene, if the images are taken over a duration where jitter effects are significant, then some blurring will be introduced when all the images are combined. This is explained in greater detail in a [later section](#combining-multiple-images-the-effect-of-time-delay-and-integration-tdi).

### Motion Blur

Motion blur is caused by the relative motion between the scene and the imager during the exposure time. It is in the direction of relative motion.

In some cases, the imager is looking at a static scene (no relative motion), but a target within the scene may be moving with respect to the imager. If the target is moving appreciably during the exposure time (for example 0.1 pixels) then the target has a motion blur. In other cases, the imager may move with respect to the scene *by design*, usually in a scanning scheme. Examples could be aircraft or low altitude satellites flying over a scene. In this case, the motion blur is strictly in the scanning direction - usually in the alongtrack direction.

Focussing on the latter case, where the entire scene is moving with respect to the imager, motion blur introduces its own MTF. However, as the resulting blurring is directional, so is the MTF.

The MTF due to motion blur is given as:

$$\text{MTF}_\text{mot blur}(f) = \frac{\sin(\pi d p f)}{\pi p f} = \text{sinc}(d p f)$$

$d$ is the blur or drift, expressed in pixel ratio (e.g., 0.1 pixel drift/blur)

only in ALT direction for sats?
practical examples!!!

--------------------------------------------------

for the satellite pushbroom case
integ_time_ratio = timings.integration_duration /  (timings.frame_duration * channel.binning)
This means

a_fx = (pixel_pitch * input_line_freq / u.lp).to_reduced_units()

blur_in_pixels = integ_time_ratio
This means that the longest blur can be one pixel (x binning) long.

Another way could be angular velocity on the image = ang_vel x integ_time / ifov

sinc((blur_in_pixels*a_fx).m)

--------------------------------------------------

### Drift/Smear

practical examples!!!

- Yaw steering error (Single Image)

### Combining Multiple Images: The Effect of Time Delay and Integration (TDI)

In imaging, multiple acquisitions (usually taken in quick succession) are often combined at the detector level to increase the Signal-to-Noise Ratio (SNR). Specifically, TDI is a common technique, where a line detector (with multiple lines) sweep over the target in the along-track direction and successive lines of the detector image the same area. If 8 lines are used to acquire 8 images of the same area, this is called an 8-stage TDI system. However, the combination of these images or acquisitions should be an exact geometric match, otherwise smearing will occur. This requires that:

- for a multiple full frame, the relative image geometry should not change, therefore no relative velocity or rotation rate of the scene or target
- for a line scanning TDI, the relative image geometry should not change, therefore no relative velocity, rotation rate or around boresight rotation of the scene or target
- for a line scanning TDI, the timing should be exactly matching the relative motion in the along-track direction, such that the successive lines image exactly the same location

Clearly, these are not fully possible in a real system. For the first and the second, any relative *rotation rate* in along-track and across-track direction will manifest itself as a drift/smear in the respective directions for all pixels in the detector line. However, a small fixed rotation in along-track or across-track directions act as biases and do not introduce any smear during the multiple image acquisitions. On the other hand, if there is a *fixed rotation* around the detector normal, then the smear will be a combination of the along-track and across-track components. The higher the number of acquisitions (or TDI stages) combined, the longer the smearing. For example, for an 8 stage TDI, the smearing may be limited to 10% of the pixel size, whereas for a 64 stage TDI the smearing will span an entire pixel, and half of the data will have been acquired from the "neighbouring" area. For example, a common occurrence with the satellites is a fixed yaw-steering error, where the satellite yaws to compensate for the rotation of the Earth. In some cases this yaw angle is not precise enough, introducing a rotation with respect to the along-track direction. And in all cases the yaw angle can be optimised for a single pixel only (usually the centre pixel), resulting in a "residual yaw steering error", increasing for the pixels towards the edges.

The latter can be called a "Speed Mismatch" error, where the imaging Frame Rate or Line Rate does not match the relative along-track motion exactly. For example, for a satellite or an aircraft with a line scanner overflying a target, the line rate should match the ground speed. If the ground speed is 7000 m/s for a satellite with a Ground Sampling Distance of 7 m, then the Line Rate should be exactly 1 ms. Any mismatch will manifest itself as an along-track smearing on the image. As with the rotation problem above, the larger the number of images acquired (or TDI stages), the longer is the smear.

For this section we will assume that the x-axis on the detector is aligned with the along-track direction and y-axis aligned with the across-track direction. An $\epsilon$ rotation the boresight axis (or z-axis on the detector) results in a $\delta y$ pixel smear in the y-axis. $N$ is the number of TDI stages (or the number of x-axis pixels accumulated during TDI).

$$ \delta y = \tan(\epsilon) \times N $$

The smear due to Speed Mismatch Error is simply given as:

$$ \delta x = \text{SME} \times N $$

where $\delta x$ is the blur in the along-track direction due to Speed Mismatch. If the line rate is 10% higher or lower than what it should be, then at each line the target area moves 10% of a spatial sampling distance, eventually advancing by the equivalent of half a pixel ($\delta x$) at 5 stages ($N$).

It should be noted that, depending on the problem and the relative geometry, some pixels may experience more smear than others (for example "residual yaw steering error").

The drift/smear is more due to more systematic (and slower) errors and manifests itself as a linear smear on the image. However, any vibration that has a higher frequency than the total frequency of the TDI column will result in random "jumps" of the target area in successive images, introducing a blurring in the combined images.

Once the drift/smear due to TDI is quantified, the MTF equations are similar to the Drift/Smear previous sections. The jitter MTF due to TDI also uses an equation similar to the Jitter in the previous sections.  

Another important point to note is that, the time scale is the entire image acquisition (or the total TDI column duration), rather than the integration or exposure duration. Therefore, jitter and drift/smear disturbances for the TDI column are different than those for a single image acquisition. For example, consider a satellite pushbroom imager that uses a 0.5 ms exposure time, 1 ms line rate and 10 TDI stages. A 100 pix/s rotation rate that drifts the image in the across-track direction can introduce a 0.05 pixel smear in a single image, but over the entire 10 ms TDI duration, the overall smear is equal to 1 pixel. Similarly, the jitter over 0.5 ms will be much smaller than the jitter over 10 ms. Consequently, in addition to the blurriness of single acquisition images, the combined image will introduce another "layer" of blurriness.

To summarise, the following effects will introduce a drift/smear for a multiple image acquisition case:

- Relative rotation rate between the imager and the target (such as an attitude rate control error)
- A fixed rotation between the along-track direction and the detector lines (such as a yaw steering error or residual yaw steering error on a satellite)
- A mismatch between the relative motion and the image acquisition rate

Furthermore, "TDI jitter" will introduce a random blurring.

Finally, if some of the images are more blurry than others (for example due to "single image jitter" varying between images), then the combined image will have the average blur of all images.

## Atmospheric Contributors to the MTF

### Atmospheric Turbulence

TBW

### Aeorosol

TBW

[^1]: A Tutorial on Electro-Optical/Infrared (EO/IR) Theory and Systems; G. M. Koretsky, J. F. Nicoll, M. S. Taylor; Institute for Defense Analyses, IDA Document D-4642, 2013.

[^2]: The Infrared & Electro-Optical Systems Handbook; J. S. Accetta, David L. Shumaker (Ed.);  Infrared Information Analysis Center, 1993.

[^3]: The Art and Science of Optical Design; R. R. Shannon;  Cambridge University Press; 1997.
