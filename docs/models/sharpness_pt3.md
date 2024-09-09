# Dynamic Contributors to the MTF

The topic of sharpness and contrast performance is continued from [here](sharpness_pt2.md).

## Overview

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

The discussion above is valid for most usecases, where single images are generated that are independent of each other. However, if the images of the same location need to be combined to create another image (e.g., Time Delay Integration (TDI) applications), then the images should match exactly and the same pixel in the consecutive images should "see" the same area on the ground. In reality, this is not possible, and the jitter as well as drift/smear within the *total combined imaging duration* (e.g., the total TDI column duration) will also impact the total imaging sharpness.

Continuing from the example above with the 10 ms exposure duration, if we introduce a 5 stage TDI, then the imager should stay stable not only for just the 10 ms exposure duration, but also within the 50 ms total TDI column duration. The jitter and drift/smear within this total TDI column duration will also introduce a loss in sharpness.

When computing the *total* MTF for the system, all the relevant contributors should be combined and must be then multiplied with the Static or Imager MTF. Therefore all the contributors should be thought together and trade-offs will have to be considered for the overall Imaging System, including the Imager and the Platform that the Imager is mounted on. As an example, the following are some practical considerations when designing an Imaging System:

- Using an excellent imager on an unstabilised platform may result in mediocre image quality.
- Using a detector with better Quantum Efficiency (or optics with better transmission) may result in shorter exposure durations and/or TDI stages, resulting in better image sharpness.
- As the jitter and drift is measured in the percentage of the detector pixel size, binning the pixels will increase the effective pixel size and reduce the relative error, yielding sharper images, albeit at reduced image resolution.  

## Jitter

As explained above, jitter is due to the high frequency "shaking" of the imager (or more specifically, the line of sight (LoS) vectors), at a frequency that is higher than that of the imaging. It is usually caused by the shaking of the imaging platform (e.g. the aircraft or the ground vehicle that carries the imager) or the mechanical cooling system of the detector where present. The shaking usually involves multiple frequencies and multiple phases, therefore the combined effect can be modelled as a Gaussian variation of all LoS vectors at the same time, and in both horizontal and vertical directions - though the vibration profiles in the horizontal and vertical directions may differ. As this introduces irradiance into a pixel from neighbouring areas, this effect manifests itself as a blur (or a reduction in MTF) on the image. As mentioned above, the shorter the imaging or exposure duration, the higher the frequency that causes blurring and therefore the smaller the effects of jitter.

The jitter MTF is computed as a Gaussian of the form (from [^2] Vol 4 pg 69 eqn 2.6):

$$\text{MTF}_\text{jitter}(f) = exp(-2 (\pi  j_{LoS} f)^2)$$

where $j_{LoS}$ is the standard deviation (or 1 sigma) of the jitter amplitude (in the angular sense) and $f$ is the input line frequency.

This equation can also be written in terms of jitter defined as distance on the image plane (in pixels):

$$\text{MTF}_\text{jitter}(f) = exp(-2 (\pi  (j_{pix} p) f)^2)$$

where $j_{pix}$ is the standard deviation (or 1 sigma) of the jitter amplitude (as percentage of the pixel), $p$ is the pixel pitch and $f$ is the input line frequency.s

When multiple images are combined (as in TDI), then not only the blurring in the individual images that make up the final image is important, but also the motion during the combination of multiple images. Even if we have extremely short duration (and therefore sharp) images of the same scene, if the images are taken over a duration where jitter effects are significant, then some blurring will be introduced when all the images are combined. This is explained in greater detail in a [later section](#combining-multiple-images-the-effect-of-time-delay-and-integration-tdi).

## Motion Blur

Motion blur is caused by the relative motion between the scene and the imager during the exposure time. It is in the direction of relative motion. For this text we will assume a simple linear motion, though for long exposures, or very fast relative motion, this may not be correct.

In some cases, the imager is looking at a static scene (no relative motion), but a target within the scene may be moving with respect to the imager. If the target is moving appreciably during the exposure time (for example 0.1 pixels) then the target has a motion blur. In other cases, the imager may move with respect to the scene *by design*, usually in a scanning scheme. Examples could be aircraft or low altitude satellites flying over a scene. In this case, the motion blur is strictly in the scanning direction - usually in the along-track direction.

Focussing on the latter case, where the entire scene is moving with respect to the imager, motion blur introduces its own MTF contributor. As the resulting blurring is directional, so is the MTF.

The MTF due to motion blur is given as ([^4] Equation 14.3.3):

$$\text{MTF}_\text{mot blur}(f) = \frac{\sin(\pi d p f)}{\pi p f} = \text{sinc}(d p f)$$

where $d$ is the blur or drift extent, expressed in pixel ratio (e.g., 0.1 pixel drift/blur).

The blur extent for a scanner can be computed by multiplying the exposure time with the motion rate on the image plane. For the ideal case of a pushbroom imager with no oversampling and square pixels in the resulting image, line duration corresponds to one pixel. Therefore the ratio of the exposure duration to the line duration defines the blur extent as a percentage of a pixel. For the more general case, the ratio of the velocity on the target plane (for example ground velocity) to the distance of the imager to the scene (for example altitude) is equal to the ratio of the velocity on the image plane to the focal length. Using this, we can compute the blur extent:

$$ d = \text{velocity on image plane} \times \text{integ time} = \left( \frac{(\text{velocity on the target plane}) x (\text{focal length})} {\text{distance to target}} \right) \times \text{integ time} $$

For a satellite with a ground velocity of 6700 m/s at an altitude of 670 km and an imager with an effective focal length of 500 mm, the scene moves in the image plane at a rate of 5 mm/s or 5 µm/ms.

For scanning or rotating imagers, the blur extent is equal to the angle scanned over integration time over the angle corresponding to a single pixel (or instantaneous field of view):

$$ d = \frac{\text{ang vel} \times \text{integ time}}{\text{ifov}}$$

When converting the blur extent from distance on image plane to number of pixels, binning should also be taken into account. For a full-frame detector where image is generated by reading the entire frame, a x2 binning would reduce a full (unbinned) pixel blur to only half a (binned) pixel. For the pushbroom it works differently, where a x2 binning in along-track direction would still require two exposure times - one for each line. If the motion blur spans an entire (unbinned) pixel, then a x2 binning would mean x2 the motion blur, and therefore a full binned pixel.

## Drift/Smear

Drift/Smear is the more general form of the Motion Blur. Motion Blur deals with the target/scene motion or the imager motion (or both) during a single exposure time, without any disturbances. In reality, there are a lot of disturbances in most imagers, usually in the form of platform shaking or drifting, in a frequency that is lower than that of imaging.

The basic equations are the same as those in Motion Blur, though this time the causes are different. And, unlike Motion Blur, the blur extent direction does not have to be in the direction of the platform relative motion with respect to the scene.

$$\text{MTF}_\text{mot blur}(f) = \frac{\sin(\pi d p f)}{\pi p f} = \text{sinc}(d p f)$$

where $d$ is the blur or drift extent, expressed in pixel ratio. In this case, the blur extent should be computed in along-track and across-track axes separately, summing all the sources that introduce a drift rate over the integration time. Each disturbance is assigned a rotation/translation rate and a rotation/translation axis and, multiplying with the exposure time, the total blur extent on the image plane is computed.

## Combining Multiple Images: The Effect of Time Delay and Integration (TDI)

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

---

The topic of sharpness and contrast performance is continued [here](sharpness_pt4.md).

[^2]: The Infrared & Electro-Optical Systems Handbook; J. S. Accetta, David L. Shumaker (Ed.); Infrared Information Analysis Center, 1993.

[^4]: A System Engineering Approach to Imaging; N. S. Kopeika; SPIE Press, 1998.
