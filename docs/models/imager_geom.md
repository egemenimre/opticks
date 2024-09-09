# Imager Geometry Models

## F-number

Computed as:

$$\frac{\text{focal length}}{\text{aperture diameter}}$$

## Full Optical FoV

This is the optical FoV. Image FoV is a function of the detector size on the focal plane. But this is the maximum possible FoV for the given optics.

$$2 \times \arctan\left( \frac{0.5 \times \text{image diameter on focal plane}}{\text{focal length}} \right) $$

## Aperture Area

Computed as:

$$\pi \times \left( \frac{ \text{aperture diameter}}{2} \right)^2$$

## Aperture Solid Angle

Aperture solid angle in steradians.

Computed as:

$$\frac{\pi}{4} \frac{\text{aperture diameter}^2}{\text{focal length}^2}$$

## Spatial Cut-off Frequency

Spatial cut-off frequency, assumes perfect incoherent optics.

Determines the theoretical limit of the optical resolution, or the smallest object resolvable by the optical system.

Computed as:

$$ 1 \over {\lambda F_\#} $$

in cycles per mm.

## IFOV

Instantaneous field of view (works in vertical and horizontal).

Assumes constant IFOV per pixel.

$$2 \times \arctan{ \left( \frac{\text{pixel pitch} }{2 \times \text{focal length}} \right)}$$

Pixel pitch (and the resulting IFOV) may or may not be binned, depending on the application.

## Pixel Solid Angle

Pixel solid angle (of a pyramid).

$$4 \times \arcsin \left(\sin(ifov / 2) \times \sin(ifov) / 2 \right) $$

IFOV (and the resulting Pixel Solid Angle) may or may not be binned, depending on the application.

## Horizontal and Vertical FoV

Assumes constant IFOV per pixel.

$$2 \times \tan \left( \frac{ifov \times pixels}{2} \right)$$

where `pixels` are the horizontal or vertical pixels used.

## Nyquist Frequency or Nyquist Limit

Nyquist frequency is defined as:

$$\frac{1 \ \text{line pair}}{2 \times \text{pix pitch}}$$

which translates to one line pair per two pixels (e.g., one pixel white, next one black). This is the absolute maximum limiting resolution of the sensor.
