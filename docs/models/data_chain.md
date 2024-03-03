# Data Chain Models

1. Pixel Read Rate: Pixel read rate out of the detector
2. Data Write Rate: Data write rate after encoding, framing and compression.

This chain is repeated (usually in parallel) for every channel or band.


## Pixel Read Rate

Pixel read rate out of the detector. (in Mpixel/s) (per channel)

Computed as:
 - Pushbroom type detector (with binning if applicable): 

    $$\text{horizontal pixels (binned)} \times \text{TDI stages} \times \text{line rate}$$

 - Full frame type detector (with binning if applicable): 

    $$\text{horizontal pixels (binned)} \times \text{vertical pixels (binned)} \times \text{frame rate}$$

The computation above can be modified for the unbinned pixels, depending on where (before or after binning) the pixel read rate is evaluated.

Note that the unused pixels are also read, this assumes that the detector does not have ROI functionality.

## Data Write Rate

Data write rate after quantisation, framing and compression.  (in Mbit/s) (per channel)

Firstly, data is converted from Analog to Digital with a certain quantisation level (e.g., 12 bit/pix):

$$ \text{encoded data rate} = \text{pixel read rate} \times \text{encoding} $$

Afterwards, the encoded data is processed (non-uniformity correction etc.) and then compressed:

$$ \text{processed data rate} = \text{encoded data rate} / \text{compression} $$

Finally, the processed data is written with some overhead:

$$ \text{write data rate} = \text{processed read rate} \times (1 + \text{overhead percentage}) $$
