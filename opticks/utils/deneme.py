from pathlib import Path

import numpy as np

from opticks import u
from opticks.imager_model.imager import Imager

if __name__ == "__main__":

    file_directory = Path("opticks", "temp", "sim")
    optics_file_path = file_directory.joinpath("optics.yaml")
    detector_file_path = file_directory.joinpath("detector.yaml")

    imager = Imager.from_yaml_file(optics_file_path, detector_file_path, None)

    optics = imager.optics
    detector = imager.detector
    timings = detector.params.timings
    # select the first channel
    channel = next(iter(detector.params.channels.all.values()))

    # ******************

    nyq_freq = channel.nyquist_freq( with_binning=True)
    print(f"Nyquist limit: {nyq_freq:.5}")
    
    # jitter is defined as 10% of binned pixel (1 sigma)
    jitter_in_pix = 0.1
    jitter_angular = jitter_in_pix * channel.ifov(optics, with_binning=True)

    print(f"allowable jitter (1 sigma): {jitter_angular.to("mdeg"):~P.5} ({jitter_angular.to("urad"):~P.5})")
    
    # smear is defined as 10% of binned pixel
    drift_in_pix = 0.1
    drift_angular = drift_in_pix * channel.ifov(optics, with_binning=True)

    print(f"allowable drift (max): {drift_angular.to("mdeg"):~P.5} ({drift_angular.to("urad"):~P.5})")

    # yaw steering error is defined as 10% of binned pixel
    yaw_steer_err_in_pix = 0.1
    yaw_steer_err_angular = np.arctan(yaw_steer_err_in_pix / channel.tdi_stages)
    
    # Ground velocity mismatch
    delta_v_err = 1/100  # percentage of effective ground velocity
    gnd_vel = 6.99808 * u.km / u.s
    altitude = 540 * u.km
    delta_v_err_ang_vel = delta_v_err * gnd_vel/altitude * u.rad

    print("-----------------------------------------------------")
    
    # compute plot data
    input_line_freq = 2.0 * u.lp / u.mm
    norm_line_freq = (input_line_freq / nyq_freq).to_reduced_units()
    
    # col AA (Optical Front-End Polychromatic)
    # Optical MTF
    mtf_optical = 97.5705941176471 /100
    print(f"MTF optical: {mtf_optical:.5%} (@ {input_line_freq:P~})")
    
    # col AB
    # Detector MTF 
    a_fx = norm_line_freq/2
    mtf_det_sampling = np.abs(np.sin(np.pi*a_fx)/(np.pi*a_fx))
    print(f"MTF detector sampling: {mtf_det_sampling:~.5%} (@ {input_line_freq:P~})")
    a_fx = (channel.pixel_pitch(with_binning=True) * input_line_freq / u.lp).to_reduced_units()
    # sinc does not receive Quantity input.
    mtf_det_sampling = np.sinc(a_fx.m)
    print(f"MTF detector sampling: {mtf_det_sampling:.5%} (@ {input_line_freq:P~})")
    
    # col AC = AA x AB
    # Imager MTF
    print(f"MTF imager: {(mtf_optical * mtf_det_sampling):.5%} (@ {input_line_freq:P~})")

    # col AD
    # Jitter (ALT and ACT)
    mtf_jitter = np.exp(-2*((np.pi*jitter_in_pix*norm_line_freq/2)**2))  
    print(f"MTF jitter: {mtf_jitter.m:.5%} (@ {input_line_freq:P~})")
    
    # col AE
    # Drift / Smear (ADCS roll, pitch, yaw errors) (ALT and ACT)
    mtf_drift = np.sin(np.pi*drift_in_pix*norm_line_freq/2)/(np.pi*drift_in_pix*norm_line_freq/2)
    print(f"MTF drift: {mtf_drift.m:.5%} (@ {input_line_freq:P~})")
    mtf_drift = np.sinc(drift_in_pix*a_fx.m)
    print(f"MTF drift: {mtf_drift:.5%} (@ {input_line_freq:P~})")
    
    # Yaw steering error (ACT)
    mtf_yaw_steer_err = np.sin(np.pi*yaw_steer_err_in_pix*norm_line_freq/2)/(np.pi*yaw_steer_err_in_pix*norm_line_freq/2)
    print(f"MTF yaw steer err: {mtf_yaw_steer_err.m:.5%} (@ {input_line_freq:P~})")
    
    # Gnd velociy mismatch (ALT)
    mtf_gnd_vel_mismatch = np.sin(np.pi* (delta_v_err * channel.tdi_stages) * norm_line_freq/2)/(np.pi* (delta_v_err * channel.tdi_stages) *norm_line_freq/2)
    print(f"MTF gnd vel mismatch: {mtf_gnd_vel_mismatch.m:.5%} (@ {input_line_freq:P~})")
    
    # col AG
    # Motion blur (Fwd linear motion) (ALT)
    integ_time_ratio =  timings.integration_duration /  (timings.frame_duration * channel.binning)
    detector.params.timings.frame_rate
    mtf_mot_blur = np.sin(np.pi*integ_time_ratio*norm_line_freq/2)/(np.pi*integ_time_ratio*norm_line_freq/2)
    print(f"MTF motion blur: {mtf_mot_blur.m:.5%} (@ {input_line_freq:P~})")
    
    print("-----------------------------------------------------")
    # col AH
    # Total ACT
    mtf_total_act = mtf_optical * mtf_det_sampling * mtf_jitter * mtf_drift 
    print(f"MTF total ACT: {mtf_total_act.m:.5%} (@ {input_line_freq:P~})")
    
    # col AI
    # Total ALT
    mtf_total_alt = mtf_optical * mtf_det_sampling * mtf_jitter * mtf_drift * mtf_gnd_vel_mismatch * mtf_mot_blur
    print(f"MTF total ALT: {mtf_total_alt.m:.5%} (@ {input_line_freq:P~})")
    
    print("-----------------------------------------------------")
    # ideal optical MTF (uniformly illuminated circ aperture, no significant aberrations)
    wavelength = 620 * u.nm # PAN
    # normalised optical frequency
    x = input_line_freq / optics.spatial_cutoff_freq(wavelength)     
    
    psi = np.arccos(x)
    mtf_ideal_optical = 2/np.pi * (psi-np.cos(psi)*np.sin(psi))
    print(f"MTF ideal optics: {mtf_ideal_optical.m:.5%} (@ {input_line_freq:P~})")
    
    # use this
    mtf_ideal_optical = 2/np.pi * (np.arccos(x)-x*np.sqrt(1-x**2))
    print(f"MTF ideal optics: {mtf_ideal_optical.m:.5%} (@ {input_line_freq:P~})")
