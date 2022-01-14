# -*- coding: utf-8 -*-

import numpy as np

from numpy import fft

from wpg import srwlib

from felpy.utils.opt_utils import ekev2k, geometric_focus
from felpy.model.backend.wpg_converters import wavefront_from_array
from felpy.model.source.coherent import modify_beam_divergence

FWHM2RMS = np.sqrt(8*np.log(2)) ### FWHM = sqrt(8ln(2))*sigma
FWHM2E2 = 1.66

def gaussian_profile(x, x0, width):
    """
    generate a gaussian envelope for modelling the integrated pulse_data

    :param x: 1D list of time/frequency positions
    :param x0: central time/frequency
    :param width: width of XFEL pulse
    """
    return np.exp(-np.power(x - x0, 2.) / (2 * np.power(width, 2.)))


def temporal_sampling_requirements(pulse_time, S = 10, VERBOSE = True):
    """
    calculate the number of samples required for the defined pulse length.

    from Partial-coherence method to model experimental free-electron laser pulse statistics - Pfeifer - 2021 -
    Optics Letters - Vol. 35 No 20.

    We expect the number of sampling intervals to satisfy
    |w_{i} - w_{i+1} << 2\pi/tau
    where tau is pulse time.

    :params pulse_time: length of the XFEL pulse
    :param S: sampling fadctor
    :returns n: number of sampling intervals req.

    NOTE: a<<b is set to satisfy a*10<b
    """
    freq_sampling = (2*np.pi)/(pulse_time)
    temporal_sampling = 1/(freq_sampling)
    n = np.ceil(S*pulse_time/(temporal_sampling))

    if VERBOSE:
        print("Frequency Sampling Interval: {:.2e} Hz".format(freq_sampling))
        print("Temporal Sampling Interval: {:.2e} s".format(temporal_sampling))
        print("Number of Req. Samples: {}".format(n))

    return int(n), temporal_sampling


def generate_temporal_SASE_pulse(pulse_time, n_samples = 100, sigma = 4):
    """
    generate a single SASE pulse

    - assumes that the spectral and temporal bandwidth are the pulse are reciprocal,
    - assumes that the spectral and temporal profiles are gaussian,
    - assumes that there is no jitter from the central value of the Gaussian (ie, the pulse profile
    is persistent).

    in future, we will extend this extned the functionality to account for non-Gaussian profiles,
    and pulse-length jitter.

    it will also be useful to trim the function to the relevant regions (and thus reduce the number of points)

    :param pulse_time: expectation value of the SASE pulse time
    :param VERBOSE: [bool] enables printing and plotting
    """
    
    t = np.linspace(-pulse_time*sigma, pulse_time*sigma, n_samples)

    temporal_envelope = (1/np.sqrt(2*np.pi))*gaussian_profile(t, 0, pulse_time)
        
    spectral_bw = 1/pulse_time
    w = 1/t
    spectral_envelope = gaussian_profile(w, 0, spectral_bw)
    
    random_phases = np.random.uniform(-np.pi, np.pi, temporal_envelope.shape)

    E_t = fft.fft(spectral_envelope*np.exp(1j*random_phases))*temporal_envelope
 

    return E_t

def wavefront_to_wavefield(spatial_wfr, temporal_profile):
    
    new_wfr = spatial_wfr.as_complex_array()[:,:,:]*temporal_profile


    dx, dy = spatial_wfr.get_spatial_resolution()
    dz = spatial_wfr.get_temporal_resolution()
    

    wfr = wavefront_from_array(new_wfr,
                             nx = spatial_wfr.params.Mesh.nx,
                             ny = spatial_wfr.params.Mesh.ny,
                             nz = len(temporal_profile),
                             dx = dx,
                             dy = dy,
                             dz = dz,
                             ekev = spatial_wfr.params.photonEnergy/1000)
    return wfr

def gaussian_envelope(nx, ny,
                      xMin, xMax,
                      yMin, yMax,
                      fwhm):
    """

    :param nx: DESCRIPTION
    :param ny: DESCRIPTION
    :param dx: DESCRIPTION
    :param dy: DESCRIPTION
    :param fwhm: DESCRIPTION
    :param divergence: DESCRIPTION
    
    :return gaussian: DESCRIPTION
    """
     
    # Initializing value of x-axis and y-axis
    # in the range -1 to 1
    x, y = np.meshgrid(np.linspace(xMin, xMax, nx), np.linspace(yMin, yMax, ny))
    dst = np.sqrt(x*x+y*y)
     
    # Initializing sigma and muu
    sigma = fwhm/FWHM2E2
    muu = 0.000
     
    # Calculating Gaussian array
    gaussian = np.exp(-( (dst-muu)**2 / (2*sigma**2 ) ) )

    return gaussian 

def spherical_phase(nx, ny,
                      xMin, xMax,
                      yMin, yMax,
                      fwhm, divergence,
                      ekev):
    
 
    
    k = ekev2k(ekev)
    f = geometric_focus(fwhm, divergence) 
    print(f,k)
    # Initializing value of x-axis and y-axis
    # in the range -1 to 1
    x, y = np.meshgrid(np.linspace(xMin, xMax, nx), np.linspace(yMin, yMax, ny))
    dst = np.sqrt(x*x+y*y)
     
    # Initializing sigma and muu
    muu = 0.000
     
    # Calculating Gaussian array
    spherical_phase_term = np.exp(1j*k*(x**2+y**2)/(2*f))

    return spherical_phase_term 
    
def spherical_phase_sampling():
    pass
def complex_gaussian_envelope(nx, ny,
                      xMin, xMax,
                      yMin, yMax,
                      fwhm, divergence,
                      ekev):
    """
    
    :param nx: DESCRIPTION
    :param ny: DESCRIPTION
    :param dx: DESCRIPTION
    :param dy: DESCRIPTION
    :param fwhm: DESCRIPTION
    :param divergence: DESCRIPTION
    
    :return gaussian: DESCRIPTION

    """

    gaussian = gaussian_envelope(nx, ny, xMin, xMax, yMin, yMax, fwhm).astype('complex128')
    
    gaussian*= spherical_phase(nx, ny, xMin, xMax, yMin, yMax,
                          fwhm, divergence,
                          ekev)
    
    return gaussian
    


def generate_coherent_source_wavefront(nx, ny, fwhm, divergence):
    
    gaussian_1d = gaussian_envelope(nx, ny, fwhm, divergence)



if __name__ == '__main__':
    
    from wpg.wpg_uti_wf import plot_intensity_map, plot_intensity_qmap
    
    nx = 512
    ny = 512
    
    xMin = yMin = -400e-06
    xMax = yMax =  400e-06
    
    for fwhm in np.linspace(50e-06, 150e-06, 10):
        
        divergence = 11e-06
        ekev = 10
        
        env = complex_gaussian_envelope(nx, ny,
                                        xMin, xMax, yMin, yMax,
                                        fwhm,
                                        divergence,
                                        ekev)
        
        pulse_time = 55e-15
        sigma = 4
        S = 4
        Seff = S/sigma
        
        n_samples, sampling_interval_t = temporal_sampling_requirements(pulse_time, VERBOSE = True, S = S)
        sampling_interval_w = 1/sampling_interval_t
    
        n_samples *= sigma 
        #n_samples = 10
        t = np.arange(-pulse_time, pulse_time, pulse_time/n_samples)
    
        temporal_profile = generate_temporal_SASE_pulse(pulse_time = pulse_time,
                                                        n_samples = n_samples,
                                                        sigma = sigma)
        env = env[:,:,np.newaxis]*temporal_profile
        
        wfr = wavefront_from_array(env, nx = nx, ny = ny,
                                  nz = n_samples, dx = (xMax-xMin)/nx,
                                  dy = (yMax-yMin)/ny, dz = (pulse_time/n_samples),
                                  ekev = ekev, pulse_duration = pulse_time)
        
        
        #modify_beam_divergence(wfr, wfr.get_fwhm()[0], divergence)
        
        plot_intensity_map(wfr)
        print(wfr.get_fwhm())
        print(wfr.get_divergence())
        plot_intensity_map(wfr)
        #srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
    # =============================================================================
    #     plot_intensity_map(wfr)
    #     print(wfr.get_divergence())
    #     plot_intensity_map(wfr)
    # 
    #     #print(wfr.get_spatial_resolution())
    #     modify_beam_divergence(wfr, fwhm, divergence)
    #     #print(wfr)
    # 
    #     #env = env[:,:,np.newaxis] * temporal_profile
    # =============================================================================
        
