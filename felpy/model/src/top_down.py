# -*- coding: utf-8 -*-

import numpy as np
import seaborn as sns

from felpy.utils.np_utils import get_mesh
from felpy.utils.vis_utils import basic_plot, simple_line_plot, colorbar_plot, signal_plot, scatter_plot
from wpg.beamline import Beamline
from wpg.optical_elements import Aperture, Drift
from felpy.model.core.wavefront import Wavefront

from felpy.model.src.coherent import construct_SA1_wavefront

from matplotlib import pyplot as plt

from numpy import fft
from felpy.utils.opt_utils import ekev2wav

from scipy.constants import c

def gaussian_envelope(x, x0, width):
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


def generate_temporal_SASE_pulse(pulse_time, n_samples = 100, sigma = 4, VERBOSE = True):
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

    temporal_envelope = (1/np.sqrt(2*np.pi))*gaussian_envelope(t, 0, pulse_time)
        
    spectral_bw = 1/pulse_time
    w = 1/t
    spectral_envelope = gaussian_envelope(w, 0, spectral_bw)
    
    random_phases = np.random.uniform(-np.pi, np.pi, temporal_envelope.shape)

    E_t = fft.fft(spectral_envelope*np.exp(1j*random_phases))*temporal_envelope

    
    if VERBOSE:
        signal_plot(t, abs(E_t)**2,
                   xlabel = "Time (s)",
                   ylabel = "Intensity (a.u.)",
                   title = 'SASE Pulse',
                   context = 'talk')

    return E_t


def wavefront_from_array(cfr,nx,ny,nz,dx,dy,dz,ekev, pulse_duration = 40e-15):


    # Initialize empty wavefront.
    wfr = Wavefront()

    # Setup E-field.
    wfr.data.arrEhor = np.zeros(shape=(nx, ny, nz, 2))
    wfr.data.arrEver = np.zeros(shape=(nx, ny, nz, 2))

    wfr.params.wEFieldUnit = 'sqrt(W/mm^2)'
    wfr.params.photonEnergy = ekev * 1000
    wfr.params.wDomain = 'time'
    wfr.params.Mesh.nSlices = nz
    wfr.params.Mesh.nx = nx
    wfr.params.Mesh.ny = ny


    
    wfr.params.Mesh.sliceMin = -pulse_duration / 2.
    wfr.params.Mesh.sliceMax = pulse_duration / 2.

    range_x = dx*nx
    range_y = dy*ny

    wfr.params.Mesh.xMin = -range_x / 2.
    wfr.params.Mesh.xMax = range_x / 2.
    wfr.params.Mesh.yMin = -range_y / 2.
    wfr.params.Mesh.yMax = range_y / 2.


    wfr.data.arrEhor = complex_to_wpg(cfr)
    
    #wfr.set_electric_field_representation('f')
        
    return wfr

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


def complex_to_wpg(arr): ### converter
    new_arr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2], 2])
    new_arr[:,:,:,0] = arr.real
    new_arr[:,:,:,1] = arr.imag
    return new_arr


if __name__ == '__main__':

    pulse_time = 100e-15
    sigma = 4 
    S = 4
    Seff = S/sigma
    
    n_samples, sampling_interval_t = temporal_sampling_requirements(pulse_time, VERBOSE = True, S = S)
    sampling_interval_w = 1/sampling_interval_t

    n_samples *= sigma 
    
    t = np.arange(-pulse_time*4, pulse_time*4, pulse_time/n_samples)

    temporal_profile = generate_temporal_SASE_pulse(pulse_time = pulse_time,
                                                    n_samples = n_samples,
                                                    sigma = sigma,
                                                    VERBOSE = True)
    print("Number of Samples: ",n_samples)
    print(temporal_profile.dtype)
    spatial_profile = construct_SA1_wavefront(128, 128, 5.0, 0.25)

    wfr = wavefront_to_wavefield(spatial_profile, temporal_profile)
    
    from wpg.wpg_uti_wf import plot_intensity_map
    plot_intensity_map(wfr)
    wfr.set_electric_field_representation('f')
    from felpy.model.beamlines.exfel_spb.methods import setup_spb
    spb = setup_spb(parameter_file = "/opt/FELpy/felpy/data/params/spb-sfx_nkb_FAST.json", theta_KB = 5e-03, theta_HOM = 3.5e-03) #bl = spb.bl
    bl = spb.bl 
    
    ## EXEC
    bl.propagate_sequential(wfr)
