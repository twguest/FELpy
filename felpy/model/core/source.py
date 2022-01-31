# -*- coding: utf-8 -*-

import numpy as np

from numpy import fft

from wpg import srwlib

from felpy.utils.opt_utils import ekev2k, geometric_focus
from felpy.model.backend.wpg_converters import wavefront_from_array
from felpy.model.source.coherent import modify_beam_divergence
from felpy.utils.opt_utils import ekev2wav, ekev2k
from felpy.utils.os_utils import timing
from wpg.wpg_uti_wf import calc_pulse_energy, plot_intensity_map

from felpy.model.source.SA1 import analytical_pulse_divergence,analytical_pulse_energy,analytical_pulse_width,analytical_pulse_duration

FWHM2RMS = np.sqrt(8*np.log(2)) ### FWHM = sqrt(8ln(2))*sigma
FWHM2E2 = 1.66

def get_pulse_energy(wfr):
    wfr.set_electric_field_representation('t')
    energy = calc_pulse_energy(wfr)[0]
    wfr.set_electric_field_representation('f')
    return energy
        
class Source:
    """
    generalised source class
    """
    def __init__(self, ekev, fwhm, divergence,
                 nx, ny, nz, xMin, xMax, yMin, yMax, pulse_duration,
                 zDomain = 'frequency', energy = 1,
                 temporal_profile = None):
        
        self.ekev = ekev
        self.wavelength = ekev2wav(ekev)
        self.k = ekev2k(ekev)
        
        self.divergence = divergence
        self.fwhm = fwhm
        
        self.pulse_duration = pulse_duration
        
        self.nx = nx
        self.ny = ny
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.zMin = -pulse_duration/2
        self.zMax = pulse_duration/2
        
        self.zDomain = 'frequency'
        
        self.energy = energy
        self.wfr = None
        
        self.temporal_profile = temporal_profile
    
    def set_beam_energy(self,  beam_energy):
        
        if self.wfr is None:
            assert("this operation acts on the wpg wavefront class")
            assert("please ensure self.wfr exists, then re-run")
        else:
            
            E0 = get_pulse_energy(self.wfr)
            self.wfr.data.arrEhor /= np.sqrt(E0/beam_energy)
            
            #print(get_pulse_energy(self.wfr))
            
            
class SA1_Source(Source):
    
    def __init__(self, ekev, q,
                 nx, ny,
                 xMin = -200e-06, xMax = 200e-06, yMin =-200e-06, yMax = 200e-06, S = 4):
        
        self.q = q
        self.S = S
        
        #divergence = analytical_pulse_divergence(ekev, 'mean')
        divergence = np.random.uniform(low = analytical_pulse_divergence(ekev, 'lower'), high = analytical_pulse_divergence(ekev, 'upper'))
        #print("Expected Divergence: {}".format(divergence))

        energy = analytical_pulse_energy(q, ekev)

        fwhm = analytical_pulse_width(ekev) 
        #print("Expected Size: {}".format(fwhm))

        pulse_duration = analytical_pulse_duration(q)
        #print("Expected Duration: {}".format(pulse_duration))

        
        
        super().__init__(ekev = ekev, fwhm = fwhm, divergence = divergence,
                         nx = nx, ny = ny,nz = 1,
                         xMin = xMin, xMax = xMax,
                         yMin = yMin, yMax = yMax,
                         pulse_duration = pulse_duration,
                         zDomain ='frequency',
                         energy = energy)
        
        self.get_wfr()
        

    #@timing
    def get_temporal_profile(self, sigma = 4, refresh = False):
        
        if self.temporal_profile is None:
            
            n_samples, sampling_interval_t = temporal_sampling_requirements(self.pulse_duration, VERBOSE = True, S = self.S)
         
            n_samples *= sigma 
            self.nz = n_samples
            self.temporal_profile = generate_temporal_SASE_pulse(pulse_time = self.pulse_duration,
                                                                 n_samples = n_samples,
                                                                 sigma = sigma)
        elif refresh == True:


            self.temporal_profile = generate_temporal_SASE_pulse(pulse_time = self.pulse_duration,
                                                                 n_samples = self.nz,
                                                                 sigma = sigma)
            
        return self.temporal_profile
        
    
    def get_wfr(self):
        
        if self.wfr is None:
        
            env = complex_gaussian_envelope(self.nx, self.ny,
                                            self.xMin, self.xMax, self.yMin, self.yMax,
                                            self.fwhm,
                                            self.divergence,
                                            self.ekev)
               

             
            #sampling_interval_w = 1/sampling_interval_t
      
            self.get_temporal_profile()
            
            env = env[:,:,np.newaxis]*self.temporal_profile
            
            self.wfr = wavefront_from_array(env, nx = self.nx, ny = self.ny,
                                      nz = len(self.temporal_profile), dx = (self.xMax-self.xMin)/self.nx,
                                  dy = (self.yMax-self.yMin)/self.ny, dz = (self.pulse_duration/len(self.temporal_profile)),
                                  ekev = self.ekev, pulse_duration = self.pulse_duration)
            
            self.set_beam_energy(self.energy)
            
            return self.wfr
        else: 
            return self.wfr
    
    
    
def gaussian_profile(x, x0, width):
    """
    generate a gaussian envelope for modelling the integrated pulse_data

    :param x: 1D list of time/frequency positions
    :param x0: central time/frequency
    :param width: width of XFEL pulse
    """
    return np.exp(-np.power(x - x0, 2.) / (2 * np.power(width, 2.)))


def temporal_sampling_requirements(pulse_time, S = 10, VERBOSE = False):
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
        #print("Frequency Sampling Interval: {:.2e} Hz".format(freq_sampling))
        #print("Temporal Sampling Interval: {:.2e} s".format(temporal_sampling))
        #print("Number of Req. Samples: {}".format(n))
        pass
    
    return int(n), temporal_sampling


def generate_temporal_SASE_pulse(pulse_time, n_samples = 100, sigma = 4, t0 = 0):
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

    temporal_envelope = (1/np.sqrt(2*np.pi))*gaussian_profile(t, t0, pulse_time)
        
    spectral_bw = 1/pulse_time
    w = 1/t
    
    if t0 == 0:
        w0 = t0
    elif t0 != 0:
        w0 = 1/t0
        
    spectral_envelope = gaussian_profile(w, w0, spectral_bw)
    
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
                      fwhm, x0 = 0, y0 = 0):
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
 
     
    # Calculating Gaussian array
    gaussian = np.exp(-(((x-x0)**2+(y-y0)**2) / (2*sigma**2 ) ) )

    return gaussian 

def spherical_phase(nx, ny,
                      xMin, xMax,
                      yMin, yMax,
                      fwhm, divergence,
                      ekev, x0 = 0, y0 = 0):
    
 
    
    k = ekev2k(ekev)
    f = geometric_focus(fwhm, divergence) 

    # Initializing value of x-axis and y-axis
    # in the range -1 to 1
    x, y = np.meshgrid(np.linspace(xMin, xMax, nx), np.linspace(yMin, yMax, ny))
 
    # Initializing sigma and muu
    
     
    # Calculating Gaussian array
    spherical_phase_term = np.exp(1j*k*((x-x0)**2+(y-y0)**2)/(2*f))

    return spherical_phase_term 
    
def spherical_phase_sampling():
    pass

def complex_gaussian_envelope(nx, ny,
                      xMin, xMax,
                      yMin, yMax,
                      fwhm, divergence,
                      ekev, x0 = 0, y0 = 0):
    """
    
    :param nx: DESCRIPTION
    :param ny: DESCRIPTION
    :param dx: DESCRIPTION
    :param dy: DESCRIPTION
    :param fwhm: DESCRIPTION
    :param divergence: DESCRIPTION
    
    :return gaussian: DESCRIPTION

    """

    gaussian = gaussian_envelope(nx, ny, xMin, xMax, yMin, yMax, fwhm, x0 = x0, y0 = y0).astype('complex128')
    
    gaussian*= spherical_phase(nx, ny, xMin, xMax, yMin, yMax,
                          fwhm, divergence,
                          ekev, x0 = x0, y0 = y0)
    
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
        
        divergence = 25e-06
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
        
        n_samples, sampling_interval_t = temporal_sampling_requirements(pulse_time, VERBOSE = False, S = S)
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
        
        print(wfr.get_divergence())
        
        #modify_beam_divergence(wfr, wfr.get_fwhm()[0], divergence)
        
        plot_intensity_map(wfr)

        
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
        
