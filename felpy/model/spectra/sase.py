import numpy as np

from numpy import fft



from scipy.constants import h,c,e


h_eV_s = h/e
hr_eV_s = h_eV_s/(2*np.pi)


def gaussian_profile(x, x0, width):
    """
    generate a gaussian envelope for modelling the integrated pulse_data

    :param x: 1D list of time/frequency positions
    :param x0: central time/frequency
    :param width: width of XFEL pulse
    """
    return np.exp(-np.power(x - x0, 2.) / (2 * np.power(width, 2.)))

def linear_SASE_spectrum(pulse_duration, E0, dE = 1e-04, t0 = 0, n_samples = 500, sigma = 3):
    """
    generate a single SASE pulse profiles

    - assumes that the spectral and temporal profiles are gaussian,
    - assumes that there is no jitter from the central value of the Gaussian (ie, the pulse profile
    is persistent).

    in future, we will extend this extned the functionality to account for non-Gaussian profiles,
    and pulse-length jitter.

    it will also be useful to trim the function to the relevant regions (and thus reduce the number of points)

    :param pulse_duration: expectation fwhm value of the SASE pulse time [in s]
    :param E0: central energy of pulse [in eV]
    :param dE: energy bandwidth / relative of the pulse
    :param t0: center of intensity distribution in time (default 0s)
    :param n_samples: number of samples in time and freq domains.
    :param sigma: number of spectral bandwidth boundaries to define energy-domain axis
    
    
    :returns t: time-axis
    :returns E_t: SASE intensity profile in time-domain
    :returns E: energy-axis
    :returns E_eV: SASE intensity profile in energy-domain    
    """
    
    ### convert fwhm to sigma
    dE = E0*dE

    
    pulse_duration = np.sqrt(2) * pulse_duration / (2 * np.sqrt(2 * np.log(2)))  
    dE = np.sqrt(2) * dE / (2 * np.sqrt(2 * np.log(2)))

    E = np.linspace(E0-dE*sigma, E0+dE*sigma, n_samples) ### define frequency/energy domain
    estep = (E[1] - E[0])/hr_eV_s ### step size of freq domain

    t = np.linspace(-np.pi / estep, np.pi / estep, n_samples) ### define time-domain

    temporal_envelope = (1/np.sqrt(2*np.pi)) * \
        gaussian_profile(t, t0, pulse_duration)
    
    spectral_envelope = gaussian_profile(E, E0, dE)

    random_phases = np.random.uniform(-np.pi, np.pi, temporal_envelope.shape)

    E_t = fft.fftshift(fft.fft(spectral_envelope*np.exp(-1j*random_phases)))*temporal_envelope
    E_eV = np.fft.ifft(E_t)
    
    return t, E_t, E, E_eV