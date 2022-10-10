# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch

def power_spectral_density(array, spatial_sampling = 1, pulse_duration = 1):
    """
    plot the power spectral density of an array w/ respect to spatial freq.
    
    :param array: 2D complex valued array
    :param spatial_sampling: pixel size in m
    :param pulse_duration: in s, for conversion from energy to power. 
    """
    q = np.fft.fftshift(np.fft.fft2(array))

    s = abs(q)**2*1e-06
    s = np.diag(np.fliplr(s))
    freq = np.linspace(-s.shape[0]//2 * spatial_sampling,
                       s.shape[0]//2 * spatial_sampling,
                       s.shape[0])
    
    return freq[s.shape[0]//2:], s[s.shape[0]//2:]

if __name__ == '__main__':
    
    from matplotlib import pyplot as plt
    from felpy.model.src.coherent import construct_SA1_pulse

    from felpy.utils.vis_utils import Grids
    from matplotlib.ticker import FormatStrFormatter
    
    wfr = construct_SA1_pulse(512, 512, 2, 5, 1)

    ef = wfr.as_complex_array().sum(-1)
    ef *= np.random.rand(512, 512)

    freq, s = power_spectral_density(ef, spatial_sampling = wfr.dx,
                                     pulse_duration = wfr.pulse_duration)
    
    grid = Grids()
    grid.create_grid(1,1)
    ax = grid.axes
    ax.semilogy(freq, s)
    
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.e'))
    
    ax.set_ylabel("Power Spectral Density (W$m^{-3}s^{-1}$)")
    grid.set_fontsize(22)