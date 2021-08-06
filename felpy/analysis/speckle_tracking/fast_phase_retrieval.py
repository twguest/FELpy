# -*- coding: utf-8 -*-

import numpy as np
from math import floor

from numpy.fft import fft2, fftshift, ifft2, ifftshift

def fast_phase_retrieval(dx,dy,wav,z):
    """
    this function calculates the phase-shift of a set of wavefields given the 
    displacement fields dx,dy between phase-reconstructions.
    
    :param dx: horizontal component of the displacement field [m]
    :param dy: vertical component of the displacement field [m]
    :param wav: image wavelength [m]
    :param z: propagation distance [m]
    
    :returns delta_phi: phase-difference between image pair generating dx,dy
    """
    
    assert dx.shape == dy.shape, "displacement fields must be the same size"
    
    k = 2*np.pi/wav
    
    Nx, Ny = dx.shape
  
    dkx = 2 * np.pi / (Nx)
    dky = 2 * np.pi / (Ny)
    
    kx, ky = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dky, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dkx)
    ftfilt = 1 / (1j * kx - ky)
    ftfilt[np.isnan(ftfilt)] = 0
    
    C = (z/k)
    
    numerator = fft2(fftshift(dx + 1j*dy))
    delta_phi = C * ifftshift(ifft2((numerator)*(ftfilt)))
    
    return delta_phi    


