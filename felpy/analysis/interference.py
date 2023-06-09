"""
Sets of utilities for interference and interferometry measurements
"""

import numpy as np


def double_slit(x, ss, sw, px, wav, z, mu = 1):
    """
    regular double slit pattern across a screen detector a distance z (m) upstream
    
    :param x: array of horizontal coordinates (m)
    :param ss: slit_seperation
    :param sw: slit width
    :param px: pixel-size
    :param wav: wavelength
    :param z: distance between slit and detector (m)
    :param mu: absolute mutual coherence b/w field between slits
    """
    b = (np.pi*px*ss)/(wav*z)
    B = sin(b)/b

    c = cos((2*np.pi*ss*x)/(wav*z))
    C = B*c
    
    A = envelope(x,sw,wav,z)
    
    return A*(1+mu*C)
