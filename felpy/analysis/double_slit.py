import numpy as np
 
from matplotlib import pyplot as plt
from felpy.utils.opt_utils import ekev2wav

from numpy import sin,cos
 

def envelope(x, sw, wav, z):
    """
    envelope function of a double-slit
    
    :param x: horizontal coordinates
    :param sw: slit-seperation
    :param wav: photon wavelength
    :param z: distance to screen
    """

    a = (np.pi*sw*x)/(wav*z)
    A = (sin(a)/a)**2

    return A

def double_slit_interference(x, d, sw, px, wav, z, mu = 1):
    """
    returns the double-slit interference pattern

    :param x: horizontal coordinates
    :param d: slit-seperation
    :param sw: slit-width
    :param px: pixel-size
    :param wav: photon wavelength
    :param z: distance to screen
    :param mu: absolute complex degree of coherence
    
    """
    b = (np.pi*px*d)/(wav*z)
    B = sin(b)/b

    c = cos((2*np.pi*d*x)/(wav*z))
    C = B*c
    
    A = envelope(x,sw,wav,z)
    
    return A*(1+mu*C)

def finite_pixel_width(px,d,wav,z):
    """
    returns the sampling function of the detector

    :param px: pixel-size
    :param d: slit-seperation
    :param wav: photon-wavelength
    :param z: distance to screen
    """
    b = (np.pi*px*d)/(wav*z)
    return sin(b)/b