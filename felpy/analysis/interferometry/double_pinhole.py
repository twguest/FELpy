import math
import scipy

import numpy as np

from scipy.optimize import curve_fit
from felpy.utils.opt_utils import ekev2wav
from functools import partial
from sklearn.metrics import r2_score

def Airy(x, w, wavelength, z, d):
    Z = (math.pi / wavelength) * w * (x - d) / z
 
    f =  (math.pi / wavelength) * w ** 2 * 1 / z * scipy.special.jv(1, Z) / (Z)
    return f

def diffraction_amplitude(wavelength, slit_width, diffraction_angle):
    """
    Calculates the diffraction amplitude for a 1D slit of width 'slit_width' at a diffraction angle 'diffraction_angle'
    for an incident wave with wavelength 'wavelength'.
    """
    sinc_arg = (np.pi * slit_width / wavelength) * np.sin(diffraction_angle)
    amplitude = (slit_width / (2 * np.sqrt(2 * np.pi))) * np.sinc(sinc_arg)
    return amplitude

def double_pinhole_interference(x, gamma, I1, I2, d, w1,w2 = None, x0 = None, wav = ekev2wav(9.3), z = 1, norm = 1):
    """ 
    Young's double pinholes experiment simulation
    Adapted from: https://github.com/ThomasWodzinski/coherenceanalysis/
    T Wodzinski et al 2020 J. Phys. Commun. 4 075014
    
    :param x: spatial coordinates
    :param x0: shift of intensity peak
    :param wav: photon wavelength
    :param z: source to screen distance
    :param w1: width of first slit
    :param w2: width of second slit
    :param I1: intensity through first slit
    :param I2: intensity through second slit
    :param x1: position of first slit
    :param x2: position of second slit
    :param gamma: mutual coherence
    :param norm: normalisation factor
    
    :returns ii: interference pattern 
    """
    
    
    x1 = -d/2
    x2 = d/2
    
    if w2 is None:
        w2 = w1
        
    k = 2 * math.pi / wav
    theta = -(k * (d * (x - x0) / z))
 
    I = (
        0
        + I1 * Airy((x - x0), w1, wav, z, x1) ** 2
        + I2 * Airy((x - x0), w2, wav, z, x2) ** 2
        + 2
        * np.sqrt(I1)
        * Airy((x - x0), w1, wav, z, x1)
        * np.sqrt(I2)
        * Airy((x - x0), w2, wav, z, x2)
        * gamma
        * np.cos(theta)
    )
    

    I_normalized = norm * I / np.max(I)

    return I_normalized


def fix_experimental_params(wav = ekev2wav(9.3), z = 3.5, norm = 1, x0 = 0):
    """
    bind experimental parameters using functools partial
    """
    fitting_function = partial(double_pinhole_interference,
                               z = z, 
                               wav = wav,
                               norm = norm,
                               x0 = x0)
    
    return fitting_function 


def fit_pinhole_interference(ii, x, f, p0 = (0.85, 1, 1, 25e-06, 7.5e-06), b0 = ([0,0.25,0.25, 5e-06, 2e-06],[1,1,1, 100e-06, 15e-06])):
    """
    fits intensity data of a double pinhole experiment described by the function f
    
    :param ii: one-dimensional double-pinhole intensity 
    :param x: intensity coordinates - 1d array the same size as ii
    :param f: pinhole fit function (python method)
    :param p0: initial guess of pinhole parameters, tuple of length equal to number of fit function argument
    :param b0: optimisation bounds of fit functions list containing two tuples of equal length to p0
    """
    p,_ = curve_fit(f, x, ii, p0=p0, maxfev = 10000, bounds = b0)
    
    y_pred = f(x, *p)
    r = r2_score(ii, y_pred)
    
    return p, r
