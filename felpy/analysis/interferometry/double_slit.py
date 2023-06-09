import numpy as np

from numpy import sin,cos
from scipy.optimize import curve_fit
from felpy.utils.opt_utils import ekev2wav
from functools import partial
from sklearn.metrics import r2_score

def envelope(x, sw, wav, z):
    a = (np.pi*sw*x)/(wav*z)
    A = (sin(a)/a)**2

    return A

def double_slit_interference(x, gamma, I1, I2, d, w1,x0, w2 = None, wav = ekev2wav(9.3), z = 1, norm = 1):
    """ 
    regular double slit pattern across a screen detector a distance z (m) upstream
    
    :param x: array of horizontal coordinates (m)
    :param d: slit_seperation
    :param w: slit width
    :param wav: wavelength
    :param p: ratio of intensity between slits
    :param z: distance between slit and detector (m)
    :param mu: absolute mutual coherence b/w field between slits
    """
    
    
    p = I1/I2
    
    px = abs(x[1]-x[0])
    b = (np.pi*px*d)/(wav*z)

    B = sin(b)/b

    c = cos((2*np.pi*d*x)/(wav*z))
    C = B*c

    A = envelope(x,w1,wav,z)
    V = (2*abs(gamma)*np.sqrt(p))/(1+p)
    ii = A*(1+abs(V)*C)
    ii /= np.max(ii)
    
    return ii


def fix_experimental_params(wav = ekev2wav(9.3), z = 3.5, w = 7.5e-06):
    """
    bind experimental parameters using functools partial
    """
    fitting_function = partial(double_slit_interference,
                               z = z, 
                               wav = wav,
                               w=w)
    
    return fitting_function 
    

def fit_slit_interference(ii, x, f, p0 = (0.85, 1, 25e-06), b0 = ([0,0.25,1e-06],[1,1,250e-06])):
    """
    """
    p,_ = curve_fit(f, x, ii, p0=p0, maxfev = 10000, bounds = b0)
    
    y_pred = f(x, *p)
    r = r2_score(ii, y_pred)
    
    return p, r
