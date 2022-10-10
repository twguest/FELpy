#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.1"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""

import sys
import numpy as np
from numpy.lib.stride_tricks import as_strided
from scipy.ndimage import gaussian_filter
from matplotlib import pyplot as plt
from scipy import optimize


import matplotlib.colors as colors


def _check_arg(x, xname):
    """
    check if the argument x is one-dimensional
    
    :param x: array argument
    :param xname: identifier for x
    """
    
    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError('%s must be one-dimensional.' % xname)
    

def window_2D(arr, nrows, ncols):
    """
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size

    If arr is a 2D array, the returned array should look like n subblocks with
    each subblock preserving the "physical" layout of arr.
    """
    
    h, w = arr.shape
    assert h % nrows == 0, "{} rows is not evenly divisble by {}".format(h, nrows)
    assert w % ncols == 0, "{} cols is not evenly divisble by {}".format(w, ncols)
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def checkRoot(n):
    """
    function to check if value n has an integer root
    
    :param n: input value
    
    :returns BOOL: [bool]
    """    

    if np.sqrt(n) % 1 == 0 and np.sqrt(n)/2 % 1 == 0:
        BOOL = True
    else:
        BOOL = False
    return BOOL



def get_windows(arr, n):
    """
    useful wrapper for window - returns n windows, where sqrt(n) % 0 must be true
    
    :param arr: array to be windowed
    :param n: number of windows
    
    :returns w: 3D array of shape [n, nrows, ncols]
    """
    
    w = None
    
    BOOL = checkRoot(n)
    
    if BOOL:
        w = window_2D(arr, int(arr.shape[0]//np.sqrt(n)), int(arr.shape[1]//np.sqrt(n)))
    elif not BOOL:
        sys.exit("sqrt(n) is not an even integer value")
    
    if w is not None:
        return w
    
def odd_even(N):
    """
    check if a number is odd or even
    """
    
    if N % 2 == 0:
        ans = "even"
    elif N % 2 == 1:
        ans = 'odd'
    return ans
 

def window_2D(N, **kwargs):
    """ 
    generate a 2D window function
    
    :param N: window width - if no other keyword value submitted, then the array is assumed to be square
    
    :kwarg M: window height, if none, assumed to be equal to N
    :kwarg window: window method, e.g. numpy.hamming, scipy.signal.windows.hamm, etc, if none, assumed uniform.
    
    :returns: 2D numpy array of NxM dimensions
    """
    if 'M' in kwargs:
        M = kwargs['M']
    else:
        M = N
    if 'window' in kwargs:
        window = kwargs['window']
    else:
        window = np.ones
        
    return np.sqrt(np.outer(np.abs(window(N)),np.abs(window(M))))