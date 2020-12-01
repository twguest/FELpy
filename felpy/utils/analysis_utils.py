# -*- coding: utf-8 -*-

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
    

def window(arr, nrows, ncols):
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
        w = window(arr, int(arr.shape[0]//np.sqrt(n)), int(arr.shape[1]//np.sqrt(n)))
    elif not BOOL:
        sys.exit("sqrt(n) is not an even integer value")
    
    if w is not None:
        return w