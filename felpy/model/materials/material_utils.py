#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.0"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""

import numpy as np


def add_extent(arr, extent):
    """
    add extent to phaseshift array
    
    :param arr: phaseshift (or otherwise) array
    :param extent: extent to be added (following srw_opt_setup_surf_height_2d)
    """
    
    x = np.ones([1, arr.shape[1]])
    x[0,:] = np.linspace(-extent[0]/2, extent[0]/2, arr.shape[0])
    
    y = np.ones([arr.shape[0]+1, 1])
    
    y[1:,0] = np.linspace(-extent[1]/2, extent[1]/2, arr.shape[0])
    
    arr = np.vstack([x,arr])
    arr = np.hstack([y,arr])
    
    arr[0,0] = 0
    return arr    

if __name__ == "__main__":
    arr = np.ones((50,50))
    ext = [5e-03,5e-03]
    
    arr = add_extent(arr, ext)