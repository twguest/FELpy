#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FELPY

@author: twguest
Created on Tue Mar 15 22:11:57 2022

__version__ = "1.0.1"
__email__ = "twguest@students.latrobe.edu.au"
"""
import numpy as np

from scipy import optimize
from felpy.math.fit_funcs import gaussian as gaussian_func

def fitExponent(acov, ax, p0 = None):
    """
    optimisation for fit of autocovariance function
    
    :param acov: autocovariance function [1d array]
    
    :returns params: parameter fit of autocov function
    """
    
    if p0 is None:
        p0 = [np.max(acov),0,np.std(acov),0]

    params, params_covariance = optimize.curve_fit(gaussian_func, ax, acov[:len(acov)],
                                                   p0=p0, maxfev = 10000)
    return params


if __name__ == '__main__':
    pass
# -*- coding: utf-8 -*-

