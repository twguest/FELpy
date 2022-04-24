#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FELPY

@author: twguest
Created on Tue Mar 15 22:11:57 2022

__version__ = "1.0.1"
__email__ = "twguest@students.latrobe.edu.au"
"""

from scipy.optimize import 

def fitExponent(acov, ax, sigmaGuess = 1):
    """
    optimisation for fit of autocovariance function
    
    :param acov: autocovariance function [1d array]
    
    :returns params: parameter fit of autocov function
    """
    

    params, params_covariance = optimize.curve_fit(gaussian_func, ax, acov[:len(acov)-1],
                                                   p0=[1,0,sigmaGuess*50,0], maxfev = 10000)
    return params


if __name__ == '__main__':
    pass
# -*- coding: utf-8 -*-

