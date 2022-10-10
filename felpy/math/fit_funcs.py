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

def gaussian(x, a, x0, sigma,c):
    """
    gaussian function for fitting
    """
    
    return a * np.exp(-(x-x0)**2/(2*sigma**2)) + c


def parabola(t, a, t0, b):
    return -a*(t-t0)**2+b


if __name__ == '__main__':
    pass

