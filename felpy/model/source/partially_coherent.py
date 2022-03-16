#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 14:21:55 2021

@author: twguest
"""

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from felpy.utils.maths.generator_funcs import complex_gaussian_2d
 
kx = np.pi/4
ky = 0

wav = 0.1e-10
k = np.pi*2/wav

z = 1
x = -np.linspace(-400, 400, 400)
y = np.linspace(-400, 400, 400)
grid = np.meshgrid(x,y)

def define_wfr_tilt(mesh, kx, ky):
    """ 
    define a two-dimensional complex tilted plane-wave, tilted along transverse pointing vectors kx and ky
    
    :param mesh: [2,nx,ny] array of coordinates, first-axis corresponds to x and y axis,
    see felpy.utils.np_utils.get_mesh for more
    :param kx: horizontal pointing vector [-pi, pi]
    :param ky: vertical pointing vector [-pi, pi]
    """
    wfr_tilt =   np.exp(1j*(kx*grid[0]+ky*grid[1])) #*np.exp(1j*k*z)
    return wfr_tilt

 


# -*- coding: utf-8 -*-

