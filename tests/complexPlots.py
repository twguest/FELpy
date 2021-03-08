#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:54:43 2020

@author: twguest
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
 
def Complex2HSV(z, rmin, rmax, hue_start=90):
    # get amplidude of z and limit to [rmin, rmax]
    amp = np.abs(z)
    amp = np.where(amp < rmin, rmin, amp)
    amp = np.where(amp > rmax, rmax, amp)
    ph = np.angle(z, deg=1) + hue_start
    # HSV are values in range [0,1]
    h = (ph % 360) / 360
    s = 0.85 * np.ones_like(h)
    v = (amp -rmin) / (rmax - rmin)
    return hsv_to_rgb(np.dstack((h,s,v)))

plt.imshow(Complex2HSV(wfr.as_complex()[0,:,:,0], 100, np.max(wfr.as_complex()[0,:,:,0].real/1000)))