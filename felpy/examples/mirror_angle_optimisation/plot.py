#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 14:12:53 2021

@author: twguest
"""


import numpy as np

from os import listdir
from matplotlib import pyplot as plt

def get_max_np(file_loc):
    
    file = np.load(file_loc)
    return file

def get_max_all(test_dir):
    pass     

if __name__ == '__main__':
   file =  get_max_np("/opt/FELpy/felpy/tests/mirror_angle_optimisation/3.0keV/HOM1_mirror_flux_data.npy")
    
   plt.plot(file[0:,0],file[0:,1])