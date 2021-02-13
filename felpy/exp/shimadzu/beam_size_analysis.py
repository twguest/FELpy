#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:32:08 2021

@author: twguest
"""


import numpy as np
from matplotlib import pyplot as plt
from felpy.utils.vis_utils import plot_fill_between


DATA_DIR = "/media/twguest/shared/data/nkb_sensing/beam_area_r0046.npy"

def load_data(DATA_DIR):
    
    beam_area = np.load(DATA_DIR)
    return beam_area

def beam_size_analysis(data):
    
    ### load beam area
    plot_fill_between(data,
                      title = "r0046 Beam Area",
                      xlabel = "Pulse No.",
                      ylabel = "Beam Area ($\mu$m)",
                      plot_max = False)
      
    ### plot beam average and variance
    plt.plot([np.mean(np.gradient(data[:,itr])) for itr in range(data.shape[-1])])
    plt.show()
    
if __name__ == '__main__':
    data = load_data(DATA_DIR)[3:103,]*1e6
    data = np.delete(data, 24, axis = -1)
    beam_size_analysis(data)