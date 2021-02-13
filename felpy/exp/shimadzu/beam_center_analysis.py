#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:32:08 2021

@author: twguest
"""


import numpy as np
from matplotlib import pyplot as plt
from felpy.utils.vis_utils import plot_fill_between


DATA_DIR = "/media/twguest/shared/data/nkb_sensing/beam_center_r0046.npy"

def load_data(DATA_DIR):
    
    beam_area = np.load(DATA_DIR)
    return beam_area

def beam_center_analysis(data, color = 'blue'):
    
    ### load beam area
    plot_fill_between(data[:,:,0] - data[:,:,0].mean(),
                      title = "r0046 Horizontal Beam Center",
                      xlabel = "Pulse No.",
                      ylabel = "Horizontal Deviation from Mean ($\mu$m)",
                      plot_max = False,
                      color = color,
                      xlim = [0,data.shape[-3]])


    plot_fill_between(data[:,:,1] - data[:,:,1].mean(),
                      title = "r0046 Vertical Beam Center",
                      xlabel = "Pulse No.",
                      ylabel = "Vertical Deviationfrom Mean ($\mu$m)",
                      plot_max = False,
                      color = color,
                      xlim = [0,data.shape[-3]])

    plot_fill_between(np.moveaxis(data, 1,0)[:,:,0] - np.moveaxis(data, 1,0)[:,:,0].mean(),
                      title = "r0046 Horizontal Beam Center",
                      xlabel = "Pulse No.",
                      ylabel = "Horizontal Deviation from Mean ($\mu$m)",
                      plot_max = False,
                      color = color,
                      xlim = [0,np.moveaxis(data, 1,0).shape[-3]])

    plt.plot(np.arange(data.std(0)[:,0].shape[0]),data.std(0)[:,0])
    
if __name__ == '__main__':
    data = load_data(DATA_DIR)[3:103,]*1e6*28.9e-06
    data = np.delete(data, 24, axis = -2)
    beam_center_analysis(data, color = 'red')