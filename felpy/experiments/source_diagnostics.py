#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FELPY

@author: %(username)s
Created on %(date)s

__version__ = "1.0.1"
__email__ = "twguest@students.latrobe.edu.au"
"""

import numpy as np
 
from felpy.model.core.source import SA1_Source
from felpy.utils.vis_utils import Grids

from wpg.wpg_uti_wf import plot_intensity_map
from tqdm import tqdm
from felpy.analysis.statistics.correlation import norm
from matplotlib import pyplot as plt
from wpg.wpg_uti_wf import calc_pulse_energy


def plot_average_spectra(src, n_spectra = 5000):
    """
    defines and fits the average spectra and temporal profiles of a given statistically defined source
    
    we expect that the source (src) has a temporal profile which can be called by get_temporal profile
    
    :param src: felpy.model.core.source.Source type object
    :param n_spectra: number of spectra over which we take the average
    """
    
    temporal_profiles = np.zeros([src.nz, n_spectra])
    freq_profiles = np.zeros([src.nz, n_spectra])

    tMin, tMax = -src.pulse_duration/2, src.pulse_duration/2
    fMin, fMax = 1/tMin, 1/tMax

    t = np.linspace(tMin, tMax, src.nz)*1e15
    
    for itr in tqdm(range(n_spectra)):
        temporal_profiles[:, itr] = src.get_temporal_profile(refresh = True, sigma = 3, S = 2)

    
    #for i in range(10):
        #plt.plot(abs(temporal_profiles[:, i])**2)
   
    plt.plot(abs(temporal_profiles[:, :10])**2)
        
    ### plotting     
    grid = Grids(global_aspect = 1.5)
    grid.create_grid(n = 1, m = 1, sharex = False, sharey = True)
    
    ax1 = grid.axes
    
    
    for i in range(n_spectra):
        if i % (n_spectra/10)  == 0:
            ax1.plot(t, norm(abs(temporal_profiles[:,i])**2), alpha = 0.25, linestyle = 'dashed')
    

    ax1.plot(t, norm((abs(temporal_profiles)**2).mean(-1)), color = 'k', linestyle = 'solid', label = "Average Temporal Profile")
    
    fontsize = 16
    ax1.set_xlabel("Time (fs)", fontsize = fontsize)
    ax1.set_ylabel("Normalised Intensity (a.u.)", fontsize = fontsize)
    ax1.tick_params(axis='both', which='major', labelsize=fontsize)
    ax1.xaxis.label.set_size(fontsize)
    ax1.yaxis.label.set_size(fontsize)
    
    
if __name__ == '__main__':
    
    src = SA1_Source(ekev = 5.2, q = 0.25, nx = 512, ny = 512,
                     xMin = -400e-06, xMax = 400e-06, yMin= -400e-06, yMax= 400e-06)
    
    #plot_average_spectra(src, n_spectra = 100000)

    wfr = src.get_wfr()
    
    #plot_intensity_map(wfr)    