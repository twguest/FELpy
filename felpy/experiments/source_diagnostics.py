#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FELPY

@author: %(username)s
Created on %(date)s

__version__ = "1.0.1"
__email__ = "trey.guest@xfel.eu"
"""

import numpy as np
 
from felpy.model.source import SA1_Source
from felpy.utils.vis_utils import Grids

from wpg.wpg_uti_wf import plot_intensity_map
from tqdm import tqdm
from felpy.analysis.statistics.correlation import norm
from matplotlib import pyplot as plt
from wpg.wpg_uti_wf import calc_pulse_energy

from felpy.utils.np_utils import memory_map

def get_pulse_energy(wfr):
    wfr.set_electric_field_representation('t')
    energy = calc_pulse_energy(wfr)[0]
    wfr.set_electric_field_representation('f')
    return energy

def plot_average_spectra(src, n_spectra = 5000, ax = None):
    """
    defines and fits the average spectra and temporal profiles of a given statistically defined source
    
    we expect that the source (src) has a temporal profile which can be called by get_temporal profile
    
    :param src: felpy.model.src.Source type object
    :param n_spectra: number of spectra over which we take the average
    """ 
    
    temporal_profiles = np.zeros([src.nz, n_spectra])
    freq_profiles = np.zeros([src.nz, n_spectra])

    tMin, tMax = -src.pulse_duration/2, src.pulse_duration/2
    fMin, fMax = 1/tMin, 1/tMax

    t = np.linspace(tMin, tMax, src.nz)*1e15
    
    for itr in tqdm(range(n_spectra)):
        temporal_profiles[:, itr] = src.get_temporal_profile(refresh = True, sigma = 3)

    
    #for i in range(10):
        #plt.plot(abs(temporal_profiles[:, i])**2)
   
       
    ### plotting     
    if ax is None:
        grid = Grids(global_aspect = 1.5)
        grid.create_grid(n = 1, m = 1, sharex = False, sharey = True)

    ax = grid.axes
    
    
    for i in range(n_spectra):
        if i % (n_spectra/10)  == 0:
            ax.plot(t, norm(abs(temporal_profiles[:,i])**2), alpha = 0.25, linestyle = 'dashed')
    

    ax.plot(t, norm((abs(temporal_profiles)**2).mean(-1)), color = 'k', linestyle = 'solid', label = "Average Temporal Profile")
    
    fontsize = 16
    ax.set_xlabel("Time (fs)", fontsize = fontsize)
    ax.set_ylabel("Normalised Intensity (a.u.)", fontsize = fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.xaxis.label.set_size(fontsize)
    ax.yaxis.label.set_size(fontsize)
    

def scan_source_energy(ekev, q):
    """ 
    objective: generate a function that records the pulse energy in J for a
    given set of input parameters - here we use ekev and q to align w/ SA1 req.
    
    :param ekev: photon energy (keV) or list of photon energies
    :param q: electron bunch charge (nC) or list of bunch charges
    
    :returns data: array of pulse-energies - len(ekev)xlen(q)
    usage
    >>> data = scan_source_energy(ekev = np.arange(1,5), q = [0.1])
    
    this is part of the source validation criterion.
    """

    
    data = np.zeros([len(ekev), len(q)])
    
    for i, e in tqdm(enumerate(ekev)):
        for k, Q in enumerate(q):
            
            src = SA1_Source(ekev = e, q = Q, nx = 512, ny = 512)

            data[i,k] = get_pulse_energy(src.wfr)
            
    return data

def scan_source_size(ekev, q = 0.25):
    """ 
    objective: generate a function that records beam area in [um] for a given
    set of input parameters - here we use ekev to align w/ empirical def of SA1
    
    :param ekev: list/array of photon energies (keV)
    :param q: arbitrary choice of electron beam charge (nC)
    
    :returns data: list of beam-sizes - len = len(ekev)
    
    usage
    >>> data = scan_source_size(ekev = np.arange(1,5))
    """
    
    data = np.zeros([len(ekev)])
    
    
    for i, e in tqdm(enumerate(ekev)):
        src = SA1_Source(ekev = e, q = q, nx = 1024, ny = 1024, S = 1)
        data[i] = src.get_fwhm()[0][0]
    
    return data

def scan_source_divergence(ekev, q, n = 10):
    
    data = np.zeros([len(ekev), 2, len(q), n])
    
    for j in tqdm(range(n)):
        for i, e in enumerate(ekev):
            for k, Q in enumerate(q):
                    
                src = SA1_Source(ekev = e, q = Q, nx = 512, ny = 512,
                                 xMin = -300e-06, xMax = 300e-06, yMin= -300e-06, yMax= 300e-06, S = 1)
                
                data[i,:,k,j] = src.get_divergence()[0]
   
    return data    



if __name__ == '__main__':
# =============================================================================
#     
#     src = SA1_Source(ekev = 5.2, q = 0.25, nx = 512, ny = 512,
#                      xMin = -400e-06, xMax = 400e-06, yMin= -400e-06, yMax= 400e-06)
#     
# =============================================================================
    #plot_average_spectra(src, n_spectra = 100000)

    #wfr = src.get_wfr()
    
    #data = scan_source_energy(ekev = np.arange(1,5), q = [0.1])
    #data = scan_source_size(ekev = np.arange(1,5))
    data = scan_source_divergence(ekev = np.arange(5,10), q = [0.1], n = 1)
    
    #print(get_pulse_energy(wfr))
    #plot_intensity_map(wfr)    