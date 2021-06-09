#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 12:35:15 2020

@author: twguest
"""


import time
import os
import shutil
import numpy as np
from matplotlib.colors import LogNorm
import wpg.srwlib as srwlib
#from wpg.wpg_uti_wf import calc_pulse_energy, get_axial_power_density, get_centroid
from matplotlib import pyplot as plt
from matplotlib import ticker, cm
from felpy.model.tools import argmax2d
from wpg import srwlib
from tqdm import tqdm
from felpy.model.tools import create_circular_mask
from felpy.analysis.optics.scalar.centroid import get_com
#from wpg.wpg_uti_wf import get_axis

from scipy.stats import norm

fwhm_area = norm.cdf(np.sqrt(2*np.log(2))) - norm.cdf(-np.sqrt(2*np.log(2)))

def finder(ii, nx, ny, c, efraction = fwhm_area, threshold = .01, iguess = None, nguesses = 10, fname = None, outdir = None, VERBOSE = True):
    """
    Modified Bisection Algortihm, iterates over a set of guesses to reduce the 
    error function 1-(0.95/(ii-ii*mask(r))) in searching for a solution for r
    
    :param wfr: wpg wavefront structure
    :param efraction: desired energy fraction for enclosed energy, 0-1 [float]
    :param threshold: desired accuracy of algorithm [float] % of efraction
    :param iguess: initial radius guess [int]
    :param nguesses: number of initial guess (only activate when iguess = None) [int]
    
    :returns r: radius in pixels
    :returns error: error of radius
    """

    i = 0
    llim = efraction - threshold*efraction
    ulim = efraction + threshold*efraction 
    
    ### make some random guesses
    if iguess is None:
        ins = np.random.randint(0, nx, nguesses)
        h = []
        for gs in ins:
            mask = create_circular_mask(nx, ny, c, r = gs)
            arr = mask*ii
            h.append(arr.sum()/ii.sum())
        
        gs = ins[np.argmin(h)]
        
        va = min(h)
        del h
    else:
        gs = iguess
        mask = create_circular_mask(nx, ny, c, r = gs)
        arr = mask*ii
        va = arr.sum()/ii.sum()[nx,ny]
        
        

    ### now do some bisection style loop
    
    if VERBOSE:
        print("\nInitial Guess: ", gs, "Initial Error:", va)
    
    
    while (ulim  <= va or va <= llim) and i < nx:
        
        if  va <= llim:
            gs = gs + 1#gs//2
        elif va >= ulim:
            gs = gs - 1#gs//2
            
        
        mask = create_circular_mask(nx, ny, c, r = gs)
        arr = mask*ii
        va = arr.sum()/ii.sum()

        i += 1

    r = gs
    err = abs(efraction-va)
    
    if VERBOSE:
        print("Enclosed Energy Radius: {}".format(r), va)
        print("Error = {:.4f}".format(err))
        print("After {} Iterations".format(i))
        
    
    return r, err


    
def plotEnclosed(ii, r, c, label = None, outdir = None):
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot()
    
    ax1.imshow(ii, cmap = 'jet')
    circle = plt.Circle(c, r, color='w', fill=False)
    ax1.add_artist(circle)
    plt.axis("off")
    ax1.plot(c[0],c[1], marker = 'x', color = 'k')
    if label is not None:
        ax1.text(15, ii.shape[1]-15, label, c = 'w')
        


def get_enclosed_energy(ii, dx, dy, efraction = fwhm_area, sdir = None, plot = False, VERBOSE = True, threshold = 0.01):
    
    nx, ny = ii.shape
    c = get_com(ii)
    results, err = finder(ii, nx, ny, c, efraction, VERBOSE = VERBOSE, threshold = threshold)
    
    
    if plot:

        plotEnclosed(ii, results, c)
    
    results = (float(results)*dx, float(results*dy))
    
    return results, err
    
 

if __name__ == '__main__':
    from scipy.ndimage import gaussian_filter
    ii = np.zeros([50,50])
    ii[20,20] = 100
    ii[30,30] = 100
    #ii = gaussian_filter(ii, 5)
    #ii[40,40] = 0.0001
    ii = gaussian_filter(ii, 6)
    dx = 1e-06
    dy = 1e-06
    
    results = get_enclosed_energy(ii, dx, dy, efraction = 0.1, plot = True, 
                                  VERBOSE = True)
