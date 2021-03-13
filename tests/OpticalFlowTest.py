#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 14:39:58 2020

@author: twguest


Let's test the geometric flow theory
"""

import sys

sys.path.append("/opt/spytlab")
sys.path.append("/opt/WPG/")
sys.path.append("/opt/spb_model")



import numpy as np

from model.materials.phaseMask import phaseMask
from model.beamline.structure import propagation_parameters
from model.tools import constructPulse
from utils.banded_utils import diagonal_form, solve_banded
from wpg.optical_elements import Drift
from felpy.model.core.beamline import Beamline
from OpticalFlow import processOneProjection
from felpy.model.core.wavefront import Wavefront
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from sklearn.preprocessing import minmax_scale as norm
from matplotlib import pyplot as plt
from model.src.coherent import construct_SA1_wavefront
from wpg import srwlib
from wpg.wpg_uti_wf import get_axis
from scipy.constants import h,c,e
from scipy.ndimage.filters import gaussian_filter
from matplotlib import colors

if __name__ == "__main__":
    
    slc = 2
    N = 2
    
    nx, ny = 512, 512
    II = np.zeros((nx,ny,N), 'float32')
    PH = np.zeros((nx,ny,N), 'float32')
    PHz = np.zeros((nx,ny,N), 'float32')
    A = np.zeros((N,N), 'float32')
    B = np.zeros((N,1), 'float32')
    
    

    
    
    wfr = construct_SA1_wavefront(nx,ny,4.96,0.25)
    wfr.store_hdf5("coherentSrc.h5")

    wav = (h*c)/(wfr.params.photonEnergy*e)
    
    pm = gaussian_filter(np.random.rand(5,5)/1000, sigma = 20)
    plt.imshow(pm)
    
    sp = phaseMask(pm,
                    [ get_axis(wfr, axis = 'x').max()-
                                            get_axis(wfr, axis = 'x').min(),
                                            get_axis(wfr, axis = 'y').max()-
                                            get_axis(wfr, axis = 'y').min()], wav) ##speckle

    

   
    
    

    bl = Beamline()
    bl.append(sp, propagation_parameters(1,1,1,1))
    bl.append(Drift(1), propagation_parameters(1,1,1,1, mode = 'normal'))
    bl.propagate(wfr)
    
    Ps = wfr.get_phase(slice_number = 0) #% np.pi*2
    Is =  wfr.get_intensity(slice_number = 0)
    
    
    #######################################
    
    
    wfr = Wavefront()
    wfr.load_hdf5("coherentSrc.h5")
    bl = Beamline()    
    bl.append(Drift(1), propagation_parameters(1,1,1,1, mode = 'normal'))
    bl.propagate(wfr)


    Pr = wfr.get_phase(slice_number = 0) #% np.pi*2
    Ir =  wfr.get_intensity(slice_number = 0)  


    for I in [Ir, Is]:    
        I[np.where(I == 0)] = 1e-10

    results = processOneProjection(Is, Ir)
    dy = results['dy']
    phi = results['phi'].real*wav
    
    phase = Ps-Pr
    plt.imshow(phase)
    plt.show()
    
    def plotNorm(phi):
        
        norm = colors.LogNorm(phi.mean() + 0.5 * phi.std(), phi.max(), clip='True')
        plt.imshow(phi, norm = norm)
        plt.show()
        
    plt.imshow(dy)