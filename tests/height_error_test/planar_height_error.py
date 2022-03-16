#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:10:46 2020

@author: twguest
"""



###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################

import time
from copy import copy
import numpy as np
from matplotlib import pyplot as plt
from felpy.model.source.coherent import construct_SA1_wavefront
from felpy.model.instrument import Instrument
from felpy.model.materials.load_refl import load_refl, get_refl
import seaborn as sns
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from scipy.ndimage import gaussian_filter

from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.utils.vis_utils import colorbar_plot
from felpy.utils.np_utils import get_mesh
from felpy.model.tools import propagation_parameters

def get_optical_path_difference(bl):

    mir = bl.propagation_options[0]['optical_elements'][0]
    opd = np.array(mir.get_data(3))
    opd = np.reshape(opd, [1580*1581])
    
    return opd

if __name__ == '__main__':
    ekev = 5.0
    spb = get_beamline_object(ekev = ekev, crop = ["HOM1", "HOM1"])
 
    
    ### plane-wave
    wfr = construct_SA1_wavefront(512, 512, ekev, .20)
    wfr.data.arrEhor[128:384,128:384,:,0] = 1
    wfr.data.arrEhor[np.where(wfr.data.arrEhor[:,:,:,0] != 1)] = 0
    wfr.data.arrEhor[:,:,:,1] = 0j
    
    ### plot
    ii = wfr.get_intensity().sum(-1)
    ph = wfr.get_phase().sum(-1)
    mesh = get_mesh(ii, *wfr.get_spatial_resolution())
    
    colorbar_plot(ii, mesh = mesh,
                  xlabel = "x (mm)",
                  ylabel = "y (mm)",
                  scale = 1e3,
                  context = 'talk',
                  clabel = "Intensity (a.u.)",
                  grid = False)

    colorbar_plot(ph, mesh = mesh,
                  xlabel = "x (mm)",
                  ylabel = "y (mm)",
                  scale = 1e3,
                  context = 'talk',
                  clabel = "Phase",
                  vmin = 0,
                  vmax = 2*np.pi,
                  cmap = 'hsv',
                  grid = False)
    
    
    ### optical path difference of first mirror
    opd = get_optical_path_difference(spb)
    
    colorbar_plot(opd, mesh = mesh,
                  xlabel = "x (mm)",
                  ylabel = "y (mm)",
                  scale = 1e3,
                  context = 'talk',
                  clabel = "Optical Path Difference",
                  cmap = 'jet',
                  grid = False)
    
    
    spb.propagate(wfr, propagation_parameters(1,1,1,1,mode = 'fresnel'))
    
    ph2 = wfr.get_phase().sum(-1)
    colorbar_plot(ph2, mesh = mesh,
                  xlabel = "x (mm)",
                  ylabel = "y (mm)",
                  scale = 1e3,
                  context = 'talk',
                  clabel = "Phase",
                  vmin = 0,
                  vmax = 2*np.pi,
                  cmap = 'hsv',
                  grid = False)
# =============================================================================
#     
#     bl = setupBL(ekev, toggle = 'on')
#     bl.propagate(wfr)
#     
#     plotIntensity(wfr)
#     
#     on = get_optical_path_difference(bl)
#     wi = wfr.data.arrEhor
#     wfr = construct_SA1_wavefront(1024, 1024, ekev, .25)
#     
#     bl = setupBL(ekev, toggle = 'off')
#     bl.propagate(wfr)
#     wo = wfr.data.arrEhor
#     plotIntensity(wfr)
#     
#     off = get_optical_path_difference(bl)
#     
# =============================================================================
