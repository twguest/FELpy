#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:45:13 2020

A scratchpad for solving for beamline propagation


    
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

import numpy as np
import matplotlib as mpl
from copy import copy

from felpy.model.core.beamline import Beamline
from wpg.srwlib import SRWLOptD

from felpy.model.src.coherent import construct_SA1_wavefront
from felpy.model.beamlines.exfel_spb.exfel_spb import Instrument, propagation_parameters
from felpy.model.materials.load_refl import load_refl, get_refl

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from wpg.wpg_uti_wf import plot_intensity_qmap as plotPhase
from felpy.utils.vis_utils import colorbar_plot
from felpy.utils.np_utils import get_mesh
from matplotlib import pyplot as plt 

from felpy.model.beamlines.exfel_spb.methods import get_beamline_object

def sliceFocus(wfr, ekev, focus = 'nano', nslices = 500, axisName = 'x', outdir = None):
    
    spb = get_beamline_object(apertures = False, surface = False)
    
    el_n = len(spb.propagation_options[0]['optical_elements'])-1
    
    spb.propagation_options[0]['optical_elements'][el_n].L *= 0.98
    slice_interval = 2.2-copy(spb.propagation_options[0]['optical_elements'][el_n].L/nslices) 

    spb.propagation_options[0]['propagation_parameters'][el_n] = propagation_parameters(5,1,5,1, mode = 'fresnel')
    
    spb.propagate(wfr)
    plotIntensity(wfr)
    
    
    if axisName == 'x':
        data_focslice = np.zeros([nslices, np.shape(wfr.data.arrEhor)[0]])
    elif axisName == 'y':
        data_focslice = np.zeros([nslices, np.shape(wfr.data.arrEhor)[1]])
    
 
    
    
    for n in range(nslices):
        print("Slice {}/{}".format(n+1, nslices))
        
        bl = Beamline()
        bl.append(SRWLOptD(slice_interval), propagation_parameters(1,1,1,1, mode = 'fresnel'))
        bl.propagate(wfr)
        plotIntensity(wfr)
        data_focslice[-n-1, :] = wfr.get_intensity()[:,:,0].sum(-1)
        #plt.plot(data_focslice[-n-1, :])
        plt.show()
        
    
    y = np.linspace(spb.bl.propagation_options[0]['optical_elements'][el_n].L,
                    spb.bl.propagation_options[0]['optical_elements'][el_n].L + nslices*slice_interval,
                    nslices)
    
 

# =============================================================================
#     ax1 =     colorbar_plot(data_focslice,
#               get_mesh(data_focslice,
#                                wfr.get_spatial_resolution()[0]*1e6,
#                                y[1]-y[0]),
#               aspect = 'auto',
#               scale = 1,
#               #norm=mpl.colors.LogNorm(),
#               xlabel = "x($\mu m$)",
#               ylabel = "z(m)",
#               return_axes = True)
#     
# =============================================================================
    return data_focslice
    
 
# =============================================================================
#     
#     ax1.plot(data_focslice[0]*1e6, np.ones(data_focslice.shape[0])*spb.bl.params["df"]["distance"], 'r--')
#     
# =============================================================================

if __name__ == '__main__':
    
    ekev = 9.2
    
    wfr = construct_SA1_wavefront(512, 512, ekev, 0.5)
    
# =============================================================================
#     sliceFocus(wfr = wfr, ekev = ekev, focus = 'micron', nslices = 500, axisName = 'x', outdir = "focus_test/")
#     sliceFocus(wfr = wfr, ekev = ekev, focus = 'micron', nslices = 500, axisName = 'y', outdir = "focus_test/")
#     
# =============================================================================
    df = sliceFocus(wfr = wfr, ekev = ekev, focus = 'nano', nslices = 50, axisName = 'x', outdir = "focus_test/")
    #sliceFocus(wfr = wfr, ekev = ekev, focus = 'nano', nslices = 500, axisName = 'y', outdir = "focus_test/")