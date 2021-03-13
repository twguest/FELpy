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
from felpy.model.src.coherent import construct_SA1_wavefront
from felpy.model.core.instrument import Instrument
from felpy.model.materials.load_refl import load_refl, get_refl

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from scipy.ndimage import gaussian_filter

def setupBL(ekev, toggle = "on"):
    
    spb = Instrument()
    
    
    spb.build_elements(focus = "nano")
<<<<<<< HEAD
    spb.build_beamline(focus = "nano")
=======
    spb.buildBeamline(focus = "nano")
>>>>>>> 108cfb9b6fc97d3841ee1db54862523eee5b184e
    
    spb.crop_beamline(element1 = "HOM2")
    
    return spb.get_beamline()

def getOPD(bl):
    mir = bl.propagation_options[0]['optical_elements'][1]
    opd = np.array(mir.get_data(3))
    opd.reshape(1580*1581)
    
    return opd

if __name__ == '__main__':
    
    ekev = 9.2
    
    wfr = construct_SA1_wavefront(1024, 1024, ekev, .25)
    
    bl = setupBL(ekev, toggle = 'on')
    bl.propagate(wfr)
    
    plotIntensity(wfr)
    
    on = getOPD(bl)
    wi = wfr.data.arrEhor
    wfr = construct_SA1_wavefront(1024, 1024, ekev, .25)
    
    bl = setupBL(ekev, toggle = 'off')
    bl.propagate(wfr)
    wo = wfr.data.arrEhor
    plotIntensity(wfr)
    
    off = getOPD(bl)
    