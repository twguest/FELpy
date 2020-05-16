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
from model.src.coherent import coherentSource
from model.beamline.structure import BeamlineModel
from model.materials.load_refl import load_refl, get_refl

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 


def setupBL(ekev, toggle = "on"):
    
    spb = BeamlineModel()
    
    spb.setupHOMs(ekev, 2.2e-03)
    spb.mirrorProfiles(toggle = toggle, overwrite = True)
    
    spb.buildElements(focus = "micron")
    spb.buildBeamline(focus = "micron")
    
    spb.cropBeamline(element1 = "HOM2")
    
    return spb.get_beamline()

def getOPD(bl):
    mir = bl.propagation_options[0]['optical_elements'][1]
    opd = mir.get_data(3)
    opd.reshape([np.sqrt(opd), np.sqrt(opd)])
    
    return opd

if __name__ == '__main__':
    
    ekev = 9.2
    
    wfr = coherentSource(1024, 1024, ekev, .25)
    
    bl = setupBL(ekev, toggle = 'on')
    bl.propagate(wfr)
    
    plotIntensity(wfr)
    
    on = getOPD(bl)
    wi = wfr.data.arrEhor
    wfr = coherentSource(1024, 1024, ekev, .25)
    
    bl = setupBL(ekev, toggle = 'off')
    bl.propagate(wfr)
    wo = wfr.data.arrEhor
    plotIntensity(wfr)
    
    off = getOPD(bl)
    