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

from model.src.coherent import coherentSource
from model.beamline.structure import BeamlineModel
from model.materials.load_refl import load_refl, get_refl
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from wpg.wpg_uti_wf import plot_intensity_qmap as plotPhase


def testMicron(ekev, q):
    
    wfr = coherentSource(1024, 1024, ekev, q)
    
    spb = BeamlineModel(overwrite_mirrors =  True)
    
    refl_data = load_refl()
    refl, ang = get_refl(refl_data, ekev, ang = 2.2e-03, limits = [1.1e-03, 3.6e-03])
    
    spb.adjustHOMs(refl, ang)
    
    spb.buildElements(focus = "micron")
    spb.buildBeamline(focus = "micron")
    bl = spb.get_beamline()
    
    bl.propagateSeq(wfr, "../../tmp/")

def testNano(ekev, q):
    
    wfr = coherentSource(1024, 1024, ekev, q)
    
    spb = BeamlineModel(overwrite_mirrors =  True)
    
    refl_data = load_refl()
    refl, ang = get_refl(refl_data, ekev, ang = 2.2e-03, limits = [1.1e-03, 3.6e-03])
    
    spb.adjustHOMs(refl, ang)
    
    spb.buildElements(focus = "nano")
    spb.buildBeamline(focus = "nano")
    bl = spb.get_beamline()
    
    bl.propagateSeq(wfr, "../../tmp/")

if __name__ == '__main__':
    
    ekev = 12.0
    q = 0.25
    
    #testMicron(ekev, q)
    testNano(ekev,1)

    