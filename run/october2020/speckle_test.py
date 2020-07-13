#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 09:27:24 2020

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
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from model.src.coherent import coherentSource
from model.beamline.structure import BeamlineModel, propParams
from model.materials.load_refl import load_refl, get_refl
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from wpg.wpg_uti_wf import plot_intensity_qmap as plotPhase
from wpg.beamline import Beamline
from wpg.wavefront import Wavefront
from wpg.srwlib import SRWLOptD as Drift
from model.materials.phaseMask import phaseMask

from matplotlib import pyplot as plt
from wpg.multisliceOptE import greyscaleToSlices



drift_to_pam = Drift(3.0)
drift_to_pam.name = "drift_to_pam"


drift_to_ehc = Drift(0.5)
drift_to_ehc.name = "drift_to_ehc"

def propThruMask(wfr, rotate):
    
    slices = greyscaleToSlices("../../data/tests/speckle_test/speckle_crop.png",
                               propagationParameters = propParams(1,1,1,1),
                               slices = 1,
                               resolution = [10e-06/43, 10e-06/43,  None],
                               thickness = 250e-06,
                               delta = 1.82243693E-05,
                               attenLength = 1e3,
                               rotate_angle = rotate)
    
    slices.propagate(wfr)


def main():
    """
    task: to gneerate the output of a speckle tracking experiment
    """
    
    wfr = coherentSource(1024, 1024, 4.96, 0.250)
    
    spb = BeamlineModel()
    spb.setupHOMs(4.96, 2.2e-03)
    spb.setupKBs(4.96, 3.5e-03)
    spb.mirrorProfiles(toggle = "on", aperture = True, overwrite = True)
    spb.buildElements(focus = "nano")
    spb.buildBeamline(focus = "nano")
   
    bl = spb.get_beamline()
           
    bl.append(drift_to_pam, propParams(1,1,1,1, mode = 'quadratic'))
    bl.propagate(wfr)
    plotIntensity(wfr)
    
    #### THEN PROP THROUGH MASK
    propThruMask(wfr, rotate = 0)
    plotIntensity(wfr)
    
    ### then propagate to detector
    bl = Beamline()
    bl.append(drift_to_ehc, propParams(1/5, 2, 1/2, 2, mode = 'quadratic'))    
    
    bl.propagate(wfr)
    wfr.save_tif("../../tmp/speckle_test_toggleOn")
    plotIntensity(wfr)
   
    return wfr


if __name__ == '__main__':
    wfr = main()