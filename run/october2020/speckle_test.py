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

def propThruMask(wfr):
    
    slices = greyscaleToSlices("../../data/tests/sandpaper-1500nm-pixels.png",
                               propagationParameters = propParams(1,1,1,1),
                               slices = 1,
                               resolution = [0.5e-06, 0.5e-06,  None],
                               thickness = 250e-06,
                               delta = 1.82243693E-05,
                               attenLength = 1e3)
    
    slices.propagate(wfr)
    


def main():
    """
    task: to gneerate the output of a speckle tracking experiment
    """
    
    wfr = coherentSource(1024, 1024, 9.0, 0.25)
    
    extent = [1e-02, 1e-02]
    req_speckle_size = 10e-06
    npix = [int(extent[0]/req_speckle_size), int(extent[1]/req_speckle_size)]
    phaseshift = np.random.uniform(-10*np.pi, 10*np.pi, size = npix)
    
    speckle = phaseMask(phaseshift, extent, wfr.params.wavelength)

    
    spb = BeamlineModel()
    
    spb.setupHOMs(4.96, 2.2e-03)
    spb.setupKBs(4.96, 3.5e-03)
    spb.mirrorProfiles(toggle = "on", aperture = True, overwrite = True)
    spb.buildElements(focus = "nano")
    spb.buildBeamline(focus = "nano")
   
    bl = spb.get_beamline()
    
    
    drift_to_pam = Drift(3.0)
    drift_to_pam.name = "drift_to_pam"
    
    speckle.name = 'speckle_mask'
    
    drift_to_ehc = Drift(0.5)
    drift_to_ehc.name = "drift_to_ehc"
    
    
    #bl.append(drift_to_pam, propParams(1/2, 1, 1/2, 1, mode = 'quadratic'))
    bl.propagate(wfr)
    plotIntensity(wfr)
    propThruMask(wfr)
    
    bl = Beamline()
    bl.append(drift_to_ehc, propParams(1, 2, 1, 2, mode = 'quadratic'))    
    
    bl.propagate(wfr)
    wfr.save_tif("../../tmp/speckle_test")
    plotIntensity(wfr)
   
    return wfr


if __name__ == '__main__':
    wfr = main()