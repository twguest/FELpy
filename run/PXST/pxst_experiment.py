#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 05:52:37 2020

@author: guestt
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""2
Created on Tue Jun 30 15:09:54 2020

Python Operations for NanoKB test w/ FAST PULSE
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

import multiprocessing

from model.beamline.structure import BeamlineModel, propParams
from model.src.coherent import coherentSource
from wpg import srwlib

from wpg.wavefront import Wavefront
from wpg.beamline import Beamline
from wpg.srwlib import SRWLOptD as Drift


from os import listdir
 
from model.tools import loadWavefront
from wpg.multisliceOptE import greyscaleToSlices

import numpy as np 

def propThruMask(wfr, _x = 0, _y = 0):
    
    slices = greyscaleToSlices("../../data/samples/AAO.png",
                               propagationParameters = propParams(1,1,1,1, mode = 'normal'),
                               slices = 1,
                               resolution = [10e-09, 10e-09,  None],
                               thickness = 60e-06,
                               delta = 3.35e-03,
                               attenLength = 20e-06,
                               rotate_angle = 0,
                               shift_x = _x,
                               shift_y = _y)
    
    slices.propagate(wfr)

def propagateToMask(wfr, foc2sample):
    
    """
    propagate from the focal plane to the mask
    
    :param wfr: wpg wavefront structure
    :param foc2sample: focus to sample distance [m] (float)
    """
    
    spb = BeamlineModel()
    
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "on", aperture = True, overwrite = False)
    
    spb.buildElements(focus)
    spb.buildBeamline(focus)
    print("Setup Beamline")
    
    spb.scale(wfr) 
    print("Source Scaled")
    bl = spb.get_beamline()
    bl.propagation_options[0]['optical_elements'][-1].L += foc2sample
    bl.propagation_options[0]['propagation_parameters'][-1] = propParams(1, 1, 1, 1, mode = 'quadratic')

    bl.propagate(wfr)
    
    return wfr

 




def propagate2detector(wfr, sample2detector):
    
    
    nx = 2560
    ny = 2160
    px = 6.5e-06
    
    fov_x = nx*px
    fov_y = ny*px
    
    
    [xMin, xMax, yMin, yMax] = wfr.get_limits()
    inx = wfr.params.Mesh.nx
    iny = wfr.params.Mesh.ny
    
    
    bl = Beamline()
    bl.append(Drift(sample2detector), propParams(1, 1, 1, 1, mode = 'quadratic'))
    bl.propagate(wfr)
    

    plotIntensity(wfr)
    

    wfr.save_tif("inte")
    return wfr

def resizeImage(wfr):
    
    bl = Beamline()
    
    fov_x = 2560*6.5e-06
    fov_y = 2560*6.5e-06
    
    [xMin, xMax, yMin, yMax] = wfr.get_limits()
    bl.append(Drift(0), propParams(fov_x/(xMax-xMin), 1, fov_y/(yMax-yMin), 1))
    bl.propagate(wfr)

if __name__ == '__main__':
    
    from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
    
    focus = "nano"
    indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/out/"
    gaussDir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/legacy/legacyOut.h5"


    fname = "NanoKB-Pulse_12.h5"
    print("loading wfr")
    #wfr = loadWavefront(gaussDir)
    
    wfr = coherentSource(2560, 2160, 4.96, 0.25)
    plotIntensity(wfr)
    print("Wavefront Loaded")
    wfr = propagateToMask(wfr, 500e-06)
    plotIntensity(wfr)
    propThruMask(wfr)
    plotIntensity(wfr)
    wfr = propagate2detector(wfr, 3.5)
    plotIntensity(wfr)
    resizeImage(wfr)
    plotIntensity(wfr)