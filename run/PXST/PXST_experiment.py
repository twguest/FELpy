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

 
from model.beamline.structure import BeamlineModel, propParams
from model.src.coherent import coherentSource
 
from wpg.beamline import Beamline
from wpg.srwlib import SRWLOptD as Drift


from model.h5_tools import obj, object2h5, h52object  
from wpg.multisliceOptE import greyscaleToSlices

import numpy as np 

from model.materials.phaseMask import phaseMask
from scipy.ndimage import gaussian_filter

def phaseScreen(wfr):
    
    
    wav = 2e-10
    ps = np.random.rand(250, 250)*wav
    ps *= gaussian_filter(ps, sigma = 75)
    
    opd = phaseMask(ps, [20e-06, 20e-06], wav)
    
    bl = Beamline()
    bl.append(opd, propParams(1,1,1,1, mode = 'normal'))
    bl.propagate(wfr)
    
    return wfr
    

def propThruMask(wfr, _x = 0, _y = 0):
    
    slices = greyscaleToSlices("../../data/samples/AAO.png",
                               propagationParameters = propParams(1,1,1,1, mode = 'normal'),
                               slices = 1,
                               resolution = [10e-09, 10e-09,  None],
                               thickness = 60e-06,
                               delta = 3.35e-03, #### via Henke
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
    
    
  
    bl = Beamline()
    bl.append(Drift(sample2detector), propParams(1, 1, 1, 1, mode = 'quadratic'))
    bl.propagate(wfr)
    
    return wfr

def resizeImage(wfr):
    
    bl = Beamline()
    
    fov_x = 2560*6.5e-06
    fov_y = 2560*6.5e-06
    
    [xMin, xMax, yMin, yMax] = wfr.get_limits()
    bl.append(Drift(0), propParams(fov_x/(xMax-xMin), 1, fov_y/(yMax-yMin), 1))
    bl.propagate(wfr)

if __name__ == '__main__':
    
    import imageio
    
    seed = sys.argv[0] ##random seed
    
    scanRange = 10e-06
    
    [pos_x, pos_y] = np.random.uniform(-scanRange//2, scanRange//2,2)
    
    data = obj()
    data.translation = [pos_x, pos_y]
    
    focus = "nano"
    outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/PXST/out/"
    
    wfr = coherentSource(1000, 1000, 4.96, 0.25)
    
    print("Wavefront Loaded")
    wfr = propagateToMask(wfr, 500e-06)
    
    wfr = phaseScreen(wfr)
    
    data.entryWavefront = wfr.toComplex()[0,:,:,0]
    data.x_pixel_size_i = wfr.pixelsize()[0]
    data.y_pixel_size_i = wfr.pixelsize()[1]
    
    
    propThruMask(wfr, _x = pos_x, _y = pos_y)
    
    
    sample2detector = 3.5
    wfr = propagate2detector(wfr, sample2detector)
    
    resizeImage(wfr)
    
    data.ComplexSolution = wfr.toComplex()[0,:,:,0]
    data.detectorIntensity = wfr.get_intensity().sum(axis = -1)
    data.x_pixel_size = wfr.get_pixel_size()[0]
    data.y_pixel_size = wfr.get_pixel_size()[1]
    data.distance = sample2detector
    data.mask(imageio.imread("../../data/samples/AAO.png"))
    data.wavelength = wfr.params.wavelength
