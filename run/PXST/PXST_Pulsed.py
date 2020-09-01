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

import sys

sys.path.append("../../")
import os
 
from model.beamline.structure import BeamlineModel, propParams
from model.src.coherent import coherentSource
 
from wpg.beamline import Beamline
from wpg.srwlib import SRWLOptD as Drift

from wpg.wavefront import Wavefront
from model.h5_tools import obj, object2h5, h52object  
 
import numpy as np 
from utils.job_utils import JobScheduler
from model.materials.phaseMask import phaseMask
from scipy.ndimage import gaussian_filter
from utils.os_utils import mkdir_p

from wpg.srwl_uti_smp import srwl_opt_setup_transm_from_file as Sample

def propThruMask(wfr, _x = 0, _y = 0):
    """
    propagate through the speckle generator
    """
    s = Sample(filepath = "../../data/samples/AAO.png", 
               rx = 10e-09, ry = 10e-09,
               thickness = 60e-06,
               delta = 3.35e-03,
               atten_len = 20e-06,
               xc = 0, yc = 0,
               shift_x = _y, shift_y = _x)


    bl = Beamline()
    bl.append(s, propParams(1,1,1,1,mode = 'normal'))
    bl.propagate(wfr)
    
    return wfr

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
    
    
    bl = spb.get_beamline()
    bl.propagation_options[0]['optical_elements'][-1].L += foc2sample
    bl.propagation_options[0]['propagation_parameters'][-1] = propParams(1, 1, 1, 1, mode = 'quadratic')

    bl.propagate(wfr)
    
    return wfr


def propagate2detector(wfr, sample2detector):
  
    bl = Beamline()
    bl.append(Drift(sample2detector), propParams(2, 1, 1, 1, mode = 'quadratic'))
    bl.propagate(wfr)
    
    return wfr

def resizeImage(wfr):
    
    bl = Beamline()
    
    fov_x = 2560*6.5e-06
    fov_y = 2160*6.5e-06
    
    [xMin, xMax, yMin, yMax] = wfr.get_limits()
    bl.append(Drift(0), propParams(fov_x/(xMax-xMin), 1, fov_y/(yMax-yMin), 1, mode = 'normal'))
    bl.propagate(wfr)
    
    return wfr


def launch():
    """
    This part launches the jobs that run in main 
    """
    
    cwd = os.getcwd()
    script = "PXSE_experiment.py"
    
    js = JobScheduler(cwd + "/" + script, logDir = "../../logs/",
                      jobName = "coherentPXST", partition = 'exfel', nodes = 4, jobType = 'single')
    
        
    js.run(test = False)
    
if __name__ == '__main__':
    
 
 
    focus = "nano"
    glob = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/PXST/"
    outdir = glob + "pulsed/"
    indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB/in/"
    
    mkdir_p(glob)
    mkdir_p(outdir)
    
    wfr = Wavefront()
    wfr.load_hdf5(indir + "NanoKB-Pulse_12.h5")
    
    print("Wavefront Loaded")
    wfr = propagateToMask(wfr, 500e-06)

    
    propThruMask(wfr, _x = 0, _y = 0)
    
    sample2detector = 3.5
    wfr = propagate2detector(wfr, sample2detector)
    

    detectorIntensity = wfr.get_intensity().sum(axis = -1)
    np.save(outdir + "detectorIntensity", detectorIntensity)
    




