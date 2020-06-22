#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 11:44:18 2020

@author: twguest
"""


#################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################

from model.src.coherent import coherentSource
from model.beamline.structure import BeamlineModel
from model.beamline.samplingCalc import getImageProperties, getRequiredDistance

from wpg.wavefront import Wavefront

def get_wfr(element, outdir, focus = 'micron'):
    """
    Return the beamline at some optical element / transverse plane (coherent case)

    :param element: beamline element [str] to return beamline from
    :param outdir: directory for .h5 saving
    """  
    
    if focus == 'micron':
        ekev = 9.2
        q = 0.25
    elif focus == 'nano':
        ekev = 4.96
        q = 0.25
    
    wfr = coherentSource(1024, 1024, ekev, q)
    
    spb = BeamlineModel()
    
    spb.setupHOMs(ekev, 2.2e-03)
    spb.setupKBs(q, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "off", aperture = True, overwrite = True)
    spb.buildElements(focus = focus)
    spb.buildBeamline(focus = focus)
    spb.cropBeamline(element)
    
    bl = spb.get_beamline()
    bl.propagate(wfr)
    
    wfr.store_hdf5(outdir)
    
    return wfr

def main():
    wfr = get_wfr("focus", "../../data/h5/focus_NKB_9-2keV_0-25nC.hdf5", focus = 'nano')
    
    #wfr = Wavefront()
    #wfr.load_hdf5("../../data/h5/MKB_MHE_9-2keV_0-25nC.hdf5")
    
    zmin = getRequiredDistance(wfr, 13.5e-06)
    getImageProperties(wfr, 1, 1e-06, [2156, 2048])
    
if __name__ == "__main__":
    main()