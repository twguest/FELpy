#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:45:13 2020

A scratchpad for testing beamline methods


@author: twguest
"""


###############################################################################
import sys
 
sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################

import time

from model.src.coherent import coherentSource
from model.beamline.structure import BeamlineModel

    
def testNano(wfr, outdir = None, toggle = 'on'):
     
    print("Testing SPB-SFX Instrument in Nano-KB Config")
    spb = BeamlineModel()

    print("Setting Up Beamline")
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    spb.mirrorProfiles(toggle = toggle, overwrite = False)
    
    print("Building Beamline")
    spb.buildElements(focus = "nano")
    spb.buildBeamline(focus = "nano")
  
    bl = spb.get_beamline()
    
    print("Propagating Beamline")
    s = time.time()
    bl.propagateSeq(wfr, outdir)
    f = time.time()
    print("Propagation Finished")
    print("Propagating Time: {} s".format(f-s))


def plotMirrorProfiles(outdir):

    spb = BeamlineModel()
    spb.mirrorProfiles(toggle = "on", overwrite = True)
    
    spb.plotMirrorProfile("HOM1", outdir = outdir)
    spb.plotMirrorProfile("HOM2", outdir = outdir)
    spb.plotMirrorProfile("NHE", outdir = outdir)
    spb.plotMirrorProfile("NVE", outdir = outdir)


    
if __name__ == '__main__': 
    
    wfr = coherentSource(1024, 1024, 4.96, 0.250)
    
    testNano(wfr, outdir = "/opt/spb_model/data/beamlineTests/prop/")
