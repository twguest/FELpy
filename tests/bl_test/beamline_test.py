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

from model.src.coherent import construct_SA1_wavefront
from model.beamline.structure import Instrument

    
def testNano(wfr, outdir = None, toggle = 'on'):
     
    print("Testing SPB-SFX Instrument in Nano-KB Config")
    spb = Instrument()

    print("Setting Up Beamline")
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    spb.mirrorProfiles(toggle = toggle, overwrite = False)
    
    print("Building Beamline")
    spb.build_elements(focus = "nano")
<<<<<<< HEAD
    spb.build_beamline(focus = "nano")
=======
    spb.buildBeamline(focus = "nano")
>>>>>>> 108cfb9b6fc97d3841ee1db54862523eee5b184e
  
    bl = spb.get_beamline()
    
    print("Propagating Beamline")
    s = time.time()
    bl.propagate_sequential(wfr, outdir)
    f = time.time()
    print("Propagation Finished")
    print("Propagating Time: {} s".format(f-s))


def plotMirrorProfiles(outdir):

    spb = Instrument()
    spb.mirrorProfiles(toggle = "on", overwrite = True)
    
    spb.plotMirrorProfile("HOM1", outdir = outdir)
    spb.plotMirrorProfile("HOM2", outdir = outdir)
    spb.plotMirrorProfile("NHE", outdir = outdir)
    spb.plotMirrorProfile("NVE", outdir = outdir)


    
if __name__ == '__main__': 
    
    wfr = construct_SA1_wavefront(1024, 1024, 4.96, 0.250)
    
    testNano(wfr, outdir = "/opt/spb_model/data/beamlineTests/prop/")
