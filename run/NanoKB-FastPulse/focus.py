#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 09:55:11 2020

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
from matplotlib import pyplot as plt
import matplotlib as mpl
from wpg.beamline import Beamline
from wpg.wavefront import Wavefront
from model.tools import constructPulscmdclass={'install': MyInstall}e
from model.beamline.structure import BeamlineModel, propParams
from wpg.srwlib import SRWLOptD as Drift
import numpy as np
from copy import deepcopy
def getSPB(wfr):
    
    focus = "nano"
    
    spb = BeamlineModel()
    
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "on", aperture = True, overwrite = False)
    
    spb.buildElements(focus)
    spb.buildBeamline(focus)
    spb.cropBeamline("NVE")
    print("Setup Beamline")
    
    spb.scale(wfr) 
    print("Source Scaled")
    
    bl = spb.get_beamline()
    print(bl)
    return bl
cmdclass={'install': MyInstall}
def propThruFocus(wfr):
    
    outdir = "../../tmp/"
    
    bl = getSPB(wfr)
    bl.propagate(wfr)
    
    
    propDist = 100e-06
    focDist = 2.2
    nFocSlices = 100
    prop  = focDist-(propDist/2)
    sliceDist = propDist/nFocSlices
    
    merProfx = np.ones([1024, nFocSlices])
    #merProfy = np.ones([wfr.params.Mesh.ny, nFocSlices])
    
    bl = Beamline()
    bl.append(Drift(prop), propParams(1,1,1,1, 'converge'))
    bl.propagate(wfr)
    wfr.store_hdf5(outdir + "testwfr.h5")

    for itr in range(nFocSlices):
        
        bl = Beamline()
        cwfr = Wavefront()
        
        print("loading wfr")
        cwfr.load_hdf5(outdir + "testwfr.h5")

        bl.append(Drift(sliceDist*(itr+1)), propParams(1,1,1,1, mode = 'normal'))
        bl.propagate(cwfr)
        
        ix, iy = getIntensity(wfr)
        merProfx[:, itr] = ix
        #merProfy[:, itr] = iy
        print(cwfr.get_limits())
        plotIntensity(cwfr)
    return merProfx#, merProfy
    
    
def getIntensity(wfr):
    
    ii = wfr.get_intensity().sum(axis = -1)
    
    ix = ii[:, ii.shape[1]//2].T #.sum(axis = 1)
    iy = ii.sum(axis = 0)

    
    return ix, iy
    
if __name__ == '__main__':
    from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
    from model.src.coherent import coherentSource
    wfr = coherentSource(1024,1024,4.9,0.25)
    merProfx= propThruFocus(wfr)
    
    plt.imshow(merProfx,
               aspect = 'auto')