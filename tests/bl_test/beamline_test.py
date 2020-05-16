#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:45:13 2020

A scratchpad for testing beamline methods


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


def testMicron(ekev, q, outdir = "../../tmp", toggle = 'on'):
    
    wfr = coherentSource(1024, 1024, ekev, q)
    
    spb = BeamlineModel()
    
    spb.setupHOMs(ekev, 2.2e-03)
    spb.setupKBs(ekev, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "on", overwrite = False)
    spb.buildElements(focus = "micron")
    spb.buildBeamline(focus = "micron")
    bl = spb.get_beamline()
    
    bl.propagateSeq(wfr, outdir)
    
def testNano(ekev, q, outdir = "../../tmp", toggle = 'on'):
    
    wfr = coherentSource(1024, 1024, ekev, q)
    
    spb = BeamlineModel()
    
    spb.setupHOMs(ekev, 2.2e-03)
    spb.setupKBs(ekev, 3.5e-03)
    spb.mirrorProfiles(toggle = "on", overwrite = False)
    spb.buildElements(focus = "nano")
    spb.buildBeamline(focus = "nano")
    bl = spb.get_beamline()
    
    bl.propagateSeq(wfr, outdir)


def testSurfaceRoughness(energyrange, q, outdir):

    for ekev in energyrange:
        testMicron(ekev, q, outdir = outdir, toggle = 'off')
    for ekev in energyrange:
        testMicron(ekev, q, outdir = outdir, toggle = 'on')
    for ekev in energyrange:
        testNano(ekev, q, outdir = outdir, toggle = 'off')
    for ekev in energyrange:
        testNano(ekev, q, outdir = outdir, toggle = 'on')
 

def plotMirrorProfiles(outdir):

    spb = BeamlineModel()
    spb.mirrorProfiles(toggle = "on", overwrite = True)
    
    spb.plotMirrorProfile("HOM1", outdir = outdir)
    spb.plotMirrorProfile("HOM2", outdir = outdir)
    spb.plotMirrorProfile("MHP", outdir = outdir)
    spb.plotMirrorProfile("MVP", outdir = outdir)

if __name__ == '__main__': 
    
    plotMirrorProfiles("../../tests/bl_test/test_surface_roughness")
    testSurfaceRoughness([6.0, 9.2, 12.0],0.25, outdir = "../../tests/bl_test/test_surface_roughness")
    

    