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

from wpg.wavefront import Wavefront

def testMicron(ekev, q, outdir = "../../tmp", toggle = 'on'):
    
    wfr = coherentSource(1024, 1024, ekev, q)
    
    spb = BeamlineModel()
    
    spb.setupHOMs(ekev, 2.2e-03)
    spb.setupKBs(ekev, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "off", aperture = False, overwrite = True)
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


def testBeamlineGeometry(ekev, q, focus = 'nano'):
    
    wfr = coherentSource(1024, 1024, ekev, q)
    
    spb = BeamlineModel()
    
    spb.setupHOMs(ekev, 2.2e-03)
    spb.setupKBs(ekev, 3.5e-03)
    spb.mirrorProfiles(toggle = "on", overwrite = False)
    
    
    
    spb.buildElements(focus = "nano")
    spb.buildBeamline(focus = "nano")
    spb.cropBeamline("NVE")
    bl = spb.get_beamline()
    
    bl.propagateSeq(wfr, outdir = "../../tmp")

def testEllipticalRoughness(ekev, q, focus = "micron"):
    wfr = coherentSource(1024, 1024, ekev, q)
    
    spb = BeamlineModel()
    
    spb.setupHOMs(ekev, 2.2e-03)
    spb.setupKBs(ekev, 3.5e-03)
    spb.mirrorProfiles(toggle = "on", overwrite = False)
    
    
    
    spb.buildElements(focus)
    spb.buildBeamline(focus)
    
    if focus == 'nano':
        spb.cropBeamline("NVE")
        
    elif focus == 'micron':
        #spb.cropBeamline("MVE")
        pass
    bl = spb.get_beamline()
    
    bl.propagateSeq(wfr, outdir = "../../tmp")
    return wfr

def getMicronFocalPlaneBeam(outdir, wfrdir = None):
    
    if wfrdir == None:    
        wfr = coherentSource(1024, 1024, 9.2, 0.250)
    else:
        wfr = Wavefront()
        wfr.load_hdf5(wfrdir)
        
    spb = BeamlineModel()
    
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "on", overwrite = False)
    
    
    
    spb.buildElements(focus = "micron")
    spb.buildBeamline(focus = "micron")
    
    spb.scale_wfr(wfr)    
    
    bl = spb.get_beamline()
    
    bl.propagate(wfr)
    
    wfr.store_hdf5(outdir)

def testOnOff():
    

    wfr = coherentSource(1024, 1024, 9.2, 0.250)

    spb = BeamlineModel()
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    
    ### no apertures w/ profile
    #spb.mirrorProfiles(toggle = "on", aperture = False, overwrite = True)
    
    ### apertures no profile
    #spb.mirrorProfiles(toggle = "off", aperture = True, overwrite = True)
    
    ### apertures and profile
    spb.mirrorProfiles(toggle = "on", aperture = True, overwrite = True)
    
    spb.buildElements(focus = "micron")
    spb.buildBeamline(focus = "micron")
    
    
    spb.cropBeamline("MVE")
    bl = spb.get_beamline()
    
    bl.propagateSeq(wfr)

def testwfrScale():
    wfr = coherentSource(500, 500, 9.2, 0.25)
    
    spb = BeamlineModel()
    spb.scale(wfr, 2000)
    print(wfr.params.Mesh.nx)
    
if __name__ == '__main__': 
    #testNano(9.2, 0.25)
    #testMicron(9.2, 0.25)
    #testOnOff()
    #getMicronFocalPlaneBeam("../../data/input/micronfoc_9-2keV_250pC.hdf5")
    #getMicronFocalPlaneBeam("../../data/input/micronfoc_8-86keV_250pC.hdf5", wfrdir = "../../data/h5/8_86keV_0250pC.h5")
    testwfrScale()