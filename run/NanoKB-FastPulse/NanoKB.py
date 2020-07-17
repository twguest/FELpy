#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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

from model.beamline.structure import BeamlineModel

from wpg import srwlib

from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import calc_pulse_energy, calculate_fwhm, getOnAxisPowerDensity, getCentroid
from wpg.misc import calcDivergence

from os import listdir
from tqdm import tqdm

focus = "nano"
indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/in/"
outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/noProfile/"

def storeWavefrontInfo(wfr):
    
    
    sz0 = getOnAxisPowerDensity(wfr, spectrum = False)
    print("Got On Axis Power Density - Temporal")
    sz1 = getOnAxisPowerDensity(wfr, spectrum = True)
    print("Got On Axis POwer Density - Spectral")
    fwhm = calculate_fwhm(wfr)
    print("Got FWHM")
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')
    pulseEn, photons_per_pulse = calc_pulse_energy(wfr)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
    divergence = calcDivergence(wfr)
    print("Got Divergence")
    centroid = getCentroid(wfr)
    print("Got Centroid")
    wfr.custom_fields['/source/t_spectrum'] = sz0
    wfr.custom_fields['/source/f_spectrum'] = sz1
    
    wfr.custom_fields['/source/xFWHM'] = fwhm['fwhm_x']
    wfr.custom_fields['/source/yFWHM'] = fwhm['fwhm_y']
    
    wfr.custom_fields['/source/divX'] = divergence[0]
    wfr.custom_fields['/source/divX'] = divergence[1]
    
    wfr.custom_fields['/source/pulseEn'] = pulseEn
    wfr.custom_fields['/source/nPhotons'] = photons_per_pulse
    
    wfr.custom_fields['/source/centroid'] = centroid
    print("Wrote History to File")

def getSPB(wfr):

    spb = BeamlineModel()
    
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "off", aperture = True, overwrite = False)
    
    spb.buildElements(focus)
    spb.buildBeamline(focus)
    print("Setup Beamline")
    
    spb.scale(wfr) 
    print("Source Scaled")
    
    bl = spb.get_beamline()
    
    return bl
 
def propagatePulses(fname):
    
    
    wfr = Wavefront()
    wfr.load_hdf5(indir + fname)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    storeWavefrontInfo(wfr)
    
 
    bl = getSPB(wfr)
    print("propagating pulse, hold on...")
    bl.propagate(wfr)
    
    wfr.store_hdf5(outdir + fname)
    print("Writing to File: Done")
def main(fname):
    
    propagatePulses(fname)

if __name__ == '__main__':
    fname = sys.argv[1]
    main(fname)