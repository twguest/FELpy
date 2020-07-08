#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 15:09:54 2020

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
indir = "../../data/h5/NanoKB-Pulse/in/"
outdir = "../../data/h5/NanoKB-Pulse/out/"


def storeWavefrontInfo(wfr):
    
    sz0 = getOnAxisPowerDensity(wfr, spectrum = False)
    sz1 = getOnAxisPowerDensity(wfr, spectrum = True)
    
    fwhm = calculate_fwhm(wfr)
    
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')
    pulseEn, photons_per_pulse = calc_pulse_energy(wfr)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
    divergence = calcDivergence(wfr)
    centroid = getCentroid(wfr)
    
    wfr.custom_fields['/source/t_spectrum'] = sz0
    wfr.custom_fields['/source/f_spectrum'] = sz1
    
    wfr.custom_fields['/source/xFWHM'] = fwhm['fwhm_x']
    wfr.custom_fields['/source/yFWHM'] = fwhm['fwhm_y']
    
    wfr.custom_fields['/source/divX'] = divergence[0]
    wfr.custom_fields['/source/divX'] = divergence[1]
    
    wfr.custom_fields['/source/pulseEn'] = pulseEn
    wfr.custom_fields['/source/nPhotons'] = photons_per_pulse
    
    wfr.custom_fields['/source/centroid'] = centroid


def getSPB(wfr):
    
    spb = BeamlineModel()
    
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "on", aperture = True, overwrite = False)
    
    spb.buildElements(focus)
    spb.buildBeamline(focus)
    spb.scale(wfr) 
    
    
    bl = spb.get_beamline()
    
    return bl

def propagatePulses(fname):
    
    
    wfr = Wavefront()
    wfr.load_hdf5(indir + fname)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    storeWavefrontInfo(wfr)
    
 
    bl = getSPB(wfr)
    bl.propagate(wfr)
    
    wfr.store_hdf5(outdir + fname)


if __name__ == '__main__':

    MPI = True
      
    f = listdir(indir)
    
    if MPI:
        
        cores = 35 ## DESY-MAXWELL ONLY
        print("Using  {} Cores".format(cores))
        print("Available Cores: {}".format(multiprocessing.cpu_count()))
        p = multiprocessing.Pool(cores)
        
        for _ in tqdm(p.imap_unordered(propagatePulses, f), total = len(f)):
            pass
        
        p.close()
        p.join()
        
    else:
        for fname in f:
            propagatePulses(fname)
            
    