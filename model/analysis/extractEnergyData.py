#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:32:40 2020

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

import os

import numpy as np

import wpg.srwlib as srwlib

from wpg.wpg_uti_wf import calc_pulse_energy, getOnAxisPowerDensity, getCentroid, get_profile_1d
from wpg.wavefront import Wavefront



indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/out/"
directoryName = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/legacyData/"


def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def loadWavefront(fname):
    
    wfr = Wavefront()
    wfr.load_hdf5(indir + fname)
    
    return wfr
    
def storeCentroid(wfr, fname):
    
    
    intDir = directoryName + "iCentroid/"
    slicedDir = directoryName + "sCentroid/"
    
    mkdir_p(intDir)
    mkdir_p(slicedDir)
    

    cI = getCentroid(wfr, mode = 'integrated')
    np.save(intDir + fname, cI)
    del cI

    cS = getCentroid(wfr, mode = 'pulse')
    np.save(slicedDir + fname, cS)
    del cS    
    
    print("Centroids Stored: {}".format(fname))
    
    
def storeSpectrum(wfr, fname):
    
    timeDir = directoryName + "tSpectrum/"
    freqDir = directoryName + "fSpectrum/"
    
    mkdir_p(timeDir)
    mkdir_p(freqDir)

    sz0 = getOnAxisPowerDensity(wfr, spectrum = False)
    np.save(timeDir + fname, sz0)
    del sz0

    sz1 = getOnAxisPowerDensity(wfr, spectrum = True)
    np.save(freqDir + fname, sz1)
    del sz1
    
    print("Spectrums Stored: {}".format(fname))
   
def storePhotonEnergy(wfr, fname):
    
    enDir = directoryName + "pulseEnergy/"
    mkdir_p(enDir)
    
    effDir = directoryName + "systemEff/"
    mkdir_p(effDir)
    
    

    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')
    pulseEn, photons_per_pulse = calc_pulse_energy(wfr)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
    effPhotons = photons_per_pulse / wfr.custom_fields['source']['nPhotons']
    effEn = pulseEn / wfr.custom_fields['source']['pulseEn']
    
    
    np.save(enDir + fname, [pulseEn, photons_per_pulse])
    np.save(effDir + fname, [effPhotons, effEn])
    
    
    print("Photon Energy Saved: {}".format(fname))

def storeProfiles(wfr, fname):
    
    profDir = directoryName + "profiles/"
    mkdir_p(profDir)
    

    profile = get_profile_1d(wfr)
    
    np.save(profDir + fname, profile)
        
    print("1D Profiles Stored: {}".format(profDir))

def integratedAnalysis(indir, outfile):
    
    data = []
    
    for f in os.listdir(indir):
        
        print("Processing: {}".format(f))
        d = []
        
        wfr = loadWavefront(f)
        cen = getCentroid(wfr, mode = 'integrated')
           
        srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')
        pulseEn, photons_per_pulse = calc_pulse_energy(wfr)
        srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
        d.append([f, cen, pulseEn, photons_per_pulse])
        data.append(d)
    
    data = np.asarray(data)
    np.save(outfile, data)
    
def main(fname):
    
    mkdir_p(directoryName)
    
    wfr = loadWavefront(fname)

    storeCentroid(wfr, fname)
    storeSpectrum(wfr, fname)
    storePhotonEnergy(wfr, fname)
    storeProfiles(wfr, fname)
    
if __name__ == '__main__':
    
    indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/out/"
    outfile = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/integratedEnergyAnalysis.npy"
    
    integratedAnalysis(indir, outfile)
        
