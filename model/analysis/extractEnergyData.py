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

from wpg.wpg_uti_wf import calc_pulse_energy, calculate_fwhm, getOnAxisPowerDensity, getCentroid
from wpg.wavefront import Wavefront
from wpg.generators import build_gauss_wavefront


fname = sys.argv[1]
indir = "/gpfs/exfel/data/user/guest/spb_model/data/h5/NanoKB-Pulse/out/"
directoryName = "/gpfs/exfel/data/user/guest/spb_model/data/h5/NanoKB-Pulse/data/"


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
    cS = getCentroid(wfr, mode = 'pulse')
    
    np.save(intDir + fname, cI)
    np.save(intDir + fname, cS)
    
    print("Centroids Stored: {}".format(fname))
    
    del cI, cS

def storeSpectrum(wfr, fname):
    
    timeDir = directoryName + "tSpectrum/"
    freqDir = directoryName + "fSpectrum/"
    
    sz0 = getOnAxisPowerDensity(wfr, spectrum = False)
    sz1 = getOnAxisPowerDensity(wfr, spectrum = True)
    
    np.save(timeDir + fname, sz0)
    np.save(freqDir + fname, sz1)
    
    print("Spectrums Stored: {}".format(fname))
   
def storePhotonEnergy(wfr, fname):
    
    enDir = directoryName + "pulseEnergy/"
    
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')
    pulseEn, photons_per_pulse = calc_pulse_energy(wfr)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
    np.save(endir + fname, [pulseEn, photons_per_pulse])
    
    print("Photon Energy Saved: {}".format(fname))

def storeProfiles(wfr, fname):
    
    profDir = directoryName + "profiles/"
    
    profile = wfr.get_profile_1d()
    
    np.save(profDir + fname, profile)
    
    print("1D Profiles Stored: {}".format(profDir))
    
def main(fname):
    
    mkdir_p(directoryName)
    
    wfr = loadWavefront(fname)

    storeCentroid(wfr, fname)
    storeSpectrum(wfr, fname)
    storePhotonEnergy(wfr, fname)

if __name__ == '__main__':
    
    main(fname)
        
