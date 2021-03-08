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

from model.beamline.structure import propagation_parameters
#from model.beamline.structure import BeamlineModel

from wpg import srwlib

from wpg.srwlib import SRWLOptD as Drift

from wpg.wavefront import Wavefront
from wpg.beamline import Beamline

from wpg.wpg_uti_wf import calc_pulse_energy, calculate_fwhm, get_axial_power_density

from wpg.misc import calcDivergence


def storeWavefrontInfo(wfr):
    
    sz0 = get_axial_power_density(wfr, spectrum = False)
    sz1 = get_axial_power_density(wfr, spectrum = True)
    
    fwhm = calculate_fwhm(wfr)
    
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')
    pulseEn, photons_per_pulse = calc_pulse_energy(wfr)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
    divergence = wfr.get_divergence()
    
    
    wfr.custom_fields['/source/t_spectrum'] = sz0
    wfr.custom_fields['/source/f_spectrum'] = sz1
    
    wfr.custom_fields['/source/xFWHM'] = fwhm['fwhm_x']
    wfr.custom_fields['/source/yFWHM'] = fwhm['fwhm_y']
    
    wfr.custom_fields['/source/divX'] = divergence[0]
    wfr.custom_fields['/source/divX'] = divergence[1]
    
    wfr.custom_fields['/source/pulseEn'] = pulseEn
    wfr.custom_fields['/source/nPhotons'] = photons_per_pulse
    
   

def getSimpleBl():
    
    bl = Beamline()    
    
    bl.append(Drift(10), propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
    
    return bl
    

def propagatePulses(indir, fname, outdir):
    
    wfr = Wavefront()
    wfr.load_hdf5(indir + fname)
    
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    storeWavefrontInfo(wfr)
    
    
    bl = getSimpleBl()
    
    bl.propagate(wfr)
    
    wfr.store_hdf5(outdir + fname)
    
def testStorage(indir, fname, N = 25):
    
    for n in range(N):
        wfr = Wavefront()
        wfr.load_hdf5(indir + fname)
    
        print(wfr.custom_fields['source'])

        
def main(N = 2):
    

    indir = "../../data/tests/pulseTests/gaussianSource/"
    
    outdir = "../../data/tests/pulseTests/gaussianOut/simple/"
    
    
    
    for n in range(N):
        fname = "gsn_{}.h5".format(n)
        print("Propagating Pulse {} of {}".format(n+1, N))    
        propagatePulses(indir, fname , outdir )
        
    testStorage(outdir, fname, N)
        
if __name__ == '__main__':
    main(2)
    