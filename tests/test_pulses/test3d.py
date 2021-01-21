#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:48:50 2020

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
from wpg import srwlib
from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import integral_intensity, calculate_fwhm, look_at_q_space
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from model.src.coherent import construct_SA1_wavefront
from model.beamline.structure import BeamlineModel, propParams
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture
from wpg.generators import build_gauss_wavefront
from wpg.beamline import Beamline

def constructPulse():
    
    wfr = Wavefront(build_gauss_wavefront(512, 512, 10, 5.0, -400e-06, 400e-06, -400e-06, 400e-06, 1e-15, 5e-06, 5e-06, 19))
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    #look_at_q_space(wfr)
    return wfr



def simpleProp(wfr):
    
    print(calculate_fwhm(wfr))
    bl = Beamline()
    #bl.append(Aperture('c','a', 500e-06),propParams(1,1,1,1,mode = 'normal'))
    bl.append(Drift(100), propParams(1,1,1,1,mode = 'quadratic'))
    
    #bl.append(Drift(100), propParams(1,1,1,1,mode = 'quadratic'))
    bl.propagate(wfr)
    plotIntensity(wfr)
    print(calculate_fwhm(wfr))

if __name__ == "__main__":
    

    
    #wfr = constructPulse()
    wfr = Wavefront()
    wfr.load_hdf5("../../data/h5/gauss.h5")
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    plotIntensity(wfr)
    
    simpleProp(wfr)
    
    #comparePulses(wfr, cwfr, "/data/")