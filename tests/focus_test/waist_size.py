#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:45:13 2020

A scratchpad for solving for beamline propagation


    
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
import matplotlib as mpl
from copy import copy

from wpg.beamline import Beamline
from wpg.srwlib import SRWLOptD

from model.src.coherent import construct_SA1_wavefront
from model.beamline.structure import BeamlineModel, propagation_parameters
from model.materials.load_refl import load_refl, get_refl

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from wpg.wpg_uti_wf import plot_intensity_qmap as plotPhase
from wpg.wpg_uti_wf import calculate_fwhm

from matplotlib import pyplot as plt 

def getFocusSize(ekev, q, focus = 'micron'):
    
    wfr = construct_SA1_wavefront(1024, 1024, ekev, q)
    
    spb = BeamlineModel(overwrite_mirrors =  True)
    
    spb.setupHOMs(ekev, 2.2e-03)
    spb.setupKBs(ekev, 3.5e-03)
    
    spb.build_elements(focus = focus)
    spb.build_beamline(focus = focus)
    bl = spb.get_beamline()
    
    bl.propagate(wfr)
    
    meas = calculate_fwhm(wfr)
    
    fwhm_x = meas['fwhm_x']
    fwhm_y = meas['fwhm_y']
    
    return fwhm_x, fwhm_y

def plot(energyrange, data, estr, qstr, focus, outdir = "../../tmp/"):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    
    if focus == "micron":
        focstr = "Micron-KB "
    elif focus == "nano":
        focstr = "Nano-KB "
        
    ax.set_title(focstr + "Focal Spot Size: q = {} nC".format(qstr))
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Focal Spot Width ($\mu$m)")

    ax.scatter(energyrange, data[0], c = 'r')
    ax.scatter(energyrange, data[1], c = 'b')
    
    ax.set_ylim()
    
    plt.legend(['Horizontal Beam Width', 'Vertical Beam Width'])
    
    fig.savefig(outdir + "FocalSpotSize_{}_{}nC.png".format(focus, qstr))
    
if __name__ == '__main__':

    focus = 'micron'
    
    energyrange = np.linspace(6.0, 15.0, 9)
    qrange = [0.1, 0.25, 0.50]
    
    for q in qrange:
        qstr = str(q)
        _X, _Y = [], []
        for ekev in energyrange:
        
            estr = str(ekev).replace(".","-")
        
        
           
            
            fwhm_x, fwhm_y = getFocusSize(ekev, q, focus = focus)
            
            _X.append(fwhm_x*1e6)
            _Y.append(fwhm_y*1e6)
            
        plot(energyrange, [_X, _Y], estr, qstr, focus, outdir = "waist_size/")
        
        
