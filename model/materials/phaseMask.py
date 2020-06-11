#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:43:15 2020

A wrapper for WPG setup_opt_surf_height_2D, which takes direct phase shift values
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
import numpy as np

from wpg.srwlib import SRWLOptD as Drift
from model.src.coherent import coherentSource
from model.beamline.structure import propParams
from wpg.beamline import Beamline
from wpg.srwlib import srwl_opt_setup_surf_height_2d as OPD
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from matplotlib import pyplot as plt
from model.materials.material_utils import add_extent

def pltPhase(wfr):
    
    phase = wfr.get_phase()[:,:,0]
    plt.imshow(phase, cmap = 'hsv')
    plt.show()

def phaseMask(phaseshift, extent, wav, _ang = 0, outdir = None, maskName = None):
    """
    :param phaseshift: 2D array of desired phase-shifts
    :param extent: extent of phase-mask array in realspace
    :param wav: radiation wavelength 
    """
    
    height_error = (phaseshift*wav)/(2*np.pi)
    height_error = add_extent(height_error, extent)

    if outdir is not None:
        
        if maskName is not None:
            outdir = outdir + maskName + ".dat"
        else:
            outdir = outdir + "phasemask.dat"
        
        np.savetxt(outdir, height_error)
                   
    return OPD(height_error,
               _dim = 'x',
               _ang = np.pi/2,
               _refl = 1,
               _x = 0, _y = 0)

if __name__ == "__main__":
    

    wfr = coherentSource(1024,1024,3, 1)
    pltPhase(wfr)
    phaseshift = np.random.rand(200,200)*2*np.pi

    opd = phaseMask(phaseshift, [5e-03, 5e-03], wfr.params.wavelength*100)
    bl = Beamline()
    
    bl.append(opd,propParams(1,1,1,1))
    bl.append(Drift(2), propParams(1,1,1,1))
    bl.propagate(wfr)
    pltPhase(wfr)
    plotIntensity(wfr)
    
    