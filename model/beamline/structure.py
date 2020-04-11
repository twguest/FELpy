#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:43:15 2020

@author: twguest
"""

import sys
sys.path.append("/opt/WPG/")
sys.path.append("/opt/spb_model")

import numpy as np

from wpg import srwlpy
 
from wpg.beamline import Beamline

from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture
from wpg.srwlib import SRWLOptMirEl as MirEl


from wpg.optical_elements import Mirror_plane_2d as MirPl

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from model.src.coherent import coherentSource

from wpg.misc import get_wavelength, sampling
from wpg.wpg_uti_wf import check_sampling
def propParams(sx, zx, sy, zy, mode = "normal"):
    """
    wrapper for propagation parameters
    
    :param zx: horizontal scaling factor 
    :param rx: horizontal zoom factor
    :param zy: vertical scaling factor
    :param ry: vertical zoom factor
    :param mode: normal, semi-analytical, converge or diverge
    
    :return propagation parameters:
    """
    
    if mode == "normal":
        m = 0
    elif mode == "quadratic":
        m = 1
    elif mode == "diverge":
        m = 3
    elif mode == "converge":
        m = 4
    
    return [0,0,1,m,0,sx,zx/sx,sy,zy/sy,0,0,0]


def genMirrorSurface(wfr, outdir, mode = 'Flat'):
    
    if mode == 'Flat':
        
        nx, ny = wfr.params.Mesh.nx, wfr.params.Mesh.ny
        
        surface = np.zeros((max(nx,ny),3))
        surface[:,0] = np.linspace(wfr.params.Mesh.xMin, wfr.params.Mesh.xMax,
                                   nx)
        surface[:,1] = np.linspace(wfr.params.Mesh.yMin, wfr.params.Mesh.yMax,
                                   ny)
        
    np.savetxt(outdir+"mir_"+ mode +".txt", surface)
        


def NyquistFreq(wfr, z, scaler, axisName = 'x'):
    """
    calculate the required zoom to satisfy nyquist-shannan sampling
    conditions in the propagated plane
    
    for calculation of req. pixel size see: wpg/misc/sampling 
    
    :param wavelength: source wavelength [m]
    :param z: propagation distance [m]
    :param scaler: field scale multiplier (ie., extent_f = scale*extent_i)
    
    :return zoom: zoom modifier for propagation step
    """
    wavelength = get_wavelength(wfr)

    if axisName == 'x':
        d = wfr.pixelsize()[0]
        extent = wfr.params.Mesh.xMax - wfr.params.Mesh.xMin
    elif axisName == 'y':
        d = wfr.pixelsize()[1]
        extent = wfr.params.Mesh.yMax - wfr.params.Mesh.yMin

    W = (wavelength*z)/(extent*scaler*d*np.sqrt(2))
    
    return W/d
    
    
def config(beamline = "micro", screens = True):
    
    bl = Beamline()
    
    
    HOM1 = MirPl('x', 1.0e-03, 1, 1, "/opt/spb_model/data/hom1_mir_Flat",
                 x0 = 0, y0 = 0, bPlot = True)

    pass
    

class spb:
    
    def __init__(self):
        """
        """
        pass
    
    def setup(self):
        """
        setup beamline elements, can be micro or nano
        """
        pass

if __name__ == '__main__':
    
    config()
# =============================================================================
#     wfr = coherentSource(1048, 1048, 9.2, 0.1)
#     print(check_sampling(wfr))
#     W = NyquistFreq(wfr, 1, 1)
#     print(W)
#     plotIntensity(wfr)
#     d1 = Drift(246.5)
#     bl = Beamline()
#     bl.append(d1, propParams(5, 2.5, 5, 2.5))
#     bl.propagate(wfr)
#     print(check_sampling(wfr))
#     plotIntensity(wfr)
# =============================================================================
    
