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

from model.src.coherent import coherentSource
from model.beamline.structure import BeamlineModel, propParams
from model.materials.load_refl import load_refl, get_refl

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from wpg.wpg_uti_wf import plot_intensity_qmap as plotPhase

from matplotlib import pyplot as plt 

def sliceFocus(wfr, ekev, focus = 'micron', nslices = 500, axisName = 'x', outdir = None):
    
    spb = BeamlineModel(overwrite_mirrors = False)
    spb.setupHOMs(ekev)
    spb.buildElements(focus = focus)
    spb.buildBeamline(focus = focus)
    
    el_n = len(spb.bl.propagation_options[0]['optical_elements'])-1
    
    slice_interval = copy(spb.bl.propagation_options[0]['optical_elements'][el_n].L/1000) 
    spb.bl.propagation_options[0]['optical_elements'][el_n].L *= 0.75
    spb.bl.propagation_options[0]['propagation_parameters'][el_n] = propParams(1/5,1,1/5,1, mode = 'quadratic')
    bl = spb.get_beamline()
    bl.propagate(wfr)
    
    if axisName == 'x':
        data_focslice = np.zeros([nslices, np.shape(wfr.data.arrEhor)[0]])
    elif axisName == 'y':
        data_focslice = np.zeros([nslices, np.shape(wfr.data.arrEhor)[1]])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    
    
    for n in range(nslices):
        print("Slice {}/{}".format(n+1, nslices))
        
        bl = Beamline()
        bl.append(SRWLOptD(slice_interval), propParams(1,1,1,1, mode = 'normal'))
        bl.propagate(wfr)
        
        data_focslice[-n-1, :] = wfr.get_profile_1d()[0]
        
    
    y = np.linspace(spb.bl.propagation_options[0]['optical_elements'][el_n].L,
                    spb.bl.propagation_options[0]['optical_elements'][el_n].L + nslices*slice_interval,
                    nslices)
    
    if axisName == 'x':
        extent = [wfr.params.Mesh.xMin*1e6, wfr.params.Mesh.xMax*1e6,
                          y.min(), y.max()]
    elif axisName == 'y':
        extent = [wfr.params.Mesh.yMin*1e6, wfr.params.Mesh.yMax*1e6,
                  y.min(), y.max()]
    ax.imshow(data_focslice,
              extent= extent,
              aspect = 'auto',
              norm=mpl.colors.LogNorm())
    
    if focus == "micron":
        ax.set_title("Micron Focus Location")
        ax.set_ylabel("Longitudinal Distance from MKB (m)")
        
        if axisName == 'x':
            x = np.linspace(wfr.params.Mesh.xMin*1e6, wfr.params.Mesh.xMax*1e6, wfr.params.Mesh.nx)
            ax.set_xlabel("x ($\mu$m)")
        elif axisName == 'y':
            x = np.linspace(wfr.params.Mesh.yMin*1e6, wfr.params.Mesh.yMax*1e6, wfr.params.Mesh.yx)
            ax.set_xlabel("y ($\mu$m)")
            
    
    ax.plot(x, np.ones(x.shape)*spb.bl.params["df"]["distance"], 'r--')
    plt.legend(["Design Focal Plane"])
    
    plt.show()
    
    estr = str(ekev).replace(".","-")
    
    if outdir is not None:
        fig.savefig(outdir + "{}_focus_slice_{}_{}".format(focus, estr, axisName))
        
if __name__ == '__main__':
    
    ekev = 9.2
    
    wfr = coherentSource(1024, 1024, ekev, 0.5)
    
    sliceFocus(wfr = wfr, ekev = ekev, focus = 'micron', nslices = 500, axisName = 'x', outdir = "focus_test/")
    sliceFocus(wfr = wfr, ekev = ekev, focus = 'micron', nslices = 500, axisName = 'y', outdir = "focus_test/")