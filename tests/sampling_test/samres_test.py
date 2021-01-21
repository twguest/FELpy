#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 12:31:27 2020

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

from model.src.coherent import construct_SA1_wavefront
from model.beamline.structure import BeamlineModel
from model.beamline.samplingCalc import samresPlot



if __name__ == "__main__":
    
    list_of_elements = ["MVE"]
    
    detector_pixel_size = 1e-06
    N_pixels = 1024
    detector_width = detector_pixel_size*N_pixels
    
    ekev = 9.2
    q = 0.25
    
    fdir = "samres_test"
    if os.path.exists(fdir) is False:
        os.mkdir(fdir)
    
    for element in list_of_elements:
            
        wfr = construct_SA1_wavefront(1024, 1024, ekev, q)
        
        spb = BeamlineModel()
        spb.setupHOMs(ekev, 2.2e-03)
        spb.setupKBs(ekev, 3.5e-03)
        spb.mirrorProfiles(toggle = "on", overwrite = False)
        
        
        
        spb.buildElements(focus = "micron")
        spb.buildBeamline(focus = "micron")
        spb.cropBeamline(element1 = element)
        
        bl = spb.get_beamline()
        
        bl.propagate(wfr)
        
        samresPlot(wfr, 12.0, detector_width, detector_pixel_size,
                   outdir = fdir + "/" + element +"_")
        
        
    f = open(fdir+"/info.txt","w")
    f.write("Energy: {} keV/n".format(ekev))
    f.write("Beam Charge: {} nC/n".format(q))
    f.write("Detector Pixel Size: {} m/n".format(detector_pixel_size))
    f.write("Detector Width: {} m".format(detector_width))
    f.close()