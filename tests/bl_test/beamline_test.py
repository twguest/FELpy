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

from model.src.coherent import coherentSource
from model.beamline.structure import config

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 


if __name__ == '__main__':

    outdir = "../../output/"
    #outdir = None
    
    print("Testing Beamline Propagation")
    
    
    wfr = coherentSource(1024,1024,6,1.0)
    
    bl = config(focus = "micron")
    
    bl.propagateSeq(wfr, outdir = outdir)
