#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 01:54:03 2020

@author: guestt
"""

###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################

from os import listdir
from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import plot_intensity_map as plotIntentisty

def main():
    
    indir = "../../data/h5/NanoKB-Pulse/out/"
    
    for fname in listdir(indir):
        wfr = Wavefront()
        wfr.load_hdf5(indir + fname)
        plotIntentisty(wfr, save = indir  +"/int/" + fname.replace("h5", "png"))
        print("fname")
        del wfr
if __name__ == "__main__":
    main()