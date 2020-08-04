#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 05:42:37 2020

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
import os
import numpy as np
from os import listdir
from wpg.wavefront import Wavefront

indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/in/"
outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/out/"
newdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/tmp/"

fname = sys.argv[1]

def checkWavefront(fname):
    

    print("trying {}".format(fname))
    wfr = Wavefront()
    wfr.load_hdf5(outdir + fname)
    
    l = wfr.get_limits()
    np.save(newdir + fname, l)
# =============================================================================
#     try:
#         wfr.load_hdf5(outdir + fname)
#         del wfr
#     except RuntimeError:
#         print("removing {}".format(fname))
#         os.remove(outdir + fname)
# 
#     wfr = Wavefront()
#     try:
#         wfr.load_hdf5(outdir + fname)
#         del wfr
#     except OSError:
#         print("removing {}".format(fname))
#         os.remove(outdir + fname)
# =============================================================================


if __name__ == '__main__':
    checkWavefront(fname)            