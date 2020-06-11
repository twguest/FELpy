#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:43:15 2020

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
import os, shutil
from os.path import exists
import json
import numpy as np
import h5py

from matplotlib import pyplot as plt
import matplotlib as mpl 

import h5py

from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import animate, plotOnAxisPowerDensity, plotTotalPower

filename = "../../data/test_h5/XFEL_S1_04.96keV_12.0GeV_0020pC_SASE_U_BLI_2014-05-01_FAST_FXY1_0002020.h5"

wfr = Wavefront()
wfr.load_hdf5(filename)



#animate(wfr, delay = 2, outdir = "../../tmp")

plotOnAxisPowerDensity(wfr, spectrum = False, outdir = "../../tmp")
plotOnAxisPowerDensity(wfr, spectrum = True, outdir = "../../tmp")
plotTotalPower(wfr, spectrum = False, outdir = "../../tmp")
plotTotalPower(wfr, spectrum = True, outdir = "../../tmp")
