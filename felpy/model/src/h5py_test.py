#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "Apache"
__version__ = "1.0.0"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

 
import os, shutil
from os.path import exists
import json
import numpy as np
import h5py

from matplotlib import pyplot as plt
import matplotlib as mpl 
from wpg.wpg_uti_wf import integral_intensity
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
print(integral_intensity(wfr))