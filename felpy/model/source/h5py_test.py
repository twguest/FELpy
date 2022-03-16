#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.0"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
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

from felpy.model.wavefront import Wavefront
from wpg.wpg_uti_wf import animate, plot_axial_power_density, plot_total_power

filename = "../../data/test_h5/XFEL_S1_04.96keV_12.0GeV_0020pC_SASE_U_BLI_2014-05-01_FAST_FXY1_0002020.h5"

wfr = Wavefront()
wfr.load_hdf5(filename)



#animate(wfr, delay = 2, outdir = "../../tmp")

plot_axial_power_density(wfr, spectrum = False, outdir = "../../tmp")
plot_axial_power_density(wfr, spectrum = True, outdir = "../../tmp")
plot_total_power(wfr, spectrum = False, outdir = "../../tmp")
plot_total_power(wfr, spectrum = True, outdir = "../../tmp")
print(integral_intensity(wfr))