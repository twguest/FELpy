#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 10:31:29 2020

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

import multiprocessing

from model.beamline.structure import propParams
from model.beamline.structure import BeamlineModel
from model.src.coherent import coherentSource
from wpg import srwlib

from wpg.srwlib import SRWLOptD as Drift

from wpg.wavefront import Wavefront
from wpg.beamline import Beamline

from wpg.wpg_uti_wf import calc_pulse_energy, calculate_fwhm, getOnAxisPowerDensity, getCentroid
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity

from wpg.misc import calcDivergence

from os import listdir

from tqdm import tqdm

if __name__ == '__main__':
    
    wfr = coherentSource(1124, 1423, 4.96, 0.25)

    centroid = getCentroid(wfr)