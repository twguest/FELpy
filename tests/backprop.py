#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 16:05:57 2020

@author: twguest

test to see if backpropagation works (it does)
"""

import sys

sys.path.append("../")
sys.path.append("/opt/WPG/")

from model.tools import constructPulse
from model.beamline.structure import propagation_parameters
from wpg.beamline import Beamline
from wpg.optical_elements import Drift
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity



if __name__ == "__main__":
    
    wfr = constructPulse(nz = 2)
    plotIntensity(wfr)
    
    bl = Beamline()
    bl.append(Drift(-10), propagation_parameters(1,1,1,1))
    
    bl.propagate(wfr)
    
    plotIntensity(wfr)
    
