#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 12:19:49 2020

@author: twguest
"""

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

from model.src.coherent import construct_SA1_wavefront
from model.beamline.structure import config

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from wpg.misc import fresnel_sampling
from wpg.srwlib import SRWLOptD
from model.beamline.structure import propagation_parameters
def quadratic_prop(outdir):
    
    wfr = construct_SA1_wavefront(1024,1024,6,1.0)
    D1 = wfr.params.Mesh.xMax - wfr.params.Mesh.xMin
    
    bl = config(focus = "micron")
    Drift = SRWLOptD(50)
    Drift.name = 'drift'
    bl.propagation_options[0]['optical_elements'] = [Drift]
    bl.propagation_options[0]['propagation_parameters'] = [bl.propagation_options[0]['propagation_parameters'][0]]
    
    bl.propagate_sequential(wfr, outdir = outdir + "quadratic/")
    D2 = wfr.params.Mesh.xMax - wfr.params.Mesh.xMin
    print(D1, D2)
    
    return wfr.pixelsize()[0], D1, D2

def fresnel_prop(outdir, pp):
    wfr = construct_SA1_wavefront(1024,1024,6,1.0)
    
    bl = config(focus = "micron")
    Drift = SRWLOptD(50)
    Drift.name = 'drift'
    bl.propagation_options[0]['optical_elements'] = [Drift]
    bl.propagation_options[0]['propagation_parameters'] = [pp]
    
    bl.propagate_sequential(wfr, outdir = outdir + "fresnel/")

    print(wfr.params.Mesh.xMax - wfr.params.Mesh.xMin)
if __name__ == '__main__':

    outdir = "../../output/sampling_test/d1/"
    wfr = construct_SA1_wavefront(1024,1024,6,1.0)

    print("Comparing Propagation Algorithms")
    dx1, D1, D2 = quadratic_prop(outdir)
    print(dx1)
    dx2 = fresnel_sampling(50, wfr.params.wavelength, dx1, D2, D1)
    
    z = D2/D1
    s = dx2/(dx1)
    print("zoom: {}".format(z))
    print("scale: {}".format(s))
    print(s/z)
    print(dx2)
    fresnel_prop(outdir, [0,0,1,0,0,z,s/z,z,s/z,0,0])
    
    #fresnel_prop(outdir)
