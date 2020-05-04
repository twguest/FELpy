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
from model.beamline.structure import load_params, propParams

from wpg.beamline import Beamline
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from wpg.srwlib import SRWLOptMirEl as MirEl
from wpg.srwlib import SRWLOptD
if __name__ == '__main__':

    outdir = "../../output/"
    #outdir = None
    
    params = load_params()
    
    print("Testing Beamline Propagation")
    
    MHE = MirEl(_p=params['MHE']["distance from source"],
            _q=params['MHE']["distance to focus"],
            _ang_graz=params['MHE']["incident angle"],
            _r_sag=10e3,
            _size_tang=1,
            _size_sag=0.025,
            _ap_shape="r",
            _sim_meth=2,
            _npt=1,
            _nps=1,
            _treat_in_out=0,
            _ext_in=0.5,#params['MHE']["_ext_in"],
            _ext_out=0.5,#params['MHE']["_ext_out"],
            _nvx=np.cos(params['MHE']["incident angle"]),
            _nvy=np.sin(params['MHE']["roll"]),
            _nvz=-np.sin(params['MHE']["incident angle"]),
            _tvx=-np.sin(params['MHE']["incident angle"]),
            _tvy=0,
            _x= 0,#np.tan(params['MHE']["yaw"])*params['MHE']["displacement"],
            _y=0,
            _refl=1,
            _n_ph_en=1,
            _n_ang=1,
            _n_comp=1,
            _ph_en_start=1000.0,
            _ph_en_fin=1000.0,
            _ph_en_scale_type="lin",
            _ang_start=0,
            _ang_fin=0,
            _ang_scale_type="lin")
    MHE.name = 'horizontal ellipse'
    
    wfr = coherentSource(1024,1024,6,1.0)
    
    Drift = SRWLOptD(5)
    Drift.name = 'drift'
    bl = Beamline()
    bl.append(MHE, propParams(1,1,1,1))
    #bl.append(Drift, propParams(1,1, 1,1))
    bl.propagateSeq(wfr)
        