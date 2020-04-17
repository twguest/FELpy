#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:41:47 2020

@author: twguest
"""
#################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################

from os.path import exists

import numpy as np

from wpg import srwlpy
 
from wpg.beamline import Beamline
from wpg import srwlib
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture

from wpg.optical_elements import Mirror_elliptical as MirEl
from wpg.optical_elements import Screen


from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from model.src.coherent import coherentSource
from model.beamline.structure import calcSampling, propParams
from wpg.misc import get_wavelength, sampling
from wpg.wpg_uti_wf import check_sampling, calculate_fwhm

from wpg.srwlib import srwl_opt_setup_surf_height_2d as MirPl

#SCRATCH PAD FOR SAMPLING CALC

if __name__ == '__main__':
    
   
    ### generate hom1 mirror surface
    if exists("../../data/hom1_mir_Flat.dat"):
        hom1_profile = "../../data/hom1_mir_Flat.dat"
    else:
        genMirrorSurface(200, 200, [1000e-03, 25e-03], "../../data/hom1_", mode = 'Flat')
     
    ### generate hom1 mirror surface
    if exists("../../data/hom2_mir_Flat.dat"):
        hom2_profile = "../../data/hom2_mir_Flat.dat"
    else:
        genMirrorSurface(200, 200, [1000e-03, 25e-03], "../../data/hom2_", mode = 'Flat')
     
    if exists("../../data/mhp_mir_Flat.dat"):
        mhp_profile = "../../data/mhp_mir_Flat.dat"
    else:
        genMirrorSurface(200, 200, [1000e-03, 25e-03], "../../data/mhp_", mode = 'Flat')
        
    if exists("../../data/mvp_mir_Flat.dat"):
        mvp_profile = "../../data/mvp_mir_Flat.dat"
    else:
        genMirrorSurface(200, 200, [25e-03, 1000e-03], "../../data/mvp_", mode = 'Flat')
    
    d1 =  Drift(246.5)
    
    HOM1 = MirPl(srwlib.srwl_uti_read_data_cols(hom1_profile, "\t"),
                 _dim = 'x',
                 _ang = 2.1e-03, 
                 _amp_coef = 1,
                 _x = 0, _y = 0) 
    
    d2 = Drift(11.36)
    
    HOM2 = MirPl(srwlib.srwl_uti_read_data_cols(hom2_profile, "\t"),
                 _dim = 'x',
                 _ang = 2.4e-03, 
                 _amp_coef = 1,
                 _x = 0, _y = 0) 
    
    d3 = Drift(634.669)
        
    MKB_pslit = Aperture(_shape="r", _ap_or_ob="a", _Dx= 0.100, _Dy= 0.100, _x=0, _y=0)
    
    d4 = Drift(1.200)
    
    ap_MHE = Aperture(_shape="r", _ap_or_ob="a", _Dx= 0.950, _Dy= 0.025, _x=0, _y=0)
    
    MHP = MirPl(srwlib.srwl_uti_read_data_cols(mhp_profile, "\t"),
                _dim = 'x',
                _ang = 1.1e-03,
                _amp_coef = 1,
                _x = 0,
                _y = 0)
    
    
    d5 =  Drift(1.050)
    
    
    MHE = MirEl(orient = 'x', p = 894.779, q = 23.905, thetaE = 0, theta0 = 0,
                length= 1.000, roll = 0, yaw = 0, _x = 0, _y = 0, _refl = 1,
                _ext_in = 0.500, _ext_out = 0.500) 
    
    
    
    ap_MVE = Aperture(_shape="r", _ap_or_ob="a", _Dx= 0.025, _Dy= 0.950, _x=0, _y=0)
    
    d6 = Drift(0.36)
    
    MKB_scr = Screen()
    
    d7 = Drift(1.320)
    
    MVE = MirEl(orient = 'y', p = 896.459, q = 22.225, thetaE = 0, theta0 = 0,
                length= 1.000, roll = 0, yaw = 0, _x = 0, _y = 0, _refl = 1,
                _ext_in = 0.500, _ext_out = 0.500) 
    
    d8 = Drift(1.050)
    
    MVP = MirPl(srwlib.srwl_uti_read_data_cols(mvp_profile, "\t"),
                _dim = 'y',
                _ang = 1.1e-03,
                _amp_coef = 1,
                _x = 0,
                _y = 0)
    
    df = Drift(21.175)
    
    ##########################################################################
    ##########################################################################
    
    # SAMPLING TEST BEGINS HERE
    
    ##########################################################################
    ##########################################################################
    
    # INIT WFR.
    wfr = coherentSource(1048, 1048, 6, 1.0)
    
    print("D1")
    bl = Beamline()
    bl.append(d1, calcSampling(wfr, d1.L, scale = 50, verbose = True))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/D1_6keV.hdf5")
    
    print("HOM1")
    bl = Beamline()
    bl.append(HOM1, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/HOM1_6keV.hdf5")

    print("D2")
    bl = Beamline()
    bl.append(d2, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/D2_6keV.hdf5")
    
    print("HOM2")
    bl = Beamline()
    bl.append(HOM2, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/HOM2_6keV.hdf5")
    
    print("D3")
    bl = Beamline()
    bl.append(d3, propParams(2,1,2,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/D3_6keV.hdf5")
    
    print("MKB_PSLIT")
    bl = Beamline()
    bl.append(MKB_pslit, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/MKBSlit_6keV.hdf5")
    
    print("D4")
    bl = Beamline()
    bl.append(d4, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/D4_6keV.hdf5")

    print("MHP")
    bl = Beamline()
    bl.append(MHP, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/MHP_6keV.hdf5")

    print("D5")
    bl = Beamline()
    bl.append(d5, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/D5_6keV.hdf5")
    
    print("ap_MHE")
    bl = Beamline()
    bl.append(ap_MHE, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/apMHE_6keV.hdf5")
    
    print("MHE")
    bl = Beamline()
    bl.append(MHE, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/MHE_6keV.hdf5")
    
    print("D6")
    bl = Beamline()
    bl.append(d6, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/D6_6keV.hdf5")
    
    
        
    print("MKB_Screen")
    bl = Beamline()
    bl.append(MKB_scr, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/MKBSCRN_6keV.hdf5")
    
    print("D7")
    bl = Beamline()
    bl.append(d7, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/D7_6keV.hdf5")
    
    print("ap_MVE")
    bl = Beamline()
    bl.append(ap_MVE, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/apMVE_6keV.hdf5")
    
    
    print("MVE")
    bl = Beamline()
    bl.append(MVE, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/MVE_6keV.hdf5")
    
    print("D8")
    bl = Beamline()
    bl.append(d8, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/D8_6keV.hdf5")

    print("MVP")
    bl = Beamline()
    bl.append(MVP, propParams(1,1,1,1))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/MVP_6keV.hdf5")
    
    print("DFocus")
    bl = Beamline()
    bl.append(d8, propParams(1/50,1,1/50,1, m ='converge'))
    bl.propagate(wfr)
    plotIntensity(wfr)
    wfr.store_hdf5("../../data/DFoc_6keV.hdf5")
