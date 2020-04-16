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

from wpg.misc import get_wavelength, sampling
from wpg.wpg_uti_wf import check_sampling

from wpg.srwlib import srwl_opt_setup_surf_height_2d as MirPl

def propParams(sx, zx, sy, zy, mode = "normal"):
    """
    wrapper for propagation parameters
    
    :param zx: horizontal scaling factor 
    :param rx: horizontal zoom factor
    :param zy: vertical scaling factor
    :param ry: vertical zoom factor
    :param mode: normal, semi-analytical, converge or diverge
    
    :return propagation parameters:
    """
    
    if mode == "normal":
        m = 0
    elif mode == "quadratic":
        m = 1
    elif mode == "diverge":
        m = 3
    elif mode == "converge":
        m = 4
    
    return [0,0,1,m,0,sx,zx/sx,sy,zy/sy,0,0,0]


def genMirrorSurface(nx, ny, mirDim, outdir, mode = 'Flat'):
    
    if mode == 'Flat':
        
        mirLen = mirDim[0]
        mirWid = mirDim[1]
        
        surface = np.zeros((nx,ny))
        
        surface[1:, 0] = np.linspace(-mirLen/2, mirLen/2, nx-1) 
        surface[0, 1:] = np.linspace(-mirWid/2, mirWid/2, ny-1) 
        
    np.savetxt(outdir+"mir_"+ mode +".dat", surface, delimiter='\t')
        
    ## genMirrorSurface(200, 200, [1000e-03, 80e-03], "/opt/spb_model/data/", mode = 'Flat')

def NyquistFreq(wfr, z, scaler, axisName = 'x'):
    """
    calculate the required zoom to satisfy nyquist-shannan sampling
    conditions in the propagated plane
    
    for calculation of req. pixel size see: wpg/misc/sampling 
    
    :param wavelength: source wavelength [m]
    :param z: propagation distance [m]
    :param scaler: field scale multiplier (ie., extent_f = scale*extent_i)
    
    :return zoom: zoom modifier for propagation step
    """
    wavelength = get_wavelength(wfr)

    if axisName == 'x':
        d = wfr.pixelsize()[0]
        extent = wfr.params.Mesh.xMax - wfr.params.Mesh.xMin
    elif axisName == 'y':
        d = wfr.pixelsize()[1]
        extent = wfr.params.Mesh.yMax - wfr.params.Mesh.yMin

    W = (wavelength*z)/(extent*scaler*d*np.sqrt(2))
    
    return W/d
    
    
def config(beamline = "micro", screens = True):
    
   
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
    
    bl = Beamline()
    bl.append(d1, propParams(30,4,30, 4))
    #bl.append(HOM1, propParams(1, 1, 1, 1))
    
    return bl
 

if __name__ == '__main__':
    
    
    wfr = coherentSource(1048, 1048, 3, 1.0)
    plotIntensity(wfr)
    
    bl = config()
    bl.propagate(wfr)
    print(check_sampling(wfr))
    plotIntensity(wfr)