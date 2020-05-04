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
import os
from os.path import exists
import json
import numpy as np

from wpg import srwlpy
from wpg.misc import calcDivergence
from wpg.beamline import Beamline
from wpg import srwlib
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture
from wpg.srwlib import SRWLOptMirEl as MirEl2

from wpg.optical_elements import Mirror_elliptical as MirEl
from wpg.optical_elements import Screen


from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from model.src.coherent import coherentSource

from wpg.misc import get_wavelength 
from wpg.wpg_uti_wf import check_sampling, calculate_fwhm

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
    elif mode == "farfield":
        m = 2
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

def calcSampling(wfr, z, scale = 1, verbose = False):
    """
    returns the propagation parameters required to satisfy sampling requirements 
    for CDI 
    
    via: Spence, J.C.H., U. Weierstall, and M. Howells. 
    “Coherence and Sampling Requirements for Diffractive Imaging.” 
    Ultramicroscopy 101, no. 2–4 (November 2004): 149–52. 
    https://doi.org/10.1016/j.ultramic.2004.05.005.
    
    Considers both the x and y dimensions:
        pixelsize = (prop-distance*wavelength)/(2*object-width(4*sigma))
    
    :param wfr: WPG wavefront structure in the unpropagated plane
    :param z: propagation distance [m]
    :param scale: wavefield scaling factor over propagation distance
    
    :returns pp: propagation params expression req. to satisfy sampling [str]
    """
    
    fwhm2rms = np.sqrt(8*np.log(2)) ### FWHM = sqrt(8ln(2))*sigma

    fwhm_x, fwhm_y = calcDivergence(wfr) 
 
    sigma_x = fwhm_x/fwhm2rms
    sigma_y = fwhm_y/fwhm2rms
    
    #sigma_x = wfr.params.Mesh.xMax-wfr.params.Mesh.xMin
    delta_x = (z*wfr.params.wavelength)/(2*4*sigma_x)
    delta_y = (z*wfr.params.wavelength)/(2*4*sigma_y)

    zoom_x = delta_x/(wfr.pixelsize()[0]*scale)
    zoom_y = delta_y/(wfr.pixelsize()[1]*scale)
    
    pp = propParams(scale, zoom_x, scale, zoom_y)
    
    if verbose == True:
        print("Pixel Size in Unpropagated Plane: {} x {}".format(wfr.pixelsize()[0], wfr.pixelsize()[1]))
        print("Pixel Size in Propagated Plane: {} x {}".format(delta_x, delta_y))
    return pp
    
    
def load_params():
    with open("../../data/input/parameters.json", "r") as read_file:
        params = json.load(read_file)
    return params

def buildBeamline(params, focus = "micron"):
        ### generate hom1 mirror surface
    if exists("../../data/hom1_mir_Flat.dat"):
        hom1_profile = "../../data/hom1_mir_Flat.dat"
    else:
        genMirrorSurface(500, 500, [10, 10], "../../data/hom1_", mode = 'Flat')
        hom1_profile = "../../data/hom1_mir_Flat.dat"
    ### generate hom1 mirror surface
    if exists("../../data/hom2_mir_Flat.dat"):
        hom2_profile = "../../data/hom2_mir_Flat.dat"
    else:
        genMirrorSurface(500, 500, [10, 10], "../../data/hom2_", mode = 'Flat')
        hom2_profile = "../../data/hom2_mir_Flat.dat"
    if exists("../../data/mhp_mir_Flat.dat"):
        mhp_profile = "../../data/mhp_mir_Flat.dat"
    else:
        genMirrorSurface(500, 500, [10, 10], "../../data/mhp_", mode = 'Flat')
        mhp_profile = "../../data/mhp_mir_Flat.dat"
    if exists("../../data/mvp_mir_Flat.dat"):
        mvp_profile = "../../data/mvp_mir_Flat.dat"
    else:
       genMirrorSurface(500, 500, [10, 100], "../../data/mvp_", mode = 'Flat')
       mvp_profile = "../../data/mvp_mir_Flat.dat"
       
    d1 =  Drift(params["d1"]['distance'])
    d1.name = params["d1"]['name']
    
    HOM1 = MirPl(np.loadtxt(params['HOM1']['mirror profile']),
                 _dim = params['HOM1']['orientation'],
                 _ang = params['HOM1']['incidence angle'], 
                 _amp_coef = params['HOM1']['transmission'],
                 _x = params['HOM1']['xc'], _y = params['HOM1']['yc']) 
    HOM1.name = params['HOM1']['name']
    
    d2 =  Drift(params["d2"]['distance'])
    d2.name = params["d2"]['name']
    
    HOM2 = MirPl(np.loadtxt(params['HOM2']['mirror profile']),
                 _dim = params['HOM2']['orientation'],
                 _ang = params['HOM2']['incidence angle'], 
                 _amp_coef = params['HOM2']['transmission'],
                 _x = params['HOM2']['xc'], _y = params['HOM2']['yc']) 
    HOM2.name = params['HOM2']['name']
    
    d3 =  Drift(params["d3"]['distance'])
    d3.name = params["d3"]['name']
    
    MKB_pslit = Aperture(_shape=params["MKB_pslit"]['shape'],
                         _ap_or_ob=params["MKB_pslit"]['type'],
                         _Dx= params["MKB_pslit"]['dx'],
                         _Dy= params["MKB_pslit"]['dy'],
                         _x=params["MKB_pslit"]['xc'],
                         _y=params["MKB_pslit"]['yc'])
    MKB_pslit.name = params["MKB_pslit"]['name']
    
    d4 =  Drift(params["d4"]['distance'])
    d4.name = params["d4"]['name']
    
    MHE_ap = Aperture(_shape=params["MHE_ap"]['shape'],
                         _ap_or_ob=params["MHE_ap"]['type'],
                         _Dx= params["MHE_ap"]['dx'],
                         _Dy= params["MHE_ap"]['dy'],
                         _x=params["MHE_ap"]['xc'],
                         _y=params["MHE_ap"]['yc'])
    MHE_ap.name = params["MHE_ap"]['name']
    
    MHP = MirPl(np.loadtxt(params['MHP']['mirror profile']),
                 _dim = params['MHP']['orientation'],
                 _ang = params['MHP']['incidence angle'], 
                 _amp_coef = params['MHP']['transmission'],
                 _x = params['MHP']['xc'], _y = params['MHP']['yc']) 
    MHP.name = params['MHP']['name']
    
    d5 =  Drift(params["d5"]['distance'])
    d5.name = params["d5"]['name']
    
    ###TEST
# =============================================================================
#     MHE = MirEl2(_p=params['MHE']["distance from source"],
#                 _q=params['MHE']["distance to focus"],
#                 _ang_graz= 1,#params['MHE']["incident angle"],
#                 _r_sag=10e3,
#                 _size_tang=0.950,#width
#                 _size_sag=0.026,
#                 _ap_shape="r",
#                 _sim_meth=1,
#                 _npt=100,
#                 _nps=100,
#                 _treat_in_out=0,
#                 _ext_in=0.5,#params['MHE']["_ext_in"],
#                 _ext_out=0.5,#params['MHE']["_ext_out"],
#                 _nvx=np.cos(params['MHE']["incident angle"]),
#                 _nvy=np.sin(params['MHE']["roll"]),
#                 _nvz=-np.sin(params['MHE']["incident angle"]),
#                 _tvx=-np.sin(params['MHE']["incident angle"]),
#                 _tvy=0,
#                 _x=0,#np.tan(params['MHE']["yaw"])*params['MHE']["displacement"],
#                 _y=0,
#                 _refl=1,
#                 _n_ph_en=1,
#                 _n_ang=1,
#                 _n_comp=1,
#                 _ph_en_start=1000.0,
#                 _ph_en_fin=1000.0,
#                 _ph_en_scale_type="lin",
#                 _ang_start=0,
#                 _ang_fin=0,
#                 _ang_scale_type="lin")
#     ###
#     
# =============================================================================
    MHE = MirEl(orient = params['MHE']["orientation"], p = params['MHE']["distance from source"], q = params['MHE']["distance to focus"],
                thetaE = 0.005, theta0 = 0.005+2.2e-06,
                _x = 0, _y = 0,#params['MHE']["displacement"],
                length = params['MHE']["length"],
                roll = params['MHE']["roll"],
                yaw = params['MHE']["yaw"],
                _refl = params['MHE']["reflectivity"],
                _ext_in = params['MHE']["_ext_in"], _ext_out = params['MHE']["_ext_out"]) 
    
    MHE.name = "MHE"
    
    d6 =  Drift(params["d6"]['distance'])
    d6.name = params["d6"]['name']

    d7 =  Drift(params["d7"]['distance'])
    d7.name = params["d7"]['name']
    
    MVE_ap = Aperture(_shape=params["MVE_ap"]['shape'],
                         _ap_or_ob=params["MVE_ap"]['type'],
                         _Dx= params["MVE_ap"]['dx'],
                         _Dy= params["MVE_ap"]['dy'],
                         _x=params["MVE_ap"]['xc'],
                         _y=params["MVE_ap"]['yc'])
    MVE_ap.name = params["MVE_ap"]['name']
    
# =============================================================================
#     MVE = MirEl(orient = params['MVE']["orientation"],
#                 p = params['MVE']["distance from source"],
#                 q = params['MVE']["distance to focus"],
#                 thetaE = params['MVE']["design angle"], theta0 = params['MVE']["incident angle"],
#                 distance = params['MVE']["displacement"], length = params['MVE']["length"],
#                 roll = params['MVE']["roll"],
#                 yaw = params['MVE']["yaw"],
#                 _refl = params['MVE']["reflectivity"],
#                 _ext_in = params['MVE']["_ext_in"], _ext_out = params['MVE']["_ext_out"]) 
#     
#     MVE.name = "MVE"
#     
#    
# =============================================================================
    d8 =  Drift(params["d8"]['distance'])
    d8.name = params["d8"]['name']
    
    
    
    MVP = MirPl(np.loadtxt(params['MVP']['mirror profile']),
                 _dim = params['MVP']['orientation'],
                 _ang = params['MVP']['incidence angle'], 
                 _amp_coef = params['MVP']['transmission'],
                 _x = params['MVP']['xc'], _y = params['MVP']['yc']) 
    MVP.name = params['MVP']['name']
    
    df =  Drift(params["df"]['distance'])
    df.name = params["df"]['name']
    
    
    
    if focus == "micron":
        bl = Beamline()
        bl.append(d1, propParams(1,1,1,1, mode = "farfield"))

        #bl.append(HOM1, propParams(1, 1, 1, 1, mode = 'normal'))
        #bl.append(d2, propParams(1, 1, 1, 1, mode = 'quadratic'))
        #bl.append(HOM2,  propParams(1, 1, 1, 1, mode = 'normal'))
        #bl.append(d3, propParams(1,1,1,1, mode = 'farfield'))
        
        #bl.append(MKB_pslit, propParams(1/5, 1, 1/5, 1, mode = 'normal'))
        #bl.append(d4, propParams(1, 1, 1, 1, mode = 'quadratic'))
        #bl.append(MHP, propParams(1, 1, 1, 1, mode = 'normal'))
        #bl.append(d5, propParams(1, 1, 1, 1, mode = 'quadratic'))
        #bl.append(MHE_ap, propParams(1, 1, 1, 1, mode = 'normal'))
        bl.append(MHE, propParams(1, 1, 1, 1, mode = 'normal'))
        #bl.append(d6, propParams(1, 1, 1, 1, mode = 'quadratic'))
     
        #bl.append(d7, propParams(1, 1, 1, 1, mode = 'quadratic'))
        #bl.append(MVE_ap, propParams(1, 1, 1, 1, mode = 'normal'))
        #bl.append(MVE, propParams(1, 1, 1, 1, mode = 'normal'))
        #bl.append(d8, propParams(1, 1, 1, 1, mode = 'quadratic'))
    
        #bl.append(MVP, propParams(1, 1, 1, 1, mode = 'normal'))
        #bl.append(df, propParams(1, 10, 1, 10, mode = 'converge'))
    

    bl.params = params
    return bl

 
def config(focus = "micron", screens = True):
    
    params = load_params()
    bl = buildBeamline(params)
    
    return bl

    
def testPropEnergies(outdir):
    
    erange = [6000, 9200, 12000, 16000]
    
    if os.path.exists(outdir):
        pass
    else: 
        os.mkdir(outdir)
    
    for energy in erange:
        if os.path.exists(outdir + "/{}eV_1nC/".format(energy).replace(".","_")) == False:
            os.mkdir(outdir + "/{}eV_1nC/".format(energy).replace(".","_"))
    
    for energy in erange:
        wfr = coherentSource(1048, 1048, energy/1000, 1)
        
        bl = config()
        bl.propagateSeq(wfr, outdir + "/{}eV_1nC/".format(energy).replace(".","_"))
        
if __name__ == '__main__':
    #testPropEnergies("../../output")

    #params = load_params()
    wfr = coherentSource(1048, 1048, 16, 1)

    plotIntensity(wfr)
    bl = config()
    bl.propagateSeq(wfr)
    plotIntensity(wfr)
    #test_mirror_setup()
