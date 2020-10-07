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
import matplotlib.patches as mpatches
from wpg import srwlpy
 
from wpg.beamline import Beamline
from wpg import srwlib
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture

from wpg.optical_elements import Mirror_elliptical as MirEl
from wpg.optical_elements import Screen


from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity 
from model.src.coherent import coherentSource
from wpg.wpg_uti_wf import check_sampling, calculate_fwhm
from wpg.misc import calcDivergence
from wpg.srwlib import srwl_opt_setup_surf_height_2d as MirPl
from matplotlib import pyplot as plt
 


def fwhm2rms(fwhm, sigma = 1):
    return (fwhm*sigma)/np.sqrt(8*np.log(2))

def getRequiredDistance(wfr, dp, axisName = 'x', buffer = 1):
    
    """
    Calculate the propagation distance required to satisfy sampling conditions.
    
    :param wfr: wpg wfr structure in unpropagated plane
    :param z: propagation distance
    :param dp: propagated plane (detector) pixel size
    :param axisName: str transverse axis
    
    :returns S: Oversampling index (S > 1 is oversampled)
    """

    if type(dp) == list:
        dp = max(dp)
        
    Wx = fwhm2rms(calculate_fwhm(wfr)['fwhm_x'], sigma = 4)*buffer
    Wy = fwhm2rms(calculate_fwhm(wfr)['fwhm_y'], sigma = 4)*buffer
    
    W = max(Wx, Wy)
    
    zreq = (W*dp)/(wfr.params.wavelength)
    print("required propagation distance: {} m".format(zreq))
    return zreq

def getImageProperties(wfr, z, d_pixelsize, d_npix):
    """
    Return degree to which beam FWHM fills detector
    
    :param wfr: wpg wfr structure in unpropagated plane
    :param z: propagation distance
    :param d_pixelsize: detector pixelsize
    :param d_npix: number of pixels at detector
    
    :returns O: overfilling factor (O > 1 is overfilled)
    """
    
    if type(d_pixelsize) == list:
        d_px = d_pixelsize[0]
        d_py = d_pixelsize[1]
    else:
        d_px = d_py = d_pixelsize
        
    if type(d_npix) == list:
        d_nx = d_npix[0]
        d_ny = d_npix[1]
    else:
        d_nx = d_ny = d_npix
        
    
    ## get fwhm and divergence of 
    Wx = fwhm2rms(calculate_fwhm(wfr)['fwhm_x'], sigma = 4)
    div_x = calcDivergence(wfr)[0]

    Wy = fwhm2rms(calculate_fwhm(wfr)['fwhm_y'], sigma = 4)
    div_y = calcDivergence(wfr)[1]
    
    ## calculate magnified beam size
    mW_x = Wx + z*np.tan(div_x)
    mW_y = Wy + z*np.tan(div_y)
    
    ## calculate detector size
    dx = d_nx*d_px
    dy = d_ny*d_py
    
    ## calculate the percentage of the detector thats filled in each direction
    fill_x = mW_x/dx
    fill_y = mW_y/dy
    
    print("% of Detector Filled By Beam in X-direction: {} %".format(fill_x*100))
    print("% of Detector Filled By Beam in Y-direction: {} %".format(fill_y*100))

    ## calculate the number of pixels the beam occupies in the detector plane
    bx = mW_x/d_px
    by = mW_y/d_py
    
    print("Number of Pixels Beam Occupies in X-direction: {} pixels".format(bx))
    print("Number of Pixels Beam Occupies in Y-direction: {} pixels".format(by))
    
    ## calculate effective pixel-resolution
    ## ie., what area of the focus beam is represented per pixel
    res_x = Wx/bx
    res_y = Wy/by
    
    print("Image Pixel Resolution in X-direction: {} um".format(res_x*1e6))
    print("Image Pixel Resolution in Y-direction: {} um".format(res_y*1e6))

 

if __name__ == '__main__':
    
    wfr = coherentSource(250, 250, 9.2, 3)    
    
    d_pixelsize = 13.5e-06
    d_npix = [2456,2058]
    
    zreq = getRequiredDistance(wfr, d_pixelsize, buffer = 1)
    
    getImageProperties(wfr, zreq, d_pixelsize, d_npix)
    
    
    