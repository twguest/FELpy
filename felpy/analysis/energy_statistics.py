#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:32:40 2020

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

import numpy as np

import wpg.srwlib as srwlib

from wpg.wpg_uti_wf import calc_pulse_energy, getOnAxisPowerDensity, getCentroid, get_profile_1d, getAxis
from wpg.wavefront import Wavefront

from felpy.model.src.coherent import coherentSource, construct_pulse
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity

from wpg.srwlib import srwl

def get_pulse_energy(wfr, mode = 'integrated'):
    """
    calculate pulse energy (in J), number of photons and fluence in the time-domain

    :param wfr: wpg.Wavefront structure
    :param mode: treat integrated or pulse wavefront intensity
    
    :return E: pulse energy statistics [(pulse energy, nphotons, fluence/mm^2), slices]
    """
    
    J2eV = 6.24150934e18
    
    if wfr.params.wDomain != "time":
        pass
    else:
        srwl.SetRepresElecField(wfr._srwl_wf, 't')
    
    x = getAxis(wfr, 'x')
    y = getAxis(wfr, 'y')
    t = getAxis(wfr, 't')
    
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    dt = t[1]-t[0]
    
    ax = (np.max(x)*1e-03) - (np.min(x)*1e-03)
    ay = (np.max(y)*1e-03) - (np.min(y)*1e-03)
    
    if mode == 'integrated':
        
        E = np.zeros([1,3])
        
        E[0,0] = wfr.get_intensity().sum()
        E[0,0] *=  dx * dy * 1e6 * dt
        E[0,1] = E[0,0] * J2eV / wfr.params.photonEnergy
        E[0,2] = E[0,1] / (ax*ay)
        E = E[0]
    elif mode == 'pulse':
        
        ii = wfr.get_intensity()
            
        E = np.zeros([ii.shape[-1],3])
            

        for slc in range(ii.shape[-1]):
            
            E[slc,0] = wfr.get_intensity()[:,:,slc].sum()
            E[slc,0] *=  dx * dy * 1e6 * dt
            E[slc,1] = E[slc,0] * J2eV / wfr.params.photonEnergy
            E[slc,2] = E[slc,1] / (ax*ay)
   
    return E



def run(wfr):
    ci = getCentroid(wfr, mode = 'integrated')
    cs = getCentroid(wfr, mode = 'pulse')
    Ei = getPulseEnergy(wfr, mode = 'integrated')
    Ep = getPulseEnergy(wfr, mode = 'pulse')
    
    return ci, cs, Ei, Ep

if __name__ == '__main__':
    
    wfr = construct_pulse(512,512,10)
    
    E = getPulseEnergy(wfr, mode = 'pulse')

    