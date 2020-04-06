#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 10:29:23 2020

@author: twguest
"""

import sys
sys.path.append("/opt/WPG/")

import scipy

import numpy as np

from wpg import srwlib
from wpg import srwlpy


from scipy.constants import c

from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import look_at_q_space, calculate_fwhm, check_sampling
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from wpg.wpg_uti_wf import plot_intensity_qmap as plotPhase

from wpg.generators import build_gauss_wavefront_xy as buildGaussian

from wpg.beamline import Beamline

from wpg.srwlib import SRWLOptL as thinLens

from wpg.optical_elements import Drift
from wpg.srwlib import SRWLOptD

fwhm2rms = np.sqrt(8*np.log(2)) ### FWHM = sqrt(8ln(2))*sigma

h = scipy.constants.physical_constants['Planck constant in eV s'][0]

def pulseEnergy(q, ekev):
    """
    Estimate of pulseEnergy from electron bunch charge and radiation energy    
    
    :param q: electron bunch charge [nC]
    :param ekev: radiation energy [keV] 
        
    :return P: pulse energy [J]
    """
    
    P = 19*q/ekev
    return P/1e3

def pulseDuration(q):
    """
    Estimate pulseDuration from electron bunch charge 
    
    :param q: electron bunch charge [nC]
        
    :return t: Duration of pulse [s]
    """
    
    t = (q*1e3)/9.8
    return t*1e-15


def pulseWidth(ekev):
    """
    Estimate pulseWidth (FWHM) from radiation energy (assumes symmetrical beam)
    
    :param ekev: radiation energy [keV]
        
    :return sig: Radiation pulse width (FWHM) [m]
    """
    
    sig = 6*np.log((7.4e03/ekev))
    return sig/1e6

def pulseDivergence(q, ekev):
    """
    Estimate of pulseDivergence (half-angle) from electron bunch charge and radiation energy    
    
    :param q: electron bunch charge [nC]
    :param ekev: radiation energy [keV] 
        
    :return dtheta: pulse divergence [rad]
    """
    
    dtheta = (17.2*6.4*np.sqrt(q))/ekev**(0.85)
    return dtheta/1e6

def tlFocus(sig, dtheta):
    """
    Thin lens focus required to correct divergence
    
    :params sig: beam fwhm [m]
    :params dtheta: pulse divergence [rad]
    
    :return f: thin lens focus [m]
    """
    f = sig/(2*np.tan(dtheta))
    return f



def coherentSource(nx, ny, nz, ekev, q, xoff = 0, yoff = 0):
    
    wavelength = (h*c)/(ekev*1e3)
        
    xMin, xMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    yMin, yMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    
    sigX, sigY = pulseWidth(ekev)/fwhm2rms, pulseWidth(ekev)/fwhm2rms
    pulseEn = pulseEnergy(q, ekev)
    
    dtheta = pulseDivergence(q, ekev)
    tau = pulseDuration(q)
    gsnBm = buildGaussian(nx, ny, ekev, xMin, xMax, yMin, yMax, sigX, sigY, 1)
    
     
 
    wfr = Wavefront(gsnBm)
    

    return wfr



if __name__ == '__main__':
    
    from matplotlib import pyplot as plt
    ## TEST FOR USAGE
    ### SET PARAMS
    
    q = 0.1 # nC
    ekev = 1 # keV
    
    ### Print Parameters as Sanity Check
    print("Electron Beam Charge: {} nC".format(q))    
    print("Electron Beam Charge: {} pC".format(q*1e3))   
    
    print("Radiation Energy: {} keV".format(ekev))   
    
    print("\n")
    ### Estimate Energy per Pulse
    print("Pulse Energy: {} Joules".format(pulseEnergy(q, ekev)))

    ### Estimate Duration of pulse
    print("Pulse Duration: {} seconds".format(pulseDuration(q)))
        
    ### Estimate FWHM of pulse (assumes symmetry)
    print("Pulse Width: {} m".format(pulseWidth(ekev)))
    
    ### Estimate FWHM of pulse (assumes symmetry)
    print("Pulse Divergence: {} rad".format(pulseDivergence(q,ekev)))
    
    ### Estimate thin lens focus for divergence correction
    print("Thin-Lens Focus: {} m".format(tlFocus(pulseWidth(ekev),pulseDivergence(q,ekev))))
    
    wfr = coherentSource(100, 100, 5, ekev, q)
    print(calculate_fwhm(wfr))
    plotIntensity(wfr)
    
    bl = Beamline()
 
    bl.append(Drift(20), [0,0,1,0,0,1,1,1,1,0,0,0])
    print(wfr.get_limits())
    bl.propagate(wfr)
    
    plotIntensity(wfr)
    print(wfr.get_limits())
    print(calculate_fwhm(wfr))
