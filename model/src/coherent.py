#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 10:29:23 2020

@author: twguest
"""

import sys
sys.path.append("/opt/WPG/")

import numpy as np

from wpg import srwlib
from wpg import srwlpy

from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import look_at_q_space, plot_intensity_map, calculate_fwhm
from wpg.beamline import Beamline

from wpg.srwlib import SRWLOptL as thinLens

from wpg.optical_elements import Drift

fwhm2rms = np.sqrt(8*np.log(2)) ### FWHM = sqrt(8ln(2))*sigma



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

def moddedGaussian(nx, ny, nz, ekev, xMin, xMax, yMin, yMax, tau, sigX, sigY, 
                   d2waist, xoff = 0, yoff = 0, xtilt = 0, ytilt = 0,
                   pulseEn=None, pulseRange=None, _mx=None, _my=None):
    """
    Build modded 3D Gaussian beam.

    :param nx: Number of point along x-axis
    :param ny: Number of point along y-axis
    :param nz: Number of point along z-axis (slices)
    :param ekev: Energy in kEv
    :param xMin: Initial Horizontal Position [m]
    :param xMax: Final Horizontal Position [m]
    :param yMin: Initial Vertical Position [m]
    :param yMax: Final Vertical Position [m]
    :param tau: Pulse duration [s]
    :param sigX: Horiz. RMS size at Waist [m]
    :param sigY:  Vert. RMS size at Waist [m]
    :param xoff: horizontal position of gaussian beam centre [m]
    :param yoff: vertical position of gaussian beam centre [m]
    :param xtilt: average angle of beam at waist [rad]
    :param ytilt: average angle of beam at waist [rad]
    :param d2waist: Distance to Gaussian waist
    :param pulseEn: Energy per Pulse [J]
    :param pulseRange: pulse duration range in sigT's
    :param _mx: transverse Gauss-Hermite mode order in horizontal direction
    :param _my: transverse Gauss-Hermite mode order in vertical direction
    
    :return: wpg.Wavefront structure
    """


    GsnBm = srwlib.SRWLGsnBm()  # Gaussian Beam structure (just parameters)
    GsnBm.x = xoff  # Transverse Coordinates of Gaussian Beam Center at Waist [m]
    GsnBm.y = yoff
    GsnBm.z = 0  # Longitudinal Coordinate of Waist [m]
    GsnBm.xp = xtilt  # Average Angles of Gaussian Beam at Waist [rad]
    GsnBm.yp = ytilt

    GsnBm.avgPhotEn = ekev * 1.e3  # 15000. #Photon Energy [eV]
    if pulseEn is not None:
        GsnBm.pulseEn = pulseEn 
    else:
        GsnBm.pulseEn = 0.001 # was 1 mJ in the Tutorial exampes as well
    GsnBm.repRate = 1  # Rep. Rate [Hz] - to be corrected
    GsnBm.polar = 1  # 1- linear hoirizontal; 2 - linear vertical
    # Far field angular divergence: 14.1e-6 ./ (ekev) .^0.75
    # 0.17712e-09/(4*Pi)/(14.1e-06/((7)^0.75)) for 7 keV, 3.55561e-06 =
    # 0.0826561e-09/(4*Pi)/(14.1e-06/((15)^0.75)) for 15 keV #Horiz. RMS size
    # at Waist [m]
    GsnBm.sigX = sigX
    GsnBm.sigY = sigY  # Vert. RMS size at Waist [m]
    # Coherence time (~ Gaussian pulse duration)           0.12 fs @ 15 keV
    # and 0.17 fs @ 7 keV
    # 0.12e-15 #0.17e-15 for 15 keV #Pulse duration [s] #To check: Is it 0.12
    # fs or 12 fs ?
    GsnBm.sigT = tau
    if _mx is not None:
        GsnBm.mx = _mx 
    else:
        GsnBm.mx = 0  # Transverse Gauss-Hermite Mode Orders
    if _mx is not None:
        GsnBm.my = _my 
    else:
        GsnBm.my = 0

    wfr = srwlib.SRWLWfr()  # Initial Electric Field Wavefront
    wfr.allocate(nz, nx, ny)
    # Numbers of points vs Photon Energy (1), Horizontal and
    # Vertical Positions (dummy)
    wfr.presFT = 1  # Defining Initial Wavefront in Time Domain
    #wfr.presFT = 0 #Defining Initial Wavefront in Frequency Domain

    wfr.avgPhotEn = GsnBm.avgPhotEn
    if pulseRange is not None:
        wfr.mesh.eStart = -pulseRange/2. * GsnBm.sigT  # Initial Time [s]
        wfr.mesh.eFin   =  pulseRange/2. * GsnBm.sigT  # Final Time [s]
    else:
        #wfr.mesh.eStart = -100 * GsnBm.sigT  # Initial Time [s]
        #wfr.mesh.eFin = 100 * GsnBm.sigT  # Final Time [s]
        wfr.mesh.eStart = -4. * GsnBm.sigT  # Initial Time [s]
        wfr.mesh.eFin = 4. * GsnBm.sigT  # Final Time [s]

    # Longitudinal Position [m] at which Electric Field has to be calculated,
    # i.e. the position of the first optical element
    wfr.mesh.zStart = d2waist
    wfr.mesh.xStart = xMin  # Initial Horizontal Position [m]
    wfr.mesh.xFin = xMax  # Final Horizontal Position [m]
    wfr.mesh.yStart = yMin  # Initial Vertical Position [m]
    wfr.mesh.yFin = yMax  # Final Vertical Position [m]

    wfr.mesh.ne = nz

    # Some information about the source in the Wavefront structure
    wfr.partBeam.partStatMom1.x = GsnBm.x
    wfr.partBeam.partStatMom1.y = GsnBm.y
    wfr.partBeam.partStatMom1.z = GsnBm.z
    wfr.partBeam.partStatMom1.xp = GsnBm.xp
    wfr.partBeam.partStatMom1.yp = GsnBm.yp

    sampFactNxNyForProp = -1  # 5 #sampling factor for adjusting nx, ny (effective if > 0)
    arPrecPar = [sampFactNxNyForProp]
    #**********************Calculating Initial Wavefront
    srwlpy.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)

    return wfr

def coherentSource(nx, ny, nz, ekev, q, xoff = 0, yoff = 0):
    
    xMin, xMax = -200e-06, 200e-06 #based on fwhm of 1 nC, 3 keV beam
    yMin, yMax = -200e-06, 200e-06 #based on fwhm of 1 nC, 3 keV beam
    
    sigX, sigY = pulseWidth(ekev)/fwhm2rms, pulseWidth(ekev)/fwhm2rms
    pulseEn = pulseEnergy(q, ekev)
    
    dtheta = pulseDivergence(q, ekev)
    tau = pulseDuration(q)
    gsnBm = moddedGaussian(nx, ny, nz, 
                          ekev, 
                          xMin, xMax, yMin, yMax, 
                          tau = tau, 
                          sigX = sigX, sigY = sigY, 
                          xoff = 0, yoff = 0,
                          xtilt = dtheta, ytilt = dtheta,
                          d2waist = 0, 
                          pulseEn = pulseEn, 
                          pulseRange = 2, 
                          _mx = 0, _my = 0)
    
 
    wfr = Wavefront(gsnBm)
    return wfr

def divergenceCorrection(wfr, f):
    """
    """
    bl = Beamline()
    tl = thinLens(_Fx = f, _Fy = f, _x = 0, _y = 0)
    
    drift = Drift(2*f)
    
    pp = [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
    
    #bl.append(tl, pp)
    bl.append(drift, pp)
    
    
    bl.propagate(wfr)
    
    wfr.params.Mesh.zCoord = 0
    return wfr

if __name__ == '__main__':
    ## TEST FOR USAGE
    ### SET PARAMS
    
    q = 1.0 # nC
    ekev = 0.1 # keV
    
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
    
    wfr = coherentSource(256, 256, 1, ekev, q)
    wfr = divergenceCorrection(wfr, tlFocus(pulseWidth(ekev),pulseDivergence(q,ekev)))
    

