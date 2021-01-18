#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 10:29:23 2020

@author: twguest
"""

import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

import scipy

import numpy as np

from wpg import srwlib
from wpg import srwlpy

from copy import deepcopy

from scipy.constants import c
from matplotlib import pyplot as plt


import wpg.srwlpy

from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import look_at_q_space, calculate_fwhm
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from wpg.wpg_uti_wf import plot_intensity_qmap as plotKspace

from wpg.misc import calcDivergence, toKspace

from wpg.generators import build_gauss_wavefront_xy as build_gaussian
from wpg.generators import build_gauss_wavefront as build_gaussian_3D

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
    f = sig/(np.tan(dtheta))
    return f


def modDiv(wfr, f):
    """
    Construct a thin-lens with a focal length defined by a desired divergence,
    see: tlFocus, and then propagate to 2f
    
    :param wfr: WPG wfr structure
    :param f: focal length of thin-lens (x = y) [m]
    """    
    fwhm_i = calculate_fwhm(wfr)
    
    tl = thinLens(f,f)
    bl = Beamline()
    bl.append(tl, [0,0,1,0,0,1,1,1,1,0,0,0])
    bl.append(Drift(2*f), [0,0,1,0,0,1,1,1,1,0,0,0])
    bl.propagate(wfr)
    
    fwhm_f = calculate_fwhm(wfr)
    
    wfr.params.Mesh.zCoord = 0
    

def coherentSource(nx, ny, ekev, q, xoff = 0, yoff = 0, modDivergence = True):
    """
    Construct a fully-coherent Gaussian source with properties related to the
    energy of the radiation and beam charge.
    
    Important points: pulseTau (coherence length) is set to tau (pulse length)
    in the fully coherent case
    
    :param nx: number of horizontal pixels [int]
    :param ny: number of vertical pixels [int]
    :param nz: number of slices [int]
    :param ekev: radiation energy [keV]
    :param q: electron beam bunch charge [nC]
    :param xoff: horizontal offset of beam maximum [m]
    :param yoff: vertical offset of beam maximum [m]
    :param modDivergence: Choose to set non-diffraction limited div [bool]
    """
    wavelength = (h*c)/(ekev*1e3)
        
    xMin, xMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    yMin, yMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    
    sigX, sigY = pulseWidth(ekev)/fwhm2rms, pulseWidth(ekev)/fwhm2rms
    pulseEn = pulseEnergy(q, ekev)
    
    dtheta = pulseDivergence(q, ekev)
    tau = pulseDuration(q)
    
    f = tlFocus(pulseWidth(ekev),pulseDivergence(q,ekev))
        
    gsnBm = build_gaussian(nx, ny, ekev, xMin, xMax, yMin, yMax, sigX, sigY, 1, xoff = xoff, yoff = yoff, pulseTau = tau)
    
    wfr = Wavefront(gsnBm)
    wfr.params.wavelength = wavelength
    modDiv(wfr,f)
    

    return wfr


def tstCoherentSrcFWHM(outdir = None):
    """
    Test that the FWHM of the source matches the analytical prediction
    
    :param outdir: output directory [str]
    """
    print("Testing Source Size vs. Energy and Resolution")
    
    fig = plt.figure(figsize = [12,8])
    ax = fig.add_subplot(111)
    ax.set_title("Source Size Radiation Dependence")
    ax.set_ylabel("FWHM [$\mu$m]")
    ax.set_xlabel("Radiation Energy [keV]")
    ax.set_xlim([2.5,17.5])
    ax.set_ylim([30,55])
    
    dfig = plt.figure()
    dax = dfig.add_subplot(111)
    dax.set_title("Average Error of Modelled and Analytical Source Size")
    dax.set_ylabel("% Error$_{FWHM}$ ($\Delta FWHM/FWHM$)")
    dax.set_xlabel("Pixels")
    dax.set_xlim([0, 5000])
    dax.set_ylim([0, 10])
    
    esteps = 10
    rsteps = 5
    
    energyrange = np.linspace(3,16,esteps)
    resolutionrange = [128, 256, 512, 1024, 2048, 4096]#np.linspace(128, 4096, rsteps)
    analytical_data = np.array([pulseWidth(a)*1e6 for a in energyrange])
    
    mean_delta = []

    for val in resolutionrange:
        val = int(val)
        nx, ny = val, val
        print("Calculating FWHM for {}x{} Array".format(nx, ny))
        fwhms = []
        for energy in energyrange:
            wfr = coherentSource(nx, ny, energy, 0.1)
            fwhms.append(calculate_fwhm(wfr)['fwhm_x']*1e6)
        
        mean_delta.append(np.mean(abs((analytical_data-np.array(fwhms))/analytical_data))*100)
        
        ax.plot(energyrange, fwhms)
    print(len(mean_delta))
    dax.scatter(resolutionrange, mean_delta)
     
    
    leg1 = ax.legend(["{}x{} pixels".format(int(val), int(val)) for val in resolutionrange])
    ax.add_artist(leg1)
    analytical, = ax.plot(energyrange, analytical_data, '--b',)

    leg2 = ax.legend([analytical], ["Analytical Model"], loc = 'lower left')
    ax.add_artist(leg2)
    plt.show(ax)
    
    if outdir is not None:
        if type(outdir) != str:
            print("outdir should be a string, saving failed")
            pass
        else:
            fig.savefig(outdir + "SourceSize_EnergyDep.png")
            dfig.savefig(outdir + "SourceError_pixels.png")
            

    
    wfr = Wavefront(build_gaussian_3D(nx, ny, nz, ekev, -400e-06, 400e-06, -400e-06, 400e-06, tau, 5e-06, 5e-06, d2waist))
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    #look_at_q_space(wfr)
    return wfr

def construct_spb_pulse(nx, ny, nz, ekev, q, modDivergence = True):
    """
    Construct a fully-coherent Gaussian source with properties related to the
    energy of the radiation and beam charge.
    
    Important points: pulseTau (coherence length) is set to tau (pulse length)
    in the fully coherent case
    
    :param nx: number of horizontal pixels [int]
    :param ny: number of vertical pixels [int]
    :param nz: number of slices [int]
    :param ekev: radiation energy [keV]
    :param q: electron beam bunch charge [nC]
    :param xoff: horizontal offset of beam maximum [m]
    :param yoff: vertical offset of beam maximum [m]
    :param modDivergence: Choose to set non-diffraction limited div [bool]
    """
    wavelength = (h*c)/(ekev*1e3)
        
    xMin, xMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    yMin, yMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    
    sigX, sigY = pulseWidth(ekev)/fwhm2rms, pulseWidth(ekev)/fwhm2rms
    pulseEn = pulseEnergy(q, ekev)
    
    dtheta = pulseDivergence(q, ekev)
    tau = pulseDuration(q)
    
    f = tlFocus(pulseWidth(ekev),pulseDivergence(q,ekev))
        
    gsnBm = build_gaussian_3D(nx = nx,
                              ny = ny,
                              nz = nz,
                              ekev = ekev,
                              xMin = xMin, xMax = xMax,
                              yMin = yMin, yMax = yMax,
                              sigX = sigX, sigY = sigY,
                              d2waist = 1,
                              tau = tau,
                              pulseRange = 1)
    
    wfr = Wavefront(gsnBm)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
    wfr.params.wavelength = wavelength
    modDiv(wfr,f)
    
    return wfr


def tstCoherentSrcDiv(outdir = None):
    """
    Test that the FWHM of the source matches the analytical prediction
    
    :param outdir: output directory [str]
    """
    
    
    fig = plt.figure(figsize = [12,8])
    ax = fig.add_subplot(111)
    ax.set_title("Source Divergence FWHM", fontsize = 22)
    ax.set_ylabel(r"d$\theta$ [$\mu$rad]", fontsize = 16)
    ax.set_xlabel("Radiation Energy [keV]", fontsize = 16)
    ax.set_xlim([2.5,17.5])
    
    
    
    nx, ny = 512, 512
    
    esteps = 5
    energyrange = np.linspace(3,16, esteps)
    
    for q in [0.01, 0.1, 0.3, 0.5, 0.8, 1.0]:
        
        divergences = []
        for energy in energyrange:
            wfr = coherentSource(nx, ny, energy, q)
            divergences.append(calcDivergence(wfr)[0])
        ax.scatter(energyrange, np.array(divergences)*1e6, marker = 'x')
        
        analytical_data = [pulseDivergence(q, energy) for energy in energyrange]
        ax.plot(energyrange, np.array(analytical_data)*1e6, '--')
        
    ax.legend(["0.01 nC", "0.1 nC", "0.3 nC", "0.5 nC", "0.8 nC", "1.0 nC"])
    
    if outdir == None:
        pass
    else:
        fig.savefig(outdir + "SourceDivergence_Energy.png")

if __name__ == '__main__':
    
    ## TEST FOR USAGE
    print("Testing Coherent Gaussian Source Module")
    
    ### SET PARAMS
    print("Testing Parameters @ q = 0.1 nC and ekev = 9.2 keV\n")    
    q = 0.1 # nC
    ekev = 9.2 # keV
    
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
    print("\n")
    print("Thin-Lens Focus: {} m".format(tlFocus(pulseWidth(ekev),pulseDivergence(q,ekev))))
    
    ### Generate Coherent Source
    wfr = coherentSource(1024, 1024, 9.2, 0.1)
    
    ### Test Source FWHM
    #tstCoherentSrcFWHM("/opt/spb_model/tests/")
    
    ### Test Source Divergence
    tstCoherentSrcDiv("/opt/spb_model/tests/")     
    
    ### Plot Source 
    #plotIntensity(wfr, "/opt/spb_model/tests/9200eV_100pC_intensity.png")
    #plotKspace(wfr, "/opt/spb_model/tests/9200eV_100pC_kspace.png")
 
 
