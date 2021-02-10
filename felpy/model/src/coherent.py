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

def analytical_pulse_energy(q, ekev):
    """
    Estimate of analytical_pulse_energy from electron bunch charge and radiation energy    
    
    :param q: electron bunch charge [nC]
    :param ekev: radiation energy [keV] 
        
    :return P: pulse energy [J]
    """
    
    P = 19*q/ekev
    return P/1e3

def analytical_pulse_duration(q):
    """
    Estimate analytical_pulse_duration from electron bunch charge 
    
    :param q: electron bunch charge [nC]
        
    :return t: Duration of pulse [s]
    """
    
    t = (q*1e3)/9.8
    return t*1e-15


def analytical_pulse_width(ekev):
    """
    Estimate analytical_pulse_width (FWHM) from radiation energy (assumes symmetrical beam)
    
    :param ekev: radiation energy [keV]
        
    :return sig: Radiation pulse width (FWHM) [m]
    """
    
    sig = 6*np.log((7.4e03/ekev))
    return sig/1e6

def analytical_pulse_divergence(q, ekev):
    """
    Estimate of analytical_pulse_divergence (half-angle) from electron bunch charge and radiation energy    
    
    :param q: electron bunch charge [nC]
    :param ekev: radiation energy [keV] 
        
    :return dtheta: pulse divergence [rad]
    """
    
    dtheta = (17.2*6.4*np.sqrt(q))/ekev**(0.85)
    return dtheta/1e6


def modify_beam_divergence(wfr, sig, dtheta):
    """
    Construct a thin-lens with a focal length defined by a desired divergence,
    see: thin_lens_mod, and then propagate to 2f
    
    :param wfr: WPG wfr structure
    :params sig: beam fwhm [m]
    :params dtheta: pulse divergence [rad]
    """    
    f = sig/(np.tan(dtheta))
    
    tl = thinLens(f,f)
    bl = Beamline()
    bl.append(tl, [0,0,1,0,0,1,1,1,1,0,0,0])
    bl.append(Drift(2*f), [0,0,1,0,0,1,1,1,1,0,0,0])
    bl.propagate(wfr)
        
    wfr.params.Mesh.zCoord = 0
    

def construct_SA1_wavefront(nx, ny, ekev, q, xoff = 0, yoff = 0, modify_beam_divergenceergence = True):
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
    :param modify_beam_divergenceergence: Choose to set non-diffraction limited div [bool]
    """
    wavelength = (h*c)/(ekev*1e3)
        
    xMin, xMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    yMin, yMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    
    sigX, sigY = analytical_pulse_width(ekev)/fwhm2rms, analytical_pulse_width(ekev)/fwhm2rms
    pulseEn = analytical_pulse_energy(q, ekev)
    
    dtheta = analytical_pulse_divergence(q, ekev)
    tau = analytical_pulse_duration(q)
    
<<<<<<< Updated upstream
    
=======
    f = tlFocus(pulseWidth(ekev),pulseDivergence(q,ekev))
        
>>>>>>> Stashed changes
    gsnBm = build_gaussian(nx, ny, ekev, xMin, xMax, yMin, yMax, sigX, sigY, 1, xoff = xoff, yoff = yoff, pulseTau = tau)
    
    wfr = Wavefront(gsnBm)
    wfr.params.wavelength = wavelength
    modify_beam_divergence(wfr,analytical_pulse_width(ekev),analytical_pulse_divergence(q,ekev))
    

    return wfr


def construct_SA1_pulse(nx, ny, nz, ekev, q, modify_beam_divergenceergence = True):
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
    :param modify_beam_divergenceergence: Choose to set non-diffraction limited div [bool]
    """
    wavelength = (h*c)/(ekev*1e3)
        
<<<<<<< Updated upstream
    xMin, xMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    yMin, yMax = -400e-06, 400e-06 #based on fwhm of 1 nC, 3 keV beam
    
    sigX, sigY = analytical_pulse_width(ekev)/fwhm2rms, analytical_pulse_width(ekev)/fwhm2rms
    pulseEn = analytical_pulse_energy(q, ekev)
    
    dtheta = analytical_pulse_divergence(q, ekev)
    tau = analytical_pulse_duration(q)

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
=======
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
            

def construct_pulse(nx = 512, ny = 512, nz = 5, ekev = 5.0, tau = 1e-06, d2waist = 10):
    
    wfr = Wavefront(build_gaussian_3D(nx, ny, nz, ekev, -400e-06, 400e-06, -400e-06, 400e-06, tau, 5e-06, 5e-06, d2waist))
>>>>>>> Stashed changes
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
    wfr.params.wavelength = wavelength
    modify_beam_divergence(wfr,analytical_pulse_width(ekev),analytical_pulse_divergence(q,ekev))
    
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


def test_coherent_pulse_divergence(outdir = None):
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
            wfr = construct_SA1_wavefront(nx, ny, energy, q)
            divergences.append(calcDivergence(wfr)[0])
        ax.scatter(energyrange, np.array(divergences)*1e6, marker = 'x')
        
        analytical_data = [analytical_pulse_divergence(q, energy) for energy in energyrange]
        ax.plot(energyrange, np.array(analytical_data)*1e6, '--')
        
    ax.legend(["0.01 nC", "0.1 nC", "0.3 nC", "0.5 nC", "0.8 nC", "1.0 nC"])
    
    if outdir == None:
        pass
    else:
        fig.savefig(outdir + "SourceDivergence_Energy.png")


if __name__ == '__main__':
    
    from felpy.utils.os_utils import mkdir_p
    
    sdir = '../../tests/coherent_source/'
    mkdir_p(sdir)
    ## TEST FOR USAGE
    print("Testing Coherent Gaussian Source Module")
    
    ### SET PARAMS
      
    q = 0.25 # nC
    ekev = 5.0 # keV
    print("Testing Parameters @ q = {} nC and ekev = {} keV\n".format(q, ekev))  
    
    ### Print Parameters as Sanity Check
    print("Electron Beam Charge: {} nC".format(q))    
    print("Electron Beam Charge: {} pC".format(q*1e3))   
    
    print("Radiation Energy: {} keV".format(ekev))   
    
    print("\n")
    ### Estimate Energy per Pulse
    print("Pulse Energy: {} Joules".format(analytical_pulse_energy(q, ekev)))

    ### Estimate Duration of pulse
    print("Pulse Duration: {} seconds".format(analytical_pulse_duration(q)))
        
    ### Estimate FWHM of pulse (assumes symmetry)
    print("Pulse Width: {} m".format(analytical_pulse_width(ekev)))
    
    ### Estimate FWHM of pulse (assumes symmetry)
    print("Pulse Divergence: {} rad".format(analytical_pulse_divergence(q,ekev)))
    
    ### Estimate thin lens focus for divergence correction
    print("\n")
    print("Thin-Lens Focus: {} m".format(thin_lens_mod(analytical_pulse_width(ekev),analytical_pulse_divergence(q,ekev))))
    
    ### Generate Coherent Source
    wfr = construct_SA1_pulse(nx = 1024, ny = 1024, nz = 5, ekev = ekev, q = q)
    
    ### Test Source FWHM
    test_coherent_pulse_fwhm(sdir)
    
    ### Test Source Divergence
    test_coherent_pulse_divergence(sdir)     
    
    ### Plot Source 
    #plotIntensity(wfr, "/opt/spb_model/tests/9200eV_100pC_intensity.png")
    #plotKspace(wfr, "/opt/spb_model/tests/9200eV_100pC_kspace.png")
 
 
