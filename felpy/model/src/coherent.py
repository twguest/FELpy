#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""


import sys


import scipy

import numpy as np

from wpg import srwlib
from wpg.srw import srwlpy

from copy import deepcopy

from scipy.constants import c
from matplotlib import pyplot as plt


import wpg.srwlpy
from felpy.utils.opt_utils import ekev2wav
from felpy.model.core.wavefront import Wavefront
from wpg.wpg_uti_wf import look_at_q_space, calculate_fwhm, calc_pulse_energy
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from wpg.wpg_uti_wf import plot_intensity_qmap as plotKspace


from wpg.generators import build_gauss_wavefront_xy as build_gaussian
from wpg.generators import build_gauss_wavefront as build_gaussian_3D

from felpy.model.core.beamline import Beamline

from wpg.srwlib import SRWLOptL as thinLens

from wpg.optical_elements import Drift
from wpg.srwlib import SRWLOptD

fwhm2rms = np.sqrt(8*np.log(2)) ### FWHM = sqrt(8ln(2))*sigma
h = scipy.constants.physical_constants['Planck constant in eV s'][0]

def get_pulse_energy(wfr):
    wfr.set_electric_field_representation('t')
    energy = calc_pulse_energy(wfr)[0]
    wfr.set_electric_field_representation('f')
    return energy
    
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

    sig = np.log((7.4e03/ekev))#*6
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
    bl.append(tl, [0,0,0,0,0,1,1,1,1,0,0,0])
    bl.append(Drift(2*f), [0,0,1,0,0,1,1,1,1,0,0,0])
    bl.propagate(wfr)

    wfr.params.Mesh.zCoord = 0

def modify_beam_energy(wfr, desired):
    
    
    actual = get_pulse_energy(wfr)
    diff = actual/desired
    
    wfr.data.arrEhor *= diff


def construct_gaussian(nx, ny, ekev, extent, sigX, sigY, divergence, xoff = 0, yoff = 0, tau = 1, mx = 0, my = 0, tiltX = 0, tiltY = 0):
    gsnBm = build_gaussian(nx, ny, ekev, xMin = extent[0], xMax = extent[1], yMin = extent[2], yMax = extent[3],
                           sigX = sigX, sigY = sigY, d2waist = 1, xoff = xoff, yoff = yoff, pulseTau = tau, _mx = mx, _my = my, tiltX= tiltX, tiltY = tiltY)

    wfr = Wavefront(gsnBm)
    wfr.params.wavelength = ekev2wav(ekev)
    
    modify_beam_divergence(wfr,sigX, divergence) 
    return wfr



def construct_SA1_wavefront(nx, ny, ekev, q, xoff = 0, yoff = 0, mx = 0, my = 0,
                            tiltX = 0, tiltY = 0):
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


    gsnBm = build_gaussian(nx, ny, ekev, xMin, xMax, yMin, yMax, sigX, sigY, 1, xoff = xoff, yoff = yoff, pulseTau = tau, _mx = mx, _my = my, tiltX= tiltX, tiltY = tiltY)

    wfr = Wavefront(gsnBm)
    wfr.params.wavelength = wavelength
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')

    modify_beam_divergence(wfr,analytical_pulse_width(ekev), analytical_pulse_divergence(q,ekev)) ## note 300 is hardcoded fix, should not stay past tue 13/01/21



    return wfr


def construct_SA1_pulse(nx, ny, nz, ekev, q, modify_beam_divergenceergence = True, mx = 0, my = 0):
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

    gsnBm = build_gaussian_3D(nx = nx,
                              ny = ny,
                              nz = nz,
                              ekev = ekev,
                              xMin = xMin, xMax = xMax,
                              yMin = yMin, yMax = yMax,
                              sigX = sigX, sigY = sigY,
                              d2waist = 1,
                              tau = tau,
                              pulseRange = 1,
                              pulseEn = 2*pulseEn,
                              _mx = mx,
                              _my = my)

    wfr = Wavefront(gsnBm)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')

    wfr.params.wavelength = wavelength
    modify_beam_divergence(wfr,analytical_pulse_width(ekev),analytical_pulse_divergence(q,ekev))
    
    modify_beam_energy(wfr, desired = 2*pulseEn)
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
            divergences.append(wfr.get_divergence()[0])
        ax.scatter(energyrange, np.array(divergences)*1e6, marker = 'x')

        analytical_data = [analytical_pulse_divergence(q, energy) for energy in energyrange]
        ax.plot(energyrange, np.array(analytical_data)*1e6, '--')

    ax.legend(["0.01 nC", "0.1 nC", "0.3 nC", "0.5 nC", "0.8 nC", "1.0 nC"])

    if outdir == None:
        pass
    else:
        fig.savefig(outdir + "SourceDivergence_Energy.png")


if __name__ == '__main__':

    pass
