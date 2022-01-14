# -*- coding: utf-8 -*-

import scipy

import numpy as np

"""
this file contains the analytical data describing the SA1 beamline
"""

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

    sig = np.log((7.4e03/ekev))*6
    return sig/1e6


def analytical_pulse_divergence(ekev, limit = 'upper'):
    
    """
    Estimate of analytical_pulse_divergence (half-angle) from electron bunch charge and radiation energy

    :param q: electron bunch charge [nC]
    :param ekev: radiation energy [keV]

    :return dtheta: pulse divergence [rad]
    """
    
    if limit == 'lower':
        return (8.76)/(1e06*ekev**0.85)
    elif limit == 'upper':
        return (14.1)/(1e06*ekev**0.75)
