#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.0"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""

import numpy as np

from scipy.constants import c,h,e

def wrap_phase(phase):
    return (np.pi * phase % (2 * np.pi) - np.pi)

def get_phase_shift(OPD, wav):
    return (2*np.pi*OPD)/wav

def ekev2wav(ekev):
    return (h*c)/(e*ekev*1000)

def ekev2k(ekev):
    return (np.pi*2)/ekev2wav(ekev)

def geometric_focus(fwhm, divergence):
    """
    calculate the focal length of a converging beam of size fwhm
    """
    return fwhm/np.tan(divergence)

def get_magnification(src2sample, sample2det):
    """
    determine the magnification (M > 1 is a demagnification)

    :param src2sample: source to sample distance [m]
    :param sample2det: sample to detector distance [m]
    """

    M = (src2sample+sample2det)/src2sample

    return M

def fresnel_criterion(W, z, wav):
    """
    determine the frensel number of a wavefield propagated over some distance z

    :param W: approx. beam/aperture size [m]
    :param z: propagation distance [m]
    :param wav: source wavelength [m]

    :returns F: Fresnel number
    """

    F = W**2/(z*wav)
    return F

def get_required_distance(W, sigma_det, wav):

    """
    Calculate the propagation distance required to satisfy sampling conditions.

    :param W: approximate feature size [m]
    :param sigma_det: propagated plane (detector) pixel size [m]
    :param wav: source wavelength [m]

    :returns zreq: required distance to satisfy sampling conditions [m]
    """

    if type(sigma_det) == list:
        sigma_det = max(sigma_det)

    zreq = (W*sigma_det)/(wav)
    return zreq


def get_image_properties(W, M, d_get_spatial_resolution, d_npix):
    """
    Return degree to which beam FWHM fills detector

    :param W: approx. beam size [m]
    :param M: magnification factor
    :param d_get_spatial_resolution: detector get_spatial_resolution [m]
    :param d_npix: number of pixels at detector
    """

    if type(d_get_spatial_resolution) == list:
        d_px = d_get_spatial_resolution[0]
        d_py = d_get_spatial_resolution[1]
    else:
        d_px = d_py = d_get_spatial_resolution

    if type(d_npix) == list:
        d_nx = d_npix[0]
        d_ny = d_npix[1]
    else:
        d_nx = d_ny = d_npix


    ## get fwhm and divergence of
    mW = M*W


    ## calculate detector size
    dx = d_nx*d_px
    dy = d_ny*d_py

    ## calculate the percentage of the detector thats filled in each direction
    fill_x = mW/dx
    fill_y = mW/dy

    print("% of Detector Filled By Beam in X-direction: {} %".format(fill_x*100))
    print("% of Detector Filled By Beam in Y-direction: {} %".format(fill_y*100))

    ## calculate the number of pixels the beam occupies in the detector plane
    bx = mW/d_px
    by = mW/d_py

    print("Number of Pixels Beam Occupies in X-direction: {} pixels".format(int(bx)))
    print("Number of Pixels Beam Occupies in Y-direction: {} pixels".format(int(by)))
