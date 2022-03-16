# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.1"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""

import numpy as np

from matplotlib import pyplot as plt
from scipy.constants import c,h,e

def ekev2wav(ekev):
    """
    convert energy to wavelength
    
    :params ekev: radiation energy in keV
    
    :returns: radiation wavelength in m
    """
    return (h*c)/(e*ekev*1000)

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

def nfs_sampling(dz, feature_size, dx, wav):
    """
    
    :param dz: distance from the speckle mask to the detector
    :param feature_size: size of the phase-object for a nfs experiment
    :param dx: detector pixel-size
    :param wav: radiation wavelength
    
    :returns S: sampling criterion (dz/zreq)
    :return F: Fresnel Criterion
    """
    
    S = dz/get_required_distance(feature_size, dx, wav) ### required prop. dist for sampling 
    F = fresnel_criterion(feature_size, dz, wav) ### fresnel criterion
    
    return S,F

if __name__ == '__main__':
    
    ### usage
    z = 1.5 ### focus to sample distance
    dz = 1.2 ### sample to detector distance
    feature_size = 4e-06 ### mask feature size
    dx = 12.5e-06 ### detector pixel size
    ekev = 9 ### energy in keV
    M = (z+dz)/z ### effective demagnification
    z_eff = (z*dz)/(z+dz) ### effective propagation distance via fresnel scaling
    print("Z effective: {}".format(z_eff))
    sigma = 50e-06 ### estimated beam size @ sample
    
    print("Estimated Beam Footprint: {}".format(sigma*M))

    S = get_required_distance(feature_size, dx, ekev2wav(ekev))-dz
    print("Sampling Criterion (req. > 1): {}".format(S))

    F = fresnel_criterion(feature_size, z_eff, ekev2wav(ekev))
    print("Fresnel Criterion (req. > 1): {}".format(F))
    
    #print("Feature Size @ Detector: {} (should be > dx) m".format(M*feature_size))
    print("Approx. Pixels per Detected Feature: {}".format(M*feature_size/dx))
    
    print("Angular Sensitivity: {}".format(np.arctan(dx/dz)))
    
    print()
    print("NANO")
    ### usage
    z = .1 ### focus to sample distances
    dz = 10 ### sample to detector distance
    feature_size = 4e-06 ### mask feature size
    dx = 25e-06 ### detector pixel size
    ekev = 9 ### energy in keV
    M = (z+dz)/z ### effective demagnification
    z_eff = (z*dz)/(z+dz) ### effective propagation distance via fresnel scaling
    print("Z effective: {}".format(z_eff))
    sigma = 100e-06 ### estimated beam size @ sample
    
    print("Estimated Beam Footprint: {}".format(sigma*M))

    S = get_required_distance(feature_size, dx, ekev2wav(ekev))-dz
    print("Sampling Criterion (req. > 1): {}".format(S))

    F = fresnel_criterion(feature_size, z_eff, ekev2wav(ekev))
    print("Fresnel Criterion (req. > 1): {}".format(F))
    
    #print("Feature Size @ Detector: {} (should be > dx) m".format(M*feature_size))
    print("Approx. Pixels per Detected Feature: {}".format(M*feature_size/dx))
    
    print("Angular Sensitivity: {}".format(np.arctan(dx/dz)))
    