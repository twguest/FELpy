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

def angular_sensitivity(dx, dz):
    return np.arctan(dx/dz)


def theta(x,r, delta = 4e-07):
    return (2*r*delta)/np.sqrt(r**2-(x-r)**2)

def shift(theta, dx, z):
    return np.tan(theta)*z/dx
    
 
if __name__ == '__main__':
    
    ### usage
    
    ### Define the magnifcation from the source aperture at each of the imaging
    ### positions
    
    M01 = 2.50 
    M12 = 3.00
    M02 = 7.50
    
    
    ### Define z1: the source to sample distance, and dz: the sample to det. distance
    z1 = 2.0 ### src to sample distance
    dz = 3.0 ### sample to detector distance
    
    r = 75e-03
    
    
    dx = 6.5e-06 ### detector pixel size
    ekev = 25.0 ### energy in keV
    
    dw = 15e-03 ### detector width 
    fov = np.linspace(dw/M02,r, 2560)  ### detector field of view of source/cyl
    x = np.linspace(0,dw,2560)

    
    ### plotting stuff here
    from felpy.utils.vis_utils import Grids
    plots = Grids(global_aspect = 2, scale = 2)
    plots.create_grid(n = 1, m = 2, sharex = False, sharey = False)
    [ax1, ax2] = plots.axes
    
    ax1.set_xlabel("x (mm)")
    ax1.set_ylabel("Phase Gradient ($\mu rad$)")
    ax2.set_xlabel("x (mm)")
    ax2.set_ylabel("Feature Displacement (pixels.)")

    for r in [.75, 1.00, 1.50]:   

        
        t = theta(fov, r)
        t[np.where(np.isnan(t))] = 0

        ax1.plot(x*1e3, t*1e6, label = "{} mm".format(r*1e2))
        ax2.plot(x*1e3, shift(t, dx, dz+z1))

    for feature_size in [6e-06, 9e-06, 15e-06, 30e-06, 60e-06]:
        print("")
        print("**** Feature Size: {}".format(feature_size))  
        
        S, F = nfs_sampling(dz, feature_size, dx, ekev2wav(ekev))
        print("Sampling Criterion (req. > 1): {}".format(S))
        
        
        print("Fresnel Criterion (req. > 1): {}".format(F))
        
        
        print("Feature Size @ Detector: {} (should be > dx) m".format(M12*feature_size))
        
        #print("Angular Sensitivity: {}".format(angular_sensitivity(dx, dz)))
    
    ax1.legend()
  