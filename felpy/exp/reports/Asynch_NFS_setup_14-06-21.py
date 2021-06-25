# -*- coding: utf-8 -*-

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



def plot_nfs_sampling(zrange, feature_size, dx, wav):

    fig, ax1 = plt.subplots()
 
    S, F = nfs_sampling(zrange, feature_size, dx, wav)
    ax1 =   simple_line_plot(zrange, S,
                             label = "Sampling Factor",
                             return_axes = True, parse_axes = ax1,
                             xlabel = "Distance from Mask (m)", ylabel = "Sampling Factor")
    #plt.legend()
    ax2 = ax1.twinx()
    
                
    ax2 =   simple_line_plot(zrange, F,
                             label = "Fresnel Criterion".format(wav*1e10),
                             return_axes = True, parse_axes = ax2,
                             xlabel = "Distance from Mask (m)", ylabel = "Fresnel Criterion",
                             color = 'red')
    #plt.legend()
    ax1.set_ylim([1,10])
    ax2.set_ylim([1,10])
    plt.show()
    

def experiment_report_in_email():
    
    for dx in [15e-06]:
        
        for ekev in [25]:
            dx = 6.5e-06
            z = 4.5
            feature_sizes = np.linspace(1e-06, 25e-06, 50)
            ekevs = np.linspace(6e-06, 12e-06)
            
            
            S,F = nfs_sampling(z, feature_sizes, dx, ekev2wav(ekev))
            
            fig, ax1 = plt.subplots()
            
            ax1.plot(feature_sizes*1e6, S, 'red', label = "Sampling Fraction")
            ax1.plot(feature_sizes*1e6, np.ones(feature_sizes.shape[0]), 'red')
            plt.legend()
            
            ax2 = ax1.twinx()
            ax2.plot(feature_sizes*1e6, F, 'blue', label = "Fresnel Criterion")
            ax2.plot(feature_sizes*1e6,  np.ones(feature_sizes.shape[0]), color = 'b')
            
            ax2.spines["left"].set_color("red")
            ax2.spines["right"].set_color("blue")
            
            ax1.set_xlabel("Feature Size ($\mu m$)")
            ax1.set_ylabel("Sampling Fraction")
            ax2.set_ylabel("Fresnel Criterion")
            
            ax2.set_title("E: {:} keV \nPixel Size {:} m".format(ekev, dx))
            plt.legend()
            #plt.savefig("/opt/admin/project_meetings/03-05-21/nfs_exp_report_{:.2f}kev_{:.2f}um.png".format(ekev, dx*1e6))

if __name__ == '__main__':
    
    ### usage
    z = 2 ### focus to sample distance
    dz = 4.5 ### sample to detector distance
    feature_size = 15e-06 ### mask feature size
    dx = 6.5e-06 ### detector pixel size
    ekev = 25.0 ### energy in keV
    
    S, F = nfs_sampling(dz, feature_size, dx, ekev2wav(ekev))
    
    print("Sampling Criterion (req. > 1): {}".format(S))
    if S>1: print("good")
    else: print("bad")
    print("Fresnel Criterion (req. > 1): {}".format(F))
    if F>1: print("good")
    else: print("bad")
    M = (z+dz)/z
    print("Feature Size @ Detector: {} (should be > dx) m".format(M*feature_size))
    if M*feature_size>dx: print("good")
    else: print("bad")
    print("Angular Sensitivity: {}".format(angular_sensitivity(6.5e-06, dz)))