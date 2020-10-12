 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 19:11:57 2020

@author: twguest
"""

import numpy as np


def get_pulse_energy(ii, dx, dy, dt, photonEnergy, VERBOSE = False):
    """
    get the energy of an image from its intensity map
    
    :param ii: intensity image [2D np array]
    :param dx: extent along horizontal axis (m) [float]
    :param dy: extent along vertical axis (m) [float]
    :param dt: extent in temporal axis (m) [float]
    
    :returns pulse_energy_J: pulse energy in Joules [float]
    :returns photons_per_pulse: number of photons per pulse
    """

    J2eV = 6.24150934e18
 
    pulse_energy = ii.sum()
    pulse_energy_J = pulse_energy * dx * dy * 1e6 * dt
    
    photons_per_pulse = pulse_energy_J * J2eV / photonEnergy
       
    if VERBOSE:
        print("Pulse Energy: {:e} J".format(pulse_energy_J))
        print("Number of photons per pulse: {:e}".format(photons_per_pulse))
    
    return pulse_energy_J, photons_per_pulse



if __name__ == '__main__':
    
    arr = np.random.rand(100,100)*1e13
    
    dx = 5e-06
    dy = 5e-06
    dt = 25e-15
    
    ekev = 5000
    
    get_pulse_energy(arr, dx, dy, dt, ekev, VERBOSE = True)