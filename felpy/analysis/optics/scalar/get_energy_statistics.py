 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 19:11:57 2020

@author: twguest
"""

import numpy as np

from multiprocessing import cpu_count, Pool
 
class get_energy_statistics():
    """
    get the energy of an wavefield from its intensity map
    
    :param ii: intensity image [2D np array]
    :param dx: extent along horizontal axis (m) [float]
    :param dy: extent along vertical axis (m) [float]
    :param dt: extent in temporal axis (m) [float]
    
    :returns pulse_energy_J: pulse energy in mJoules [float]
    :returns photons_per_pulse: number of photons per pulse [int]
    """
    
    def __init__(self, ii, dx, dy, dt, ekev,
                 mpi = False,
                 VERBOSE = False):
        
        self.ii = ii
        self.dx = dx
        self.dy = dy
        self.dt = dt
        self.ekev = ekev
        self.mpi = mpi
        self.VERBOSE = VERBOSE
        
        self.J2eV = 6.24150934e18
    

    def get_pulse_energy(self, ii):
        
        pulse_energy = self.ii.sum()
        pulse_energy = pulse_energy * self.dx * self.dy * self.dt * 1e6 #### mJ
        
        photons_per_pulse = pulse_energy * self.J2eV // (self.ekev/1000)
       
        if self.VERBOSE:
            print("Pulse Energy: {:e} mJ".format(pulse_energy))
            print("Number of photons per pulse: {:e}".format(photons_per_pulse))
        
        return pulse_energy, photons_per_pulse
    
    def run(self):
       
        if self.mpi:
            processes = cpu_count()//2
            
            p = Pool(processes = processes)
            energy_statistics = p.map(self.get_pulse_energy, [self.ii[:,:,itr] for itr in range(self.ii.shape[-1])])
            energy_statistics = np.asarray(energy_statistics).T
            #energy_statistics.reshape([2, self.ii.shape[-1]])
        else:
            
            if self.ii.ndim == 2:
                self.get_pulse_energy(self.ii)
    
            elif self.ii.ndim == 3:
                energy_statistics = np.zeros([2, self.ii.shape[-1]])
                for itr in range(self.ii.shape[-1]):
                        energy_statistics[:,itr] = self.get_pulse_energy(
                            self.ii[:,:,itr])
        
        return energy_statistics



# =============================================================================
# import unittest
# ### unit test 
# 
# class Testing(unittest.TestCase):
# 
#     def test_one(self):
#         array = np.ones([10,10,10])*np.arange(10)
#         
#         dx = 5e-06 ### arb
#         dy = 5e-06 ### arb
#         dt = 25e-15 ### arb
#         
#         ekev = 5000
#         
#         a = get_energy_statistics(array, dx, dy, dt, ekev, mpi = True, VERBOSE = True).run()
#         
# if __name__ == '__main__':
#     unittest.main()
# =============================================================================
