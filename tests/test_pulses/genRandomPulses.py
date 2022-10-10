#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 15:09:54 2020

@author: twguest
"""

###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################

import numpy as np

from felpy.model.wavefront import Wavefront
from wpg.generators import build_gauss_wavefront


def generatePulse(n = 1):

    
    wfr = Wavefront(build_gauss_wavefront(nx = 512, ny = 512, nz = 5,
                                          ekev = 4.96,
                                          xMin = -400e-06, xMax = 400e-06,
                                          yMin = -400e-06, yMax = 400e-06,
                                          tau = 1e-05,
                                          sigX = np.random.random()*125e-06, sigY = np.random.random()*125e-06,
                                          d2waist = 1))
                    
    wfr.store_hdf5("../../data/tests/pulseTests/gaussianSource/gsn_{}.h5".format(n))

def testGeneratedPulses(indir, fname, n):
    
    wfr = Wavefront()
    wfr.load_hdf5(indir + fname + "_{}.h5".format(n))

        
def main():
    
    N = 25
    
    for n in range(N):
        
        if n % 5 == 0:
            print("Generating Test Pulse: {} of {}".format(n,N))
        
        generatePulse(n)
    
    for n in range(N):
        
        testGeneratedPulses("../../data/tests/pulseTests/gaussianSource/", "gsn", n)
             
        if n % 5 == 0:
            print("Wavefront {} Loaded Succesfully".format(n))


if __name__ == "__main__":
    main()