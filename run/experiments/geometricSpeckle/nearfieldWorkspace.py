#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 09:18:44 2020

@author: twguest
"""

import sys

sys.path.append("../../../")

import numpy as np

from model.tools import constructPulse, getCoord
from model.materials.phaseMask import phaseMask

from scipy.signal import argrelextrema
from matplotlib import pyplot as plt



def featureLocations(mask, dx = 1, dy = 1):
    """
    get the location of 
    """
    
    ### get coordinate system
    idx, idy = getCoord(mask)
    idx *= dx
    idy *= dy
    
    ## feature positions
    ## canonical axis repr. is rev (r,c) not (x,y)
    ## NOTE: We may need to change this condition for np.where
    fx, fy = argrelextrema(mask, np.greater) 
    
    
    
    #fx = fx.astype('float64')
    #fy = fy.astype('float64')
    
    #fx = dx * fx
    #fy = dy * fy
    
    fpos = []    

    for i in range(idx.shape[0]):
        fpos.append((fx[i], fy[i]))
     

    ## DEBUG t = 19
    ## DEBUG print("x: {}".format(fx[t]))
    ## DEBUG print("y: {}".format(fy[t]))
    ## DEBUG print(mask[fy[t], fx[t]]) 
     
    return fpos

def constructMask(nRepeats = 9):
    """
    we need to construct a virtual mask:
        - a periodic array of random uniform elements w/ a known size
        and position in optical space.
    """
    
    if (nRepeats%2) == 0:
        raise(ValueError("Number of Rpeates Should Be Odd"))
        
    mask = np.zeros([3,3])
    mask[1,1] = 1
    
    mask = np.tile(mask, [nRepeats, nRepeats]) 
    mask *= np.random.rand(*mask.shape)
    idx, idy = getCoord(mask)
    
    return mask


    
def run():
    """
    this runs the virtual experiment
    """
    nFeatures = 100
    fWindow = 9
    
    wfr = constructPulse(nx = nFeatures*fWindow, ny = nFeatures*fWindow,
                         nz = 10)
    
    dx, dy = wfr.pixelsize()
    
    constructMask(nRepeats = 9)
    featureLocations(mask, dx = dx, dy = dy)

if __name__ == '__main__':

    mask = constructMask(nRepeats = 9)
    plt.imshow(mask); plt.show()
    fpos = featureLocations(mask)