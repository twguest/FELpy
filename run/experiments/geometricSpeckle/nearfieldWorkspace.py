#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 09:18:44 2020

@author: twguest
"""

import sys

sys.path.append("../../../")

import numpy as np

from model.tools import constructPulse
from model.materials.phaseMask import phaseMask

from matplotlib import pyplot as plt


def featureLocations(mask, dx, dy):
    ## get index of positions
    ## translate index by wfr.pixelsize
    pass

def constructMask():
    """
    we need to construct a virtual mask:
        - a periodic array of random uniform elements w/ a known size
        and position in optical space.
    """
    
    mask = np.zeros([3,3])
    mask[1,1] = 1
    mask = np.tile(mask, [10, 10])
    
    
    
    return mask

def run():
    """
    this runs the virtual experiment
    """
    pass

if __name__ == '__main__':
    mask = constructMask()
    plt.imshow(mask)
    print(mask[mask.shape[0]//2], mask[mask.shape[1]//2])