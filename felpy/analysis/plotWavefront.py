#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 19:58:13 2020

@author: twguest

A set of high-level functions 
"""
 
import sys
from os import listdir
sys.path.append("../../") # LOCAL PATH
 
import numpy as np
from utils.vis_utils import triple


def pick3(indir):
    """
    from this directory, pick 3 wavefronts at random
    """
   
    l = listdir(indir)
    rints = np.random.randint(low = 0, high = len(l), size = 3)
    
    arr = []
    
    for r in rints:
        arr.append(np.load(indir + l[r]))
    
    arr = np.rollaxis(np.stack(arr), 0, 3)
    
    return arr
    
def plotIntensity(iii):
    """
    plot a set of 3 intensities
    """    
    triple(iii,
           cmap = 'hot',
           clabel = "Intensity (W/mm$^2$)")
    

def plotPhase(ph):
    
    triple(ph,
           cmap = 'hsv',
           clabel = "Phase (rad.)",
           vmin = -np.pi,
           vmax = np.pi, 
           cticks = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi], 
           clabels = ["$-\pi$", "", "0",  "", "$\pi$"])
    
        
def main():
    
    iii = np.random.uniform(-np.pi, np.pi, size = [100,100,3])  
    
    plotIntensity(iii)
    plotPhase(iii)
if __name__ == '__main__':
    
    main()