#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 04:00:27 2020

@author: guestt
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
from os import listdir 
from ViewPy.scatters import scatter3

def compileData(indir):
    
    data = []
    
    for f in listdir(indir):
        
        raw = np.load(indir + f)
        
        data.append(raw)
    
    
    return data
    
    
if __name__ == '__main__':
    
    
    indir =  "../../data/NanoKB/measBeamSize/"
    outfile = "../../data/proc/measBeamSize"
    
    r = compileData(indir)
    
    
    rad = [float(a[1])*1e6 for a in r] # in nm
    xpos = [float(a[3])*1e6 for a in r]
    ypos = [float(a[4])*1e6 for a in r]
    
    
    print("Average Beam Radius: {}".format(np.mean(rad)))
    print("Std. Dev Beam Radius: {}".format(np.std(rad)))
    

    #scatter3(xpos,ypos,rad)