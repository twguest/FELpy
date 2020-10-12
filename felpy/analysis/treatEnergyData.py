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
    
    
    indir =  "../../data/NanoKB/pulseEnergy/"
    outfile = "../../data/proc/pulseEnergy"
    

    data = compileData(indir)
    
    en = [d[0] for d in data]
    
    print("Average Photon Energy: {} mJ".format(np.mean(en)*1e3))
    print("Std. Dev in Photon Energy: {} mJ".format(np.std(en)*1e3))
    ph = [d[1] for d in data]
    
    print("Average nPhontons: {:.4e}".format(np.mean(ph)))
    print("Std. Dev in nPhotons: {:.4e}".format(np.std(ph)))