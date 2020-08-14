#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 13:42:08 2020

compile the energy statistics of the propagated pulses

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


from os import listdir, getcwd

import numpy as np


def getSystemEff(indir):
    """ 
    Python wrapper to get the system efficiency files from ind. measurements
    
    :param indir: path to npy files containing system efficiency vals
    :returns v: system efficiency values
    """    
    v = []
    
    for f in listdir(indir):
        v.append(np.load(indir + f)[0])
    
    v = np.asarray(v)
    
    return v
    
    
    
def main(indir):
    
    v = getSystemEff(indir)
    
    print("Average System Efficiency: {}".format(np.mean(v)))
    print("Std. Dev of System Efficiency: {}".format(np.std(v)) )
    print("Percentage Error of System Eff: {}".format(np.std(v)/np.mean(v)*100))


if __name__ == '__main__':
    
    indir = "../../data/NanoKB/systemEff/"
    
    main(indir)
    
    


