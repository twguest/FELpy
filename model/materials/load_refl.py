#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 14:32:58 2020

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

import time

import numpy as np
import pandas as pd

def treat_refl(material = 'B4C', indir = None):
    """
    Used for converting excel data obtained via henke to dat file
    
    :param material: mirror material which data is defined for
    """
    
    if indir is None:
        excel = pd.read_excel("../../data/hom_refl/{}_refl.xlsx".format(material),header = None)
    else:
        excel = pd.read_excel(indir + "{}_refl.xlsx".format(material),header = None)
    df = pd.DataFrame(excel)
    
    refl = np.asarray(df.values)
    
    if indir is None: 
        np.save("../../data/hom_refl/{}_refl.npy".format(material), refl)
    else:
        np.save(indir + "{}_refl.npy".format(material), refl)
    
def load_refl(material = 'B4C', indir = None):
    """
    load data from .npy file generate by treat_refl method
    
    :param material: mirror material which data is defined for
    """
    
    if indir is None:
        refl = np.load("../../data/hom_refl/{}_refl.npy".format(material))
    else:
        refl = np.load(indir + "/{}_refl.npy".format(material))
    return refl

def get_refl(refl, ekev, ang = 'max', limits = [0,2*np.pi]):
    
    e_idx = (np.abs(refl[0,1:] -ekev*1e3)).argmin() + 1
    refl_idx = refl[:, e_idx][1: ]
    
    angs = refl[1:,0]
    
    if ang == 'max':
        
        llim = np.min(np.where(angs >= limits[0]))
        ulim = np.max(np.where(angs <= limits[1]))
        
        refl_max = np.max(refl_idx[llim+1:ulim+1])
        a_idx = np.where(refl_idx == refl_max)[0][0]
        
        ang = angs[a_idx]
         
    
    else:
        a_idx = (np.abs(refl[1:,0] - ang)).argmin() + 1
 
        refl_max = refl[a_idx, e_idx]
    
    return refl_max, ang
 
    
if __name__ == '__main__':

    treat_refl("B4C")
    refl = load_refl("B4C", "../../data/hom_refl/")
    get_refl(refl, 9.2, ang = 'max',  limits = [1.1e-03, 3.6e-03])
