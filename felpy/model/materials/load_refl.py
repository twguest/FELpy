#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.1"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""

import time

import numpy as np
import pandas as pd

def treat_refl(material = 'B4C', indir = None):
    """
    Used for converting excel data obtained via henke to dat file
    
    :param material: mirror material which data is defined for
    """
    
    if indir is None:
        excel = pd.read_excel("../../data/spb/hom_refl/{}_refl.xlsx".format(material),header = None)
    else:
        excel = pd.read_excel(indir + "{}_refl.xlsx".format(material),header = None)
    df = pd.DataFrame(excel)
    
    refl = np.asarray(df.values)
    
    if indir is None: 
        np.save("data/hom_refl/spb/{}_refl.npy".format(material), refl)
    else:
        np.save(indir + "{}_refl.npy".format(material), refl)
    
def load_refl(material = 'B4C', indir = None):
    """
    load data from .npy file generate by treat_refl method
    
    :param material: mirror material which data is defined for
    """

    if indir is None:
        import os
        fpath = os.path.dirname(os.path.realpath(__file__)) ## function path
        fpath = os.path.join(fpath, "../../data/spb/hom_refl/")
        refl = np.load(fpath +"{}_refl.npy".format(material))
    else:
        refl = np.load(indir + "/{}_refl.npy".format(material))
    return refl

def get_refl(refl, ekev, ang):
    
    e_idx = (np.abs(refl[0,1:] -ekev*1e3)).argmin() + 1
    refl_idx = refl[:, e_idx][1: ]
    
    angs = refl[1:,0]

    a_idx = (np.abs(refl[1:,0] - ang)).argmin() + 1
 
    refl_max = refl[a_idx, e_idx]
    
    return refl_max
 
    
if __name__ == '__main__':

    treat_refl("B4C")
    refl = load_refl("B4C", "../../data/spb/hom_refl/")
    get_refl(refl, 9.2, ang = 'max',  limits = [1.1e-03, 3.6e-03])
