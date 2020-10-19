# -*- coding: utf-8 -*-

import numpy as np

def mean_intensity(ii, mode = 'train'):
    """
    get the mean intensity of the detector data
    """
    
    if mode == 'all':
        mi = np.mean(ii, axis = -1).mean(axis = -1)
    if mode == 'train':
        mi = np.mean(ii, axis = -1)
    if mode == 'pulse':
        mi = np.mean(ii, axis = -2)
    if mode == 'shape':
        mi = ii.mean(axis = 0).mean(axis = 0).mean(axis = -1)
    return mi