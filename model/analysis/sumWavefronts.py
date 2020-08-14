#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 12:22:33 2020

@author: twguest
"""

import numpy as np
from os import listdir
def sumWavefronts(indir, outfile, mode = 'both'):
    
    if mode == 'complex':
        pass
    if mode == 'intensity':
        pass
    if mode == 'both':
        
        for f in listdir(indir):
            wfr = Wavefront()