#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:17:32 2020

@author: twguest
"""

#################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################
from matplotlib import pyplot as plt
from poppy import zernike
from model.src.coherent import construct_SA1_wavefront
import numpy as np
from wpg.wpg_uti_wf import calculate_fwhm

 

wfr = construct_SA1_wavefront(1000, 1000, 12, 0.1)
ii = wfr.get_intensity()[:,:,0]
ph = wfr.get_phase()[:,:,0]


aperture = np.ones(ii.shape)
aperture[np.where(ii < ii.max()/10)] = 0
zc = zernike.opd_expand(ph, aperture = aperture)

opd = zernike.opd_from_zernikes(zc)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.bar(np.linspace(0,len(zc), len(zc)),zc)
ax1.set_xticks(np.linspace(0+1,len(zc)+1, len(zc)+1))
ax1.set_xticklabels(labels = ["$Z_{}$".format(int(z)) for z in np.linspace(0+1,len(zc)+1, len(zc)+1)])

plt.show()

for i in range(15):
    if i == 0:
        pass
    else:
        print("index: {},".format(i-1), "j = {}".format(i), zernike.zern_name(i))