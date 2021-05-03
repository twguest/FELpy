#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 09:55:11 2020

@author: twguest
"""

DOC_DIR = "/opt/FELpy/docs/coherent_focus_test/"
SLC_DIR = DOC_DIR + "/slices/"
PLT_DIR = DOC_DIR + "/plots/"

from tqdm import tqdm
import numpy as np
from wpg.wpg_uti_wf import plot_intensity_map

from felpy.model.beamlines.exfel_spb.exfel_spb import Instrument
from wpg.optical_elements import Drift
from wpg.optical_elements import Mirror_elliptical as MirEl
from felpy.model.core.beamline import Beamline
from felpy.model.tools import propagation_parameters
from felpy.model.core.wavefront import Wavefront
from felpy.model.src.coherent import construct_SA1_wavefront

wfr = construct_SA1_wavefront(512, 512, 9, 0.1)
plot_intensity_map(wfr)
spb = Instrument()
params = spb.params

NHE = MirEl(orient = params['NHE']["orientation"], p = 1, q = params['NHE']["distance to focus"],
        thetaE = params['NHE']["design angle"], theta0 = params['NHE']["incidence angle"],
        _x = params["NHE"]["xc"],
        _y = params["NHE"]["yc"],
        length = params['NHE']["length"],
        roll = params['NHE']["roll"],
        yaw = params['NHE']["yaw"],
        _refl = params['NHE']["reflectivity"],
        _ext_in = params['NHE']["_ext_in"], _ext_out = params['NHE']["_ext_out"]) 

 
NVE = MirEl(orient = params['NVE']["orientation"], p = 1, q = params['NVE']["distance to focus"],
        thetaE = params['NVE']["design angle"], theta0 = params['NVE']["incidence angle"],
        _x = params["NVE"]["xc"],
        _y = params["NVE"]["yc"],
        length = params['NVE']["length"],
        roll = params['NVE']["roll"],
        yaw = params['NVE']["yaw"],
        _refl = params['NVE']["reflectivity"],
        _ext_in = params['NVE']["_ext_in"], _ext_out = params['NVE']["_ext_out"]) 


drift = Drift(1.0)

bl = Beamline()
bl.append(NHE, propagation_parameters(1, 1, 1, 1))
bl.append(drift, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
bl.append(NVE, propagation_parameters(2, 1, 2, 1))
bl.append(Drift(2.), propagation_parameters(1/2, 1, 1/2, 1, mode = 'quadratic'))

bl.propagate(wfr)

plot_intensity_map(wfr)

SFILE = "/opt/FELpy/tmp/post_beamline"
wfr.store_hdf5(SFILE)
focus_distance = 0.2
slices = 25
dz = focus_distance/slices

beam_size = np.zeros(slices)

for itr in tqdm(range(slices)):
    
    wfr = Wavefront()
    wfr.load_hdf5(SFILE)
    
    bl = Beamline()
    bl.append(Drift(dz*(itr+1)), propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
    bl.propagate(wfr)
    
    plot_intensity_map(wfr)
    
    beam_size[itr] = wfr.get_fwhm()[0]

from matplotlib import pyplot as plt
plt.plot(beam_size)