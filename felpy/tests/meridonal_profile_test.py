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
from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.src.coherent import construct_SA1_wavefront
from felpy.model.tools import propagation_parameters
from wpg.srwlib import SRWLOptD as Drift
from felpy.utils.np_utils import get_mesh
import numpy as np
from felpy.utils.vis_utils import animate,colorbar_plot
from felpy.utils.os_utils import mkdir_p
from felpy.model.core.beamline import Beamline

def get_beamline(ekev):
    bl = get_beamline_object(ekev = ekev, crop = ["d1", "NVE"])
    return bl


def propThruFocus():
    
    mkdir_p(DOC_DIR)
    mkdir_p(PLT_DIR)
    
    nx, ny = 512, 512
    wfr = construct_SA1_wavefront(nx, ny, 5.0, 0.25)
    
    dof = 100e-06
    foc= 2.2
    nz = 100
    start  = foc-(dof/2)
    dz = dof/nz
    
    ekev = 5.0
    
    spb = get_beamline(ekev)
    spb.append(Drift(start), propagation_parameters(1,1,1,1, 'quadratic'))
    spb.propagate(wfr)
    
    mx = np.ones([nx, nz])
    my = np.ones([ny, nz])

 

    for itr in tqdm(range(nz)):
        
        bl= Beamline()
        bl.append(Drift(dz*(itr+1)), propagation_parameters(1,1,1,1, mode = 'fresnel'))
        bl.propagate(wfr)
        
# =============================================================================
#         wfr.plot(scale = "um",
#                  label = "$\Delta$z = {} m".format(-1.1 + itr*dz), 
#                  sdir = PLT_DIR + "{:04d}.png".format(itr))
#         
# =============================================================================
        mx[:,itr], my[:,itr] = wfr.get_profile_1d()
        
        
    px, py = wfr.get_spatial_resolution()
    mer_mesh_x = get_mesh(mx, px, dz)
    mer_mesh_y = get_mesh(my, py, dz)
    
    colorbar_plot(mx, mer_mesh_x,
                  xlabel = "x (m)",
                  ylabel = "z (m)",
                  clabel = "Intensity (a.u.)",
                  cmap = 'jet',
                  scale = 1,
                  normalise = True,
                  sdir = DOC_DIR + "x_mer.png")
    
    animate(PLT_DIR, DOC_DIR, "through_focus")
 
if __name__ == '__main__':
    propThruFocus()
 