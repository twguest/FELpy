#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 15:18:05 2020

@author: twguest
"""
import os
import sys
sys.path.append("../../")
from felpy.model.source.coherent import construct_SA1_wavefront
from felpy.model.tools import mkdir_p
from felpy.model.wavefront import Wavefront
import numpy as np
from matplotlib import pyplot as plt

import shutil

def get_pointing_vector(wfr, propagated_wfr, px, py, pz, sdir = None):
    """
    Calculates the variance in pointing vector of a wavefield
    
    :param wfr: complex wfr structure
    :param mode: 'slice' or 'integrated':
        if slice - returns array of pointing angles of dim [nx,ny,nt]
        if 'integrated' - returns array of pointing angles of dim [nx,ny]
    
    :returns: pvx x-component of pointing vector from slice-to-slice
    :returns: pvy y-component of pointing vector from slice-to-slice
    :returns: pvz z-component of pointing vector from slice-to-slice
    """

    x_gradient = np.gradient(wfr, px, axis = 0)
    y_gradient = np.gradient(wfr, px, axis = 1)
    print(x_gradient.shape)
    z_gradient = np.gradient(np.moveaxis([wfr, propagated_wfr], 0, -1), pz)[-1][:,:,0]
    
    
    jx = np.arctan(x_gradient/z_gradient)
    jy = np.arctan(y_gradient/z_gradient)
    
    plt.imshow(jy.real)
    return jx, jy
    
 
# =============================================================================
# def plot_pointing_angle(vdx, vdy, vdz):
#     
#     fig, axs = plt.subplots(1,3, figsize = (18, 6), dpi = 1020)
#     
#     for ax in axs:
#         ax.set_aspect('equal')
#         ax.set_xticks([])
#         ax.set_yticks([])
#     
#     
#     
#     im = axs[0].imshow(np.var(vdx, axis = -1), vmin = 0, vmax = np.pi)
#     im = axs[1].imshow(np.var(vdy, axis = -1), vmin = 0, vmax = np.pi)
#     im = axs[2].imshow(np.var(vdz, axis = -1), vmin = 0, vmax = np.pi)
#     
#     cbar = fig.colorbar(im, ax=axs.ravel().tolist(), fraction = 0.0725, aspect = 3.5)
#     cbar.set_ticks([0, np.pi/2, np.pi])
#     cbar.set_ticklabels( ["0", r"$\frac{\pi}{2}$", "$\pi$"])
#     cbar.ax.tick_params(labelsize=25)  # set your label size here
#     #cbar.ax.set_aspect(2.5)
# =============================================================================
    
if __name__ == '__main__':
    x = -np.linspace(-400, 400, 400)
    y = np.linspace(-400, 400, 400)
    grid = np.meshgrid(x,y)
    px = 1
    py = 1
    wav = 1e-10
    z = 1e-08
    from felpy.model.source.partially_coherent import define_wfr_tilt
    from felpy.model.fresnel_propagator import frensel_propagator
    
    wfr = define_wfr_tilt(grid, kx = 0, ky = 1)
    
    propagated_wfr = frensel_propagator(wfr, px, py, wav, z) 
    
    jx, jy = get_pointing_vector(wfr, propagated_wfr, px, py, pz = z)