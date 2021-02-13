#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 15:18:05 2020

@author: twguest
"""
import os
import sys
sys.path.append("../../")
from felpy.model.src.coherent import construct_SA1_wavefront
from felpy.model.tools import mkdir_p
from wpg.wavefront import Wavefront
import numpy as np
from matplotlib import pyplot as plt

import shutil

def get_pointing_vector(wfr, sdir = None):
    """
    Calculates the variance in pointing vector of a wavefield
    
    :param wfr: wpg wfr structure
    :param mode: 'slice' or 'integrated':
        if slice - returns array of pointing angles of dim [nx,ny,nt]
        if 'integrated' - returns array of pointing angles of dim [nx,ny]
    
    :returns: pvx x-component of pointing vector from slice-to-slice
    :returns: pvy y-component of pointing vector from slice-to-slice
    :returns: pvz z-component of pointing vector from slice-to-slice
    """


    pvx = np.zeros(wfr.shape)
    pvy = np.zeros(wfr.shape)
    pvz = np.zeros(wfr.shape)

                
    jx, jy = (np.matmul(wfr,np.gradient(wfr.conjugate()))
              -np.matmul(np.gradient(wfr),wfr.conjugate()))
        
    jx *= 1j
    jy *= 1j

            
    pvx[:,:] = np.arctan(jx.imag/jx.real)
    pvx[:,:] = np.arctan(jy.imag/jy.real)
    pvz[:,:] = np.sqrt(pvx**2 + pvy**2)
    
    pvx[np.where(pvx == np.nan)] = 0
    pvy[np.where(pvy == np.nan)] = 0
    pvz[np.where(pvz == np.nan)] = 0
        
    return pvx, pvx, pvz


def test_pointing_vector(mode = 'intra'):
    
    if mode == 'intra':
        wfr = constructslice(512, 512, 5)
        vdx, vdy, vdz = pointingVarInter(wfr)
        
    elif mode == 'inter':
        sdir = "../../data/testslices/"
        generateTestslices(sdir, N = 5)
        vdx, vdy, vdz = pointingVarIntra(sdir)    
    
    return vdx, vdy

def plot_pointing_angle(vdx, vdy, vdz):
    
    fig, axs = plt.subplots(1,3, figsize = (18, 6), dpi = 1020)
    
    for ax in axs:
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
    
    
    
    im = axs[0].imshow(np.var(vdx, axis = -1), vmin = 0, vmax = np.pi)
    im = axs[1].imshow(np.var(vdy, axis = -1), vmin = 0, vmax = np.pi)
    im = axs[2].imshow(np.var(vdz, axis = -1), vmin = 0, vmax = np.pi)
    
    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), fraction = 0.0725, aspect = 3.5)
    cbar.set_ticks([0, np.pi/2, np.pi])
    cbar.set_ticklabels( ["0", r"$\frac{\pi}{2}$", "$\pi$"])
    cbar.ax.tick_params(labelsize=25)  # set your label size here
    #cbar.ax.set_aspect(2.5)
    
if __name__ == '__main__':
    
    sdir = "../../data/testslices/"
    generateTestslices(sdir, nx = 512, ny = 512, N = 5)
    ddir = "../../data/pointing/"
    mkdir_p(ddir)
    for n in range(len(os.listdir(sdir))):
        
        wfr = Wavefront()
        wfr.load_hdf5(sdir + os.listdir(sdir)[n])
        px, py, pz = get_pointing_vector(wfr, mode = 'slice', sdir = ddir, ID = n)
        
    shutil.rmtree(sdir)