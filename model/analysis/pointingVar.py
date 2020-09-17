#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 15:18:05 2020

@author: twguest
"""
import os
import sys
sys.path.append("../../")
from model.tools import constructPulse, generateTestPulses
from wpg.wavefront import Wavefront
import numpy as np
from matplotlib import pyplot as plt


def pointingVarIntra(wfr):
    """
    intra-pulse
    """

    wfr = wfr.toComplex()[0,:,:,:]
    
    dx = np.zeros(wfr.shape)
    dy = np.zeros(wfr.shape)
    
    
    for t in range(wfr.shape[2]):
        
        wfr[:,:,t] *= np.random.rand(512,512)*1j*10
        
        jx, jy = (np.matmul(wfr[:,:,t],np.gradient(wfr[:,:,t].conjugate()))
                  -np.matmul(np.gradient(wfr[:,:,t]),wfr[:,:,t].conjugate()))
        
        jx *= 1j
        jy *= 1j
        
        
        dx[:,:,t] = np.arctan(jx.imag/jx.real)
        dy[:,:,t] = np.arctan(jy.imag/jy.real)
        
        plt.imshow(dx[:,:,t])
        plt.show()
        

    return vdx, vdy, vdz


def pointingVarInter(wfrdir):
    """ 
    inter-pulse pointing vector variance
    """
    
    N = len(os.listdir(wfrdir))
    
    wfr = Wavefront()
    wfr.load_hdf5(wfrdir + os.listdir(wfrdir)[0])
    wfr = wfr.toComplex()[0,:,:,:].sum(axis = -1)
    
    dx = np.zeros([*wfr.shape, N])
    dy = np.zeros([*wfr.shape, N])

    for n in range(N):
        
        wfr = Wavefront()
        wfr.load_hdf5(wfrdir + os.listdir(wfrdir)[n])
        wfr = wfr.toComplex()[0,:,:,:].sum(axis = -1)
        
        jx, jy = (np.matmul(wfr,np.gradient(wfr.conjugate()))
                  -np.matmul(np.gradient(wfr),wfr.conjugate()))
            
        jx *= 1j
        jy *= 1j
    
                
        dx[:,:,n] = np.arctan(jx.imag/jx.real)
        dy[:,:,n] = np.arctan(jy.imag/jy.real)
        

    vdx = np.var(dx, axis = -1)
    vdy = np.var(dy, axis = -1)
    vdz = np.sqrt(vdx**2 + vdy**2)
        
    return vdx, vdy, vdz


def testpointingVar(mode = 'intra'):
    
    if mode == 'intra':
        wfr = constructPulse(512, 512, 5)
        vdx, vdy, vdz = pointingVarInter(wfr)
        
    elif mode == 'inter':
        sdir = "../../data/testPulses/"
        generateTestPulses(sdir, N = 5)
        vdx, vdy, vdz = pointingVarIntra(sdir)    
    
    return vdx, vdy

def plotPointingAngle(vdx, vdy, vdz):
    
    fig, axs = plt.subplots(1,3, figsize = (18, 6), dpi = 1020)
    
    for ax in axs:
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
    
    
    
    im = axs[0].imshow(vdx, vmin = 0, vmax = np.pi)
    im = axs[1].imshow(vdy, vmin = 0, vmax = np.pi)
    im = axs[2].imshow(vdz, vmin = 0, vmax = np.pi)
    
    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), fraction = 0.0725, aspect = 3.5)
    cbar.set_ticks([0, np.pi/2, np.pi])
    cbar.set_ticklabels( ["0", r"$\frac{\pi}{2}$", "$\pi$"])
    cbar.ax.tick_params(labelsize=25)  # set your label size here
    #cbar.ax.set_aspect(2.5)
if __name__ == '__main__':
    
    sdir = "../../data/testPulses/"
    generateTestPulses(sdir, nx = 512, ny = 512, N = 5)
    vdx, vdy, vdz = pointingVarInter(sdir)    

    plotPointingAngle(vdx,vdy,vdz)
    