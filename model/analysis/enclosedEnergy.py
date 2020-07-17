#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 12:35:15 2020

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

import os
import shutil
import numpy as np
from matplotlib.colors import LogNorm
import wpg.srwlib as srwlib
from model.src.coherent import coherentSource
from wpg.wpg_uti_wf import calc_pulse_energy, getOnAxisPowerDensity, getCentroid
from wpg.wavefront import Wavefront
from matplotlib import pyplot as plt
from matplotlib import ticker, cm
from wpg.generators import build_gauss_wavefront
from wpg import srwlib

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def create_circular_mask(nx, ny, c = None, r = None):
    """
    create a circular mask numpy array
    
    :param nx: number of pixels of array (horizontal) [int]
    :param ny: number of pixels of array (vertical) [int]
    :param c: center of array in pixels [list]from wpg import srwlib
    :param r: radius of array (in pixels) [int]
    
    :returns mask: binary mask [np array]
    """ 
    if c is None: # use the middle of the image
        c = (int(nx/2), int(ny/2))
    if r is None: # use the smallest distance between the center and image walls
        r = min(c[0], c[1], nx-c[0], ny-c[1])

    X, Y = np.ogrid[:nx, :ny]
    dist_from_center = np.sqrt((X - c[0])**2 + (Y-c[1])**2)

    mask = dist_from_center <= r
    
    return mask

def parseWavefront(wfr, itr = None):
    """
    Handler function for parsing wpg wavefront information to enclosed energy
    
    :param wfr: wpg wavefront object
    """
    
    if itr is None:
        ii = wfr.get_intensity().sum(axis = -1)  
    elif itr is not None:    
        ii = wfr.get_intensity()[:,:,itr]
    
    nx, ny = wfr.params.Mesh.nx, wfr.params.Mesh.ny
    
    ### ii = np.random.rand(nx,ny) ### DEBUG FOR RANDOM ARR. 17/07/20 TWG
    c = np.where(ii == np.amax(ii)) ### get position of maximum
    c = (c[0][0], c[1][0])
    
    return ii, nx, ny, c          

def finder(ii, nx, ny, c, efraction = 0.5, threshold = .01, iguess = None, nguesses = 5, fname = None, outdir = None):
    """
    Modified Bisection Algortihm, iterates over a set of guesses to reduce the 
    error function 1-(0.95/(ii-ii*mask(r))) in searching for a solution for r
    
    :param wfr: wpg wavefront structure
    :param efraction: desired energy fraction for enclosed energy, 0-1 [float]
    :param threshold: desired accuracy of algorithm [float] % of efraction
    :param iguess: initial radius guess [int]
    :param nguesses: number of initial guess (only activate when iguess = None) [int]
    
    :returns r: radius in pixelsi, nx, ny, c = parseWavefront(wfr)
    :returns error: error of radius
    """

    i = 0
    from wpg import srwlib
    
    ### make some random guesses
    if iguess is None:
        ins = np.random.randint(100, nx//2, nguesses)
        h = []
        for gs in ins:
            mask = create_circular_mask(nx, ny, c, r = gs)
            arr = mask*ii
            h.append(arr.sum()/ii.sum())
        
        gs = ins[np.argmin(h)]
        va = min(h)
        del h
    else:
        gs = iguess
        mask = create_circular_mask(nx, ny, c, r = gs)
        arr = mask*ii
        va = arr.sum()/ii.sum()[nx,ny]
        
        
    llim = efraction - threshold*efraction
    ulim = efraction + threshold*efraction 
    ### now do some bisection style loop
    print("Initial Guess: ", gs, "Initial Error:", va)
    while ulim  <= va or va <= llim and i < nx:
        
        if  va <= llim:
            gs = gs + gs//2
        elif va >= ulim:
            gs = gs - gs//2
            
        
        mask = create_circular_mask(nx, ny, c, r = gs)
        arr = mask*ii
        va = arr.sum()/ii.sum()
        
        i += 1

    r = gs
    err = abs(efraction-va)
    print("Enclosed Energy Radius: {}".format(r), va)
    print("Error = {:.4f}".format(err))
    print("After {} Iterations".format(i))
    
    
    return r, err

def plotEnclosed(ii, r, c, label = None, outdir = None, fname = "pulse"):
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot()
    
    im1 = ax1.imshow(ii, cmap = 'hot')
    circle = plt.Circle(c, r, color='w', fill=False)
    ax1.add_artist(circle)
    plt.axis("off")
    
    if label is not None:
        print()
        ax1.text(15, ii.shape[1]-15, label, c = 'w')
        
    if outdir is not None:
        mkdir_p(outdir)
        plt.savefig(outdir + str(fname) +".png")
    
    plt.show()

def getEnclosedEnergy(wfr, mode = 'integrated', outdir = None, fname = None, **kwargs):
    
    R = []
    if mode == 'integrated':
        ii, nx, ny, c = parseWavefront(wfr)
        r, err = finder(ii, nx, ny, c, **kwargs)
    
        R = (wfr.pixelsize()[0]*r, c, err)
        print(R*1e6, " $\mu m$")
 
        
        plotEnclosed(ii, r, c, outdir, fname)


    elif mode == 'pulse':
        n = wfr.get_intensity().shape[2]
        
        for itr in range(n):
            ii, nx, ny, c = parseWavefront(wfr, itr)

            r, err = finder(ii, nx, ny, c, **kwargs)
            
            R.append(wfr.pixelsize()[0]*r, c, err)
            print(R*1e6, " $\mu m$")
 
            srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')
            tax = np.linspace(wfr.params.Mesh.sliceMin, wfr.params.Mesh.sliceMax, wfr.params.Mesh.nSlices)
            srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
            plotEnclosed(ii, r, c, label = "{:.2f} fs".format(tax[itr]*1e15), outdir = outdir, fname = itr)
    
    np.save(oudir + fname, R)
    
def animate(indir, outdir, fname, delay = 0.1, rmdir = False):
    """
    create gif from pngs in directory
    """        
    os.system("convert -delay {} {}/*.png {}{}.gif".format(delay, indir, outdir, fname) )
    
    if rmdir == True:
        shutil.rmtree(indir)

if __name__ == '__main__':
    
    fname = sys.argv[1]    
        
    mkdir_p("/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/tmp/")
    mkdir_p("/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/enclosedEnergy/")

     wfr = Wavefront()
     wfr.load_hdf5("/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/out/{}".format(fname))
    
    getEnclosedEnergy(wfr, mode = 'pulse', efraction = 0.95, threshold = .01, outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/tmp/")
    animate(indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/tmp/", outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/enclosedEnergy/", fname = fname, delay = .1, rmdir = True)