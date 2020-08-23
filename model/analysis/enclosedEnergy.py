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
from model.tools import constructPulse
from wpg import srwlib

from model.tools import create_circular_mask

from wpg.wpg_uti_wf import getAxis

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)



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
    c = (c[1][0], c[0][0])
    
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

def plotEnclosed(ii, r, c, mode = 'integrated', label = None, outdir = None, fname = "pulse", itr = None):
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot()
    
    ax1.imshow(ii, cmap = 'hot')
    circle = plt.Circle(c, r, color='w', fill=False)
    ax1.add_artist(circle)
    plt.axis("off")
    ax1.plot(c[0],c[1], marker = 'x', color = 'k')
    if label is not None:
        print()
        ax1.text(15, ii.shape[1]-15, label, c = 'w')
        
    if outdir is not None:
        
        if mode == 'pulse':
            mkdir_p(outdir + "tmp/")
            plt.savefig(outdir + "/tmp/{}/{}.png".format(fname, itr))
        elif mode == 'integrated':
            plt.savefig(outdir + "{}.png".format(fname))
            plt.show()

def getEnclosedEnergy(wfr, mode = 'integrated', outdir = None, fname = None, **kwargs):
    
    R = []
    if mode == 'integrated':
        ii, nx, ny, c = parseWavefront(wfr)
        r, err = finder(ii, nx, ny, c, **kwargs)
    
        R = [fname, wfr.pixelsize()[0]*r, c[0], c[1], err]
        print(R[1]*1e6, " um")
 
        
        plotEnclosed(ii, r, c, outdir = outdir, fname = fname, label = fname)


    elif mode == 'pulse':
        n = wfr.get_intensity().shape[2]
        
        for itr in range(n):
            ii, nx, ny, c = parseWavefront(wfr, itr)

            r, err = finder(ii, nx, ny, c, **kwargs)
            
            R.append([fname, itr, wfr.pixelsize()[0]*r, c[0], c[1], err])
 
            srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')
            tax = np.linspace(wfr.params.Mesh.sliceMin, wfr.params.Mesh.sliceMax, wfr.params.Mesh.nSlices)
            srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
            
            if mode == 'integrated':
                plotEnclosed(ii, r, c, mode = 'integrated', label = None, outdir = outdir, fname = fname, itr = None)
            elif mode == 'pulse':
                plotEnclosed(ii, r, c, mode = 'pulse',
                             label = "{:.2f} fs".format(tax[itr]*1e15), outdir = outdir, fname = fname, itr = itr)
    
    np.save(outdir + fname, R, allow_pickle= True)
    

def animate(indir, outdir, fname, delay = 0.1, rmdir = False):
    """
    create gif from pngs in directory
    """        
    os.system("convert -delay {} {}/*.png {}{}.gif".format(delay, indir, outdir, fname) )
    
    if rmdir == True:
        shutil.rmtree(indir)


def speedTest():
    
    wfr = constructPulse(500, 500, 10)
    

    if __name__ == '__main__':
    
    fname = sys.argv[1]    
    mode = sys.argv[2]
    #mode = 'integrated'
    #fname = "NanoKB-Pulse_54.h5"    
    outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/enclosedEnergy/SliceBySlice/"
    
    
    if mode == 'pulse':
        tmpdir = outdir + "/tmp/"
        ftmp = tmpdir + "/fname/"
        mkdir_p(tmpdir)
        mkdir_p(ftmp)

    mkdir_p(outdir)
    
    wfr = Wavefront()
    wfr.load_hdf5("/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/out/{}".format(fname))
    
    getEnclosedEnergy(wfr, mode = mode, efraction = 0.50, threshold = .01, outdir = outdir, fname = fname)
    
    if mode == 'pulse':
        animate(indir = ftmp, outdir = outdir, fname = fname, delay = .1, rmdir = True)