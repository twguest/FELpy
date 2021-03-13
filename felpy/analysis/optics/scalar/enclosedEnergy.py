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
import time
import os
import shutil
import numpy as np
from matplotlib.colors import LogNorm
import wpg.srwlib as srwlib
from model.src.coherent import construct_SA1_wavefront
from wpg.wpg_uti_wf import calc_pulse_energy, get_axial_power_density, get_centroid
from felpy.model.core.wavefront import Wavefront
from matplotlib import pyplot as plt
from matplotlib import ticker, cm
from model.tools import constructPulse, argmax2d
from wpg import srwlib
from tqdm import tqdm
from model.tools import create_circular_mask

from wpg.wpg_uti_wf import get_axis

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)



def finder(ii, nx, ny, c, efraction = 0.5, threshold = .01, iguess = None, nguesses = 10, fname = None, outdir = None, VERBOSE = True):
    """
    Modified Bisection Algortihm, iterates over a set of guesses to reduce the 
    error function 1-(0.95/(ii-ii*mask(r))) in searching for a solution for r
    
    :param wfr: wpg wavefront structure
    :param efraction: desired energy fraction for enclosed energy, 0-1 [float]
    :param threshold: desired accuracy of algorithm [float] % of efraction
    :param iguess: initial radius guess [int]
    :param nguesses: number of initial guess (only activate when iguess = None) [int]
    
    :returns r: radius in pixels
    :returns error: error of radius
    """

    i = 0
    llim = efraction - threshold*efraction
    ulim = efraction + threshold*efraction 
    
    ### make some random guesses
    if iguess is None:
        ins = np.random.randint(0, nx, nguesses)
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
        
        

    ### now do some bisection style loop
    
    if VERBOSE:
        print("\nInitial Guess: ", gs, "Initial Error:", va)
    
    
    while (ulim  <= va or va <= llim) and i < nx:
        
        if  va <= llim:
            gs = gs + 1#gs//2
        elif va >= ulim:
            gs = gs - 1#gs//2
            
        
        mask = create_circular_mask(nx, ny, c, r = gs)
        arr = mask*ii
        va = arr.sum()/ii.sum()

        i += 1

    r = gs
    err = abs(efraction-va)
    
    if VERBOSE:
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
        ax1.text(15, ii.shape[1]-15, label, c = 'w')
        
    if outdir is not None:
        
        if mode == 'pulse':
            mkdir_p(outdir + "tmp/")
            plt.savefig(outdir + "/tmp/{}/{}.png".format(fname, itr))
        elif mode == 'integrated':
            plt.savefig(outdir + "{}.png".format(fname))
            plt.show()



def getEnclosedEnergy(wfr, mode = 'integrated', outdir = None, fname = None, nSlc = 20, VERBOSE = True, bPlot = False, **kwargs):
    
    dx, dy = wfr.get_spatial_resolution()
    
    if mode == 'integrated':
        
        ii = wfr.get_intensity().sum(axis = -1)  
        

        nx, ny = ii.shape
        c = argmax2d(ii)
        
        results, err = finder(ii, nx, ny, c, VERBOSE = VERBOSE, **kwargs)
        results = (float(results)*dx, float(results*dy))
        
        if bPlot:
            if outdir is None:
                raise(ValueError, "Output Directory Not Specified")
            plotEnclosed(ii, results, c, outdir = outdir, fname = fname, label = fname)
        

    elif mode == 'pulse':
        
        timeAx = get_axis(wfr, 't')


        ii = wfr.get_intensity()
        sz0 = ii.shape
        
        idx = np.flip(np.argsort(ii.sum(axis = 0).sum(axis = 0)))[:nSlc]
        
        ii = ii[:,:,idx]
        
        sz1 = ii.shape
        
        ### MOVED OUTSIDE VERBOSE
        print("Wavefront Reduced from {} to {}".format(sz0, sz1))
            
        nx, ny, nz = ii.shape
        
        results = np.zeros([sz0[-1], 4])
        results[:,0] = np.arange(sz0[-1])
        
        for itr in tqdm(range(nz)):
            

            slc = ii[:,:,itr]
            c = argmax2d(slc)
            R, err = finder(slc, nx, ny, c, VERBOSE = VERBOSE, **kwargs)
            
            
           
            
            results[idx[itr], 1] = timeAx[idx[itr]]
            results[idx[itr], 2] = float(R)*dx
            results[idx[itr], 3] = float(R)*dy
            
            if bPlot:
                
                if outdir is None:
                    raise(ValueError, "Output Directory Not Specified")
                plotEnclosed(slc, results, c, mode = 'pulse',
                             label = "{:.2f} fs".format(timeAx[idx[itr]]*1e15), outdir = outdir, fname = fname, itr = itr)
        
    
        
    return results
    

def animate(indir, outdir, fname, delay = 0.1, rmdir = False):
    """
    create gif from pngs in directory
    """        
    os.system("convert -delay {} {}/*.png {}{}.gif".format(delay, indir, outdir, fname) )
    
    if rmdir == True:
        shutil.rmtree(indir)


def speedTest():
    nz = 6
    wfr = constructPulse(1024, 1024, nz = nz, tau = 1e-12)
    
    start = time.time()
    getEnclosedEnergy(wfr, mode = 'pulse', VERBOSE = True)
    stop = time.time()
    time1 = stop-start
    print("Time Taken per Slice: {} s".format(time1/nz))



def run(wfr, mode, nSlc = 250, VERBOSE = True):
    
    results = getEnclosedEnergy(wfr, mode = mode, nSlc = 250, VERBOSE = VERBOSE)
    return results

if __name__ == '__main__':
    
    nz = 6
    wfr = constructPulse(1024, 1024, nz = nz, tau = 1e-12)
    
    results = run(wfr, mode = 'pulse')