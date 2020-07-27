#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:24:28 2020

@author: twguest
"""


#################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
####################################################

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema
from wpg.wavefront import Wavefront
import matplotlib.patches as mpatches

from wpg import srwlib

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def extractComplex(wfr):
    """
    extract the wavefront from wpg in the complex form
    
    :param wfr: wpg wfr structure
    """
    
    cmplx_wfr = wfr.toComplex()[0,:,:,:]    
    cmplx_wfr = np.moveaxis(cmplx_wfr, 2, 0)
    return cmplx_wfr

def getMC1d(wfr, temporal = False):
    """
    returns mutual coherence function in each of the transverse dimensions
    
    :param wfr: wpg wfr structure
    """
    
    if temporal == True:
        wfr = wfr.toComplex()[0,:,:,:]   
    else:
        wfr = extractComplex(wfr)
    
    xslc = wfr[:,wfr.shape[1]//2,:]
    yslc = wfr[wfr.shape[0]//2,:,:]

    Jx = np.dot(xslc.T.conjugate(), xslc) / wfr.shape[2]
    Jy = np.dot(yslc.T.conjugate(), yslc) / wfr.shape[2]
    
    return Jx,Jy

def getMC2d(wfr):
    """
    WARNING: SIZE OF J = [nx,ny]**2
    SUPER COMPUTATIONALLY EXPENSIVE
    
    returns mutual coherence function of the wavefield
    
    :param wfr: wpg wfr structure
    """
    
    wfr = extractComplex(wfr)
    
    
    slc = wfr[:,:,:]
        
    J = np.dot(slc.T.conjugate(), slc) / wfr.shape[2]
    return J

def getTDOC(Jx,Jy):
    """
    get transverse degree of coherence of the wavefront across each of the
    transverse dimensions slices
    """
    tdoc_x = np.diag(np.dot(Jx, Jx)).sum() / np.diag(Jx).sum()**2
    tdoc_y = np.diag(np.dot(Jy, Jy)).sum() / np.diag(Jy).sum()**2
    
    print("Horizontal TDOC: {}".format(tdoc_x.real))
    print("Vertical TDOC: {}".format(tdoc_y.real))
    
    return tdoc_x.real, tdoc_y.real

def calcCoherence(wfr, mode = 'temporal'):
    
    """
    calculate the temporal or spatial coherence of an FEL pulse
    
    :param wfr: wpg wfr structure
    :param mode: temporal or spatial [str]

    """
    
    
    if mode == 'temporal':
        Jx, Jy = getMC1d(wfr, temporal = True)
    
        axis_x = np.linspace(wfr.params.Mesh.sliceMin, wfr.params.Mesh.sliceMax, wfr.params.Mesh.nSlices) ### temporal axis
        axis_y = np.linspace(wfr.params.Mesh.sliceMin, wfr.params.Mesh.sliceMax, wfr.params.Mesh.nSlices)
    
    elif mode == 'spatial':
        Jx, Jy = getMC1d(wfr, temporal = False)
    
        axis_x = np.linspace(wfr.params.Mesh.xMin, wfr.params.Mesh.xMax, wfr.params.Mesh.nx) ### temporal axis
        axis_y = np.linspace(wfr.params.Mesh.yMin, wfr.params.Mesh.yMax, wfr.params.Mesh.ny) ### spatial axis
                             
    ii_x = np.abs(np.diag(Jx))
    ii_y = np.abs(np.diag(Jy))
    
    Jx /= ii_x**0.5 * ii_x[:, np.newaxis]**0.5
    Jy /= ii_y**0.5 * ii_y[:, np.newaxis]**0.5
    
    Jdx = np.abs(np.diag(np.fliplr(Jx)))
    Jdy = np.abs(np.diag(np.fliplr(Jy)))
    varI_x = (ii_x * axis_x**2).sum() / ii_x.sum()
    varI_y = (ii_y * axis_y**2).sum() / ii_y.sum()
    
    axisEx_x = axis_x*2
    axisEx_y = axis_y*2
    
    
    lmx = argrelextrema(Jdx, np.less)[0]  # local min
    lmy = argrelextrema(Jdy, np.less)[0]  # local min
    
    # for variance up to the 1st local minimum which is < 0.5:
    lmx = lmx[(axisEx_x[lmx] > 0) & (Jdx[lmx] < 0.50)]
    lmy = lmy[(axisEx_y[lmy] > 0) & (Jdx[lmy] < 0.50)]
    
    if len(lmx) > 0:
        cond = np.abs(axisEx_x) <= axisEx_x[lmx[0]]
        limJdx = axisEx_x[lmx[0]]
    else:
        cond = slice(None)  # for unconstrained variance calculation
        limJdx = None
    
    if len(lmy) > 0:
        cond = np.abs(axisEx_y) <= axisEx_y[lmy[0]]
        limJdy = axisEx_y[lmy[0]]
    else:
        cond = slice(None)  # for unconstrained variance calculation
        limJdy = None


    varJdx = (Jdx * axisEx_x**2)[cond].sum() / Jdx[cond].sum()
    varJdy = (Jdy * axisEx_y**2)[cond].sum() / Jdy[cond].sum()
    
 
    
    if mode == 'temporal':
# =============================================================================
#         fig = plt.figure()
#         ax1 = fig.add_subplot()
#         ax1.plot(axis_x*1e15, ii_x)
#         ax1.set_title("On-Axis Power Density")
#         ax1.set_ylabel("Intensity (W/$mm^2$)")
#         ax1.set_xlabel("Time (fs)")
# =============================================================================
        
        print("Hor. Temporal Coherence Length: {} fs".format(varJdx*1e15))
        print("Ver. Temporal Coherence Length: {} fs".format(varJdy*1e15))
    
    elif mode == 'spatial': 
        
# =============================================================================
#         fig, axs = plt.subplots(1,2, gridspec_kw={'hspace': 0.25, 'wspace': 0.25})
#         (ax1, ax2) = axs
#         
#         fig.suptitle('On-Axis Coherent Power Density')
#         
#         ax1.plot(axis_x*1e6, ii_x, color = 'r')
#         ax2.plot(axis_y*1e6, ii_y, color = 'b')
#         
#         ax1.set_ylabel("Intensity (W/$mm^2$)")
#         ax1.set_xlabel("Position ($\mu$m)")
#         
#         ax3 = ax1.twinx()
#         ax3.plot(axisEx_x*1e6, Jdx, color = 'r', linestyle = 'dashed')
#         ax3.set_yticks([])
#         ax4 = ax2.twinx()
#         ax4.plot(axisEx_y*1e6, Jdy, color = 'b', linestyle = 'dashed')
#         
#         rectx = mpatches.Rectangle((-limJdx, 0), width=2*limJdx, height=1.00, color='r', alpha=0.10)
#         ax3.add_patch(rectx)
#         coherence_analysis
#         recty = mpatches.Rectangle((-limJdy, 0), width=2*limJdy, height=1.00, color='b', alpha=0.10)
#         ax4.add_patch(recty)
#         
#         ax1.set_xlim([0, max(axis_x*1e6)])
#         ax2.set_xlim([0, max(axis_x*1e6)])
#         
# =============================================================================
        print("Hor. Coherence Length: {} $\mu$m".format(varJdx*1e6))
        print("Ver. Coherence Length: {} $\mu$m".format(varJdy*1e6))
    
    MCx = [Jx, axis_x]
    MCy = [Jy, axis_y]
    measurement = [varJdx, varJdy]
    
    return MCx, MCy, measurement 


def getCoherenceData(fname):
    
    wfr = Wavefront()
    wfr.load_hdf5(indir + fname)
    
    mcfdir = outdir + "MCF/"
    mcfx = mcfdir + "x/"
    mcfy = mcfdir + "y/"
    
    tdocDir = outdir + "TDOC/"
  
    cTimeDir = outdir + "cTime/"
    cLenDir = outdir + "cLen/"
    
    mkdir_p(mcfdir)
    mkdir_p(mcfx)
    mkdir_p(mcfy)
    mkdir_p(tdocDir)
    mkdir_p(cTimeDir)
    mkdir_p(cLenDir)
    
    Jx, Jy = getMC1d(wfr)
    tdoc = getTDOC(Jx, Jy)
    
    MCx, MCy, Ctime = calcCoherence(wfr, mode = 'temporal')
    MCx, MCy, Clen = calcCoherence(wfr, mode = 'spatial')
    
    np.save(mcfx + fname, Jx)
    np.save(mcfy + fname, Jy)
    
    np.save(tdocDir + fname, tdoc)
    np.save(cTimeDir + fname, Ctime)
    np.save(cLenDir + fname, Clen)
    
if __name__ == "__main__":
    
# =============================================================================
#     import os
#     
#     indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/out/"
#     outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/coherence/"
#     
#     mkdir_p(outdir)
#     
#     fname = sys.argv[1]    
#         
#     getCoherenceData(fname)
# =============================================================================

    from model.tools import constructPulse
    
    wfr = constructPulse()
    calcCoherence(wfr, mode = 'temporal')