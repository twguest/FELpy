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
import os
import numpy as np
from time import time
from wpg import srwlib
from model.tools import radial_profile, binArray
from model.tools import constructPulse ## for testings
from wpg.wpg_uti_wf import getAxis
from wpg.wavefront import Wavefront 
from tqdm import tqdm




def coherenceTime(wfr, tStep, bins = 1, VERBOSE = True):
    """
    Calculate the coherence time of complex wavefield of shape
    [nx, ny, nt].
    
    ref: Coherence properties of the radiation from X-ray free electron laser
    
    :param wfr: complex wavefield
    :param tstep: temporal step between slices
    
    :returns tau: coherence time [s]
    """
    
    nz0 = wfr.shape[-1]
    
    if bins == 1:
        pass
    else:
        wfr = binArray(wfr, axis = -1, binstep = nz0//bins, binsize = 1 )
        
    nz1 = wfr.shape[-1]
    tStep *= (nz0/nz1)
    
    b = np.zeros([*wfr.shape])

    for i in tqdm(range(wfr.shape[-1])):
        
        A = np.roll(wfr, -i, axis = 2)
        B = np.repeat(wfr[:,:,i][:, :, np.newaxis], wfr.shape[-1], axis=-1)
    
        ## DEGUB print(A[:,:,0] == wfr[:,:,i])
        ## DEBUG print([B[:,:,k] == wfr[:,:,i] for k in range(wfr.shape[-1])])
        
        b[:,:,i] = ((A*B.conjugate()).mean(axis = -1))/np.sqrt(
            (abs(A)**2).mean(axis=-1)*(abs(B)**2).mean(axis = -1))

    
    tau = (abs(b)**2).sum(axis = -1)[0,0]
    
    if VERBOSE:
        print("Time Step: {}".format(tStep))
        print("Coherence Time: {:.2e} fs".format(tau*1e15*tStep))
        
    return tau*tStep


def coherenceLen(wfr, dx, dy, VERBOSE = True):
    """
    Calculate coherence length of a complex wavefield of shape
    [nx, ny. nz]
    
    :param wfr: complex wavefield
    :param dx: horizontal pixel size
    :param dy: vertical pixel size
    
    :returns Jd: complex degree of coherence
    :returns clen: coherence length [m]
    """
    
    
    profile, r = complexRadialProfile(wfr)

    
    nt = wfr.shape[-1]
    
    J = np.dot(profile, profile.T.conjugate())/ nt
    II = np.abs(np.diag(J))  # intensity as the main diagonal
    
    J /= II**0.5 * II[:, np.newaxis]**0.5
    Jd = np.abs(np.diag(np.fliplr(J)))  # DoC as the cross-diagonal
    
    lm = np.arange(Jd.shape[0])

    lm = lm[(lm >= Jd.shape[0]//2) & (Jd[lm] < 0.5)]

    rstep = np.sqrt((dx)**2 + (dy)**2)

    
    try:
        lm = lm[0] - Jd.shape[0]//2 
    except(IndexError):
        lm = np.inf
    
    clen = lm*rstep
    
    if VERBOSE: 
        print("Coherence Length: {:.2f} um".format(clen*1e6))
    return clen

def transverseDOC(wfr, VERBOSE = True):
    """
    get transverse degree of coherence of the wavefront across each of the
    transverse dimensions slices
    """
    
    p, r =  complexRadialProfile(wfr)
    nt = wfr.shape[-1]
    J = np.dot(p, p.T.conjugate())/nt
    
    
    tdoc = np.diag(np.dot(J, J)).sum() / np.diag(J).sum()**2
    
    if VERBOSE:
        print("Transverse Degree of Coherence: {:.4f}".format(tdoc.real))
    
    return tdoc

def complexRadialProfile(wfr):
    """
    Calculate the radial profile of a complex array by azimuthal averaging:
    
        I_{radial}(R) = \int_0^R \frac{I(r)2\pi r}{\pi R^2} dr
    
    :param wfr: complex wavefield [np array]
    
    :returns prof: radial profile
    """
        
    r = radial_profile(wfr[:,:,0].real, [wfr.shape[0]//2,wfr.shape[1]//2])[1]
    
    r = np.diag(r).copy()
    r[:r.shape[0]//2] *= -1
    
    rp = np.stack([radial_profile(wfr[:,:,i].real,
                                  [wfr.shape[0]//2,wfr.shape[1]//2])[0]
                   + radial_profile(wfr[:,:,i].imag,
                                    [wfr.shape[0]//2,wfr.shape[1]//2])[0]*1j
                   for i in range(wfr.shape[-1])])
    
    prof = np.moveaxis(rp, 0, -1), r
    return prof
  
    
def speedTest(nx = 1024, ny = 1024, nz = 5, N = 2000):
    """
    Estimate the approximate completion time for measuremetn of the 
    coherence time of a wavefront of size [nx, ny, nz]
    
    :param nx: wfr shape [0]
    :param ny: wfr shape [1]
    :param ny: wfr shape [2]
    :paran N: desired number of slices in full-scale test
    """
    
    
    print("Performing Coherence Time Speed Test\n")
 
    
    for t in range(2, nz+1):
        
        wfr = constructPulse(nx = nx, ny = ny, nz = t, tau = 0.5e-25)
        srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')
        wfr = wfr.toComplex()[0,:,:,:] 

        start = time()
        coherenceTime(wfr, tStep = 1, VERBOSE = False)
        fin = time()
        
        print("Time Per-Slice in minutes for {} Slices: {}".format(t, (fin-start)/60))


def MSS(arr1, arr2):
    """ 
    calculate the mean-squared similarity (1-error) of two arrays
    
    :param arr1: expected array
    :param arr2: measured array
    
    :returns sim: similarity metric [0-1]
    """
    # the 'Mean Squared Similarity' between the two images is the
    # sum of the squared difference between the two images;
    # NOTE: the two images must have the same dimension
    err = np.sum((arr1.astype("float") - arr2.astype("float")) ** 2)
    err /= float(arr2.shape[0] * arr1.shape[1])
    sim = 1-err
    # return the MSE, the lower the error, the more "similar"
    # the two images are
    return sim


def testUsage():
    
    wfr = constructPulse(nx = 500, ny = 500, nz = 15, tau = 0.5e-05)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 't')

    tstep = getAxis(wfr, axis = 't')
    tstep = tstep[1]-tstep[0]
    
    xstep, ystep = wfr.pixelsize()
    
    wfr = wfr.toComplex()[0,:,:,:]
    wfr += np.random.rand(*wfr.shape)*1j*100
    
    for b in range(1, 10):
        print("bins {}".format(b))
        tau = coherenceTime(wfr, tstep, bins = b)
        
    clen = coherenceLen(wfr, xstep, ystep)
    tdoc = transverseDOC(wfr)

    del tau, clen, tdoc
    
    
def run(wfr):
    
    tstep = getAxis(wfr, axis = 't')
    tstep = tstep[1]-tstep[0]
    
    xstep, ystep = wfr.pixelsize()
    
    wfr = wfr.toComplex()[0,:,:,:]
    
    
    tau = coherenceTime(wfr, tstep, VERBOSE=True)
    clen = coherenceLen(wfr, xstep, ystep, VERBOSE=True)
    tdoc = transverseDOC(wfr, VERBOSE=True)
    
    return tau, clen, tdoc

def binTest():
    
    wfr = Wavefront()
    wfr.load_hdf5("")
    
    tstep = getAxis(wfr, axis = 't')
    tstep = tstep[1]-tstep[0]
    
    for bins in np.arange(10, 1, -1):
        tau = coherenceTime(wfr, tstep, VERBOSE = True)
        print("{} bins: {}".format(bins, tau))
    

def launch():
    pass    


if __name__ == "__main__":
    
    binTest()
    #testUsage()
    #speedTest(nz = 15)