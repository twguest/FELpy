#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:24:28 2020

@author: twguest
"""



import os
import numpy as np
from time import time
from wpg import srwlib
from felpy.model.tools import radial_profile, binArray
from felpy.model.src.coherent import construct_SA1_pulse ## for testings
from wpg.wpg_uti_wf import get_axis
from wpg.wavefront import Wavefront 
from tqdm import tqdm
from felpy.utils.job_utils import JobScheduler
#import wpg.srwlib as srwl
from wpg import srwlpy as srwl
from felpy.utils.np_utils import memory_map, readMap
import multiprocessing as mp
from functools import partial




def get_longitudinal_coherence(slice_no, cfr, map_loc = None, bins = 1, VERBOSE = True):
    """
    Calculate the longitudinal correlation of each slice of a complex wavefront
    of shape [nx, ny, nz] against a single slice of shape [nx,ny] at longitudinal
    interval defined by the slice_no
    
    :param cfr: complex wavefield
    :param slice_no: longitudinal index [int]
    
    :returns g: complex degree of coherence
    """
    
    A = np.roll(cfr, -slice_no, axis = 2)
    B = np.repeat(cfr[:,:,slice_no][:, :, np.newaxis], cfr.shape[-1], axis=-1)

    ## DEGUB print(A[:,:,0] == wfr[:,:,i])
    ## DEBUG print([B[:,:,k] == wfr[:,:,i] for k in range(wfr.shape[-1])])
    
  
    if map_loc is not None:
        
        mmap = memory_map(map_loc,
                  shape = cfr.shape,
                  dtype = 'complex64')
            
        mmap[:,:,slice_no] = ((A*B.conjugate()).mean(axis = -1))/np.sqrt(
        (abs(A)**2).mean(axis=-1)*(abs(B)**2).mean(axis = -1))

        
    else:
        return ((A*B.conjugate()).mean(axis = -1))/np.sqrt(
        (abs(A)**2).mean(axis=-1)*(abs(B)**2).mean(axis = -1))



def get_coherence_time(cfr, tStep, mpi = False, map_loc = "/tmp/coherence_map",
                       bins = 1, VERBOSE = True):
    """
    Calculate the coherence time of complex wavefield of shape
    [nx, ny, nt].
    
    Relevant for statistically stationary sources. 
    
    ref: Coherence properties of the radiation from X-ray free electron laser
    
    :param cfr: complex wavefield
    :param tstep: temporal step between slices
    
    :returns tau: coherence time [s]
    """
    


    
    mmap = memory_map(map_loc = map_loc,
                      shape = cfr.shape,
                      dtype = 'complex64')
    
    nz0 = cfr.shape[-1]
    
    if bins == 1:
        nz1 = nz0
    else:
        cfr = binArray(cfr, axis = -1, binstep = nz0//bins, binsize = 1 )
        
        nz1 = cfr.shape[-1]
        tStep *= (nz0/nz1)
    
    g = np.zeros([*cfr.shape], dtype = 'complex64')
    
    if VERBOSE: 
        print("Calculating Coherence Time")

    if mpi:
        processes = mp.cpu_count()//2
        pool = mp.Pool(processes)
        pool.map(partial(get_longitudinal_coherence, cfr = cfr, map_loc = map_loc),
                 range(cfr.shape[-1]))
        g = readMap(map_loc, cfr.shape, dtype = 'complex64')        
    else:        
        for i in tqdm(range(cfr.shape[-1])):            
            g[:,:,i] = get_longitudinal_coherence(slice_no = i, cfr = cfr)
            
    tau = (abs(g)**2).sum(axis = -1)[0,0]
    
    if VERBOSE:
        print("\n")
        print(tau)
        print("Time Step: {} fs".format(tStep*1e15))
        print("Coherence Time: {:.2e} fs".format(tau*1e15*tStep))
    
    del mmap
    os.remove(map_loc)

    return tau*tStep

def get_coherence_time_wpg(wfr, mpi = False, VERBOSE = True):
    
    srwl.SetRepresElecField(wfr._srwl_wf, 't')
    time_step = (wfr.params.Mesh.sliceMax - wfr.params.Mesh.sliceMin)/wfr.params.Mesh.nSlices
    return get_coherence_time(wfr.as_complex(), time_step, mpi = mpi)


def get_coherence_len(wfr, dx, dy, VERBOSE = True):
    """
    Calculate coherence length of a complex wavefield of shape
    [nx, ny. nz]
    
    :param wfr: complex wavefield
    :param dx: horizontal pixel size
    :param dy: vertical pixel size
    
    :returns Jd: complex degree of coherence
    :returns clen: coherence length [m]
    """
    
    
    profile, r = get_complex_radial_profile(wfr)

    
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
        print("Radial Coherence Length: {:.2f} um".format(clen*1e6))
    return clen


def get_coherence_len_wpg(wfr, VERBOSE = True):
    srwl.SetRepresElecField(wfr._srwl_wf, 't')
    return get_coherence_len(wfr.as_complex(), wfr.get_spatial_resolution()[0], wfr.get_spatial_resolution()[1])

def get_transverse_doc(wfr, VERBOSE = True):
    """
    get transverse degree of coherence of the wavefront across each of the
    transverse dimensions slices
    """
    
    p, r =  get_complex_radial_profile(wfr)
    nt = wfr.shape[-1]
    J = np.dot(p, p.T.conjugate())/nt
    
    
    tdoc = np.diag(np.dot(J, J)).sum() / np.diag(J).sum()**2
    
    if VERBOSE:
        print("Transverse Degree of Coherence: {:.4f}".format(tdoc.real))
    
    return tdoc

def get_complex_radial_profile(wfr):
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
    
    prof = np.moveaxis(rp, 0, -1)
    return prof, r
  
    
 

 
def run(wfr):
    
    tstep = get_axis(wfr, axis = 't')
    tstep = tstep[1]-tstep[0]
    
    xstep, ystep = wfr.get_spatial_resolution()
    
    wfr = wfr.as_complex()[0,:,:,:]
    
    
    tau = get_coherence_time(wfr, tstep, VERBOSE=True)
    clen = get_coherence_len(wfr, xstep, ystep, VERBOSE=True)
    tdoc = get_transverse_doc(wfr, VERBOSE=True)
    
    return tau, clen, tdoc




if __name__ == "__main__":
    wfr = construct_SA1_pulse(512,512,10,5.0,0.25)
    
    a = get_coherence_time_wpg(wfr)    
    b = get_coherence_time_wpg(wfr, mpi = True)    
    #get_coherence_len_wpg(wfr)

