#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 14:39:58 2020

@author: twguest


Let's test the geometric flow theory
"""

import sys

sys.path.append("/opt/spytlab")
sys.path.append("/opt/WPG/")
sys.path.append("/opt/spb_model")



import numpy as np

from model.materials.phaseMask import phaseMask
from model.beamline.structure import propParams
from model.tools import constructPulse
from utils.banded_utils import diagonal_form, solve_banded
from wpg.optical_elements import Drift
from wpg.beamline import Beamline
from OpticalFlow import processOneProjection

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from sklearn.preprocessing import minmax_scale as norm
from matplotlib import pyplot as plt
from model.src.coherent import coherentSource
from wpg import srwlib
from wpg.wpg_uti_wf import getAxis
from scipy.constants import h,c,e
if __name__ == "__main__":
    
    slc = 2
    N = 10
    
    nx, ny = 128, 128
    II = np.zeros((nx,ny,N))
    PH = np.zeros((nx,ny,N))
    PHz = np.zeros((nx,ny,N))
    A = np.zeros((N,N))
    B = np.zeros((N,1))
    
    
    val = nx//4
    
    wfr = coherentSource(nx,ny,4.96,0.25)
    wav = (h*c)/(wfr.params.photonEnergy*e)
    sp = phaseMask(np.random.rand(50,50), [ getAxis(wfr, axis = 'x').max()-
                                            getAxis(wfr, axis = 'x').min(),
                                            getAxis(wfr, axis = 'y').max()-
                                            getAxis(wfr, axis = 'y').min()], wav) ##speckle
    slc = 2
    N = 10
    
    nx, ny = 128, 128
    II = np.zeros((nx,ny,N))
    PH = np.zeros((nx,ny,N))
    PHz = np.zeros((nx,ny,N))
    A = np.zeros((N,N))
    B = np.zeros((N,1))
    
    
    val = nx//4
    
    for i in range(N):
   
        wfr = coherentSource(nx,ny,4.96,0.25)
        pm = np.random.rand(nx,ny)*1e-2
        print(pm[val,val])
        srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')    
        ps = phaseMask(pm, [ getAxis(wfr, axis = 'x').max()-
                                                getAxis(wfr, axis = 'x').min(),
                                                getAxis(wfr, axis = 'y').max()-
                                                getAxis(wfr, axis = 'y').min()], wav) ##speckle
            
        bl = Beamline()
        bl.append(sp, propParams(1,1,1,1))
        bl.propagate(wfr)
        
        
        bl = Beamline()
        bl.append(ps, propParams(1,1,1,1))
        bl.propagate(wfr)
        
        
        PH[:,:,i] = wfr.get_phase()[:,:,0] #% np.pi*2
        ## PH = (PH + np.pi) % (2 * np.pi) - np.pi


        
        bl = Beamline()
        bl.append(Drift(0.10), propParams(1,1,1,1, mode = 'normal'))
        bl.propagate(wfr)
        II[:,:,i] =  wfr.get_intensity()[:,:,0]
        
        plotIntensity(wfr)
    
    II = np.random.rand(*II.shape)*1e-100
    
    print("")
    
    for i in range(N):
        
        a = II[:,:,i]
        
        if i+1 != N:
            b = II[:,:,i+1]
        else:
            b = a
            
        results = processOneProjection(a,b)
        phi = results['phi'].real
       
        print("Phase Diff: {}".format((phi[val,val] + np.pi) % (2 * np.pi) - np.pi))
        #if i+1 != N:
            #print("M Diff: {}".format((PH[val,val,i] - PH[val,val, i+1])/((phi[val,val] + np.pi) % (2 * np.pi) - np.pi)))
# =============================================================================
#         plt.imshow(norm(phi))
#         plt.show()
# =============================================================================
    
    print("\nactual phase difference")
    for i in range(N):
        
        if i+1 != N:
            print(PH[val,val,i] + PH[val,val, i+1])
        else:
            print(0)
    
# =============================================================================
#         
#         A[i,i] = 1
#         B[i,:] = (phi[val,val])
#         
# 
#         
#         if i+1 != N:
#             A[i,i+1] = -1
#             
#     
#     ab = diagonal_form(A)
#     
#     x = solve_banded((1,1), ab, B)
#     
#     
# =============================================================================
        
        
    #print(phi[val,val] - np.matmul(A,x))
    #print(np.matmul(A,x))
    #print(x)
    #print(PH[val,val,:])


    

