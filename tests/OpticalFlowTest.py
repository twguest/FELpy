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
from scipy.ndimage.filters import gaussian_filter


if __name__ == "__main__":
    
    slc = 2
    N = 2
    
    nx, ny = 512, 512
    II = np.zeros((nx,ny,N), 'float32')
    PH = np.zeros((nx,ny,N), 'float32')
    PHz = np.zeros((nx,ny,N), 'float32')
    A = np.zeros((N,N), 'float32')
    B = np.zeros((N,1), 'float32')
    
    

    
    
    wfr = coherentSource(nx,ny,4.96,0.25)
    cc = wfr.data.arrEhor
    a = wfr.get_intensity()[:,:,0]
    wav = (h*c)/(wfr.params.photonEnergy*e)
    pm = gaussian_filter(np.random.rand(nx,ny)*1, sigma = 20)

    
    sp = phaseMask(gaussian_filter(np.random.rand(50,50)*1000, sigma = 20)
                   , [ getAxis(wfr, axis = 'x').max()-
                                            getAxis(wfr, axis = 'x').min(),
                                            getAxis(wfr, axis = 'y').max()-
                                            getAxis(wfr, axis = 'y').min()], wav) ##speckle


    ps = phaseMask(pm, [getAxis(wfr, axis = 'x').max()-
                        getAxis(wfr, axis = 'x').min(),
                        getAxis(wfr, axis = 'y').max()-
                        getAxis(wfr, axis = 'y').min()], wav) ##speckle

   
    
   

   
    PH[:,:,0] += wfr.get_phase(slice_number = 0) #% np.pi*2
    
    #bl = Beamline()
    #bl.append(sp, propParams(1,1,1,1))
    #bl.append(Drift(.1), propParams(1,1,1,1, mode = 'normal'))
    #bl.propagate(wfr)
    
  
    II[:,:,0] =  wfr.get_intensity(slice_number = 0)
    
    #######################################
    
    
    wfr = coherentSource(nx,ny,4.96,0.25)
    wfr.data.arrEhor = cc
    b = wfr.get_intensity()[:,:,0]
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')    

    bl = Beamline()
    bl.append(ps, propParams(1,1,1,1))
    bl.propagate(wfr)
    PH[:,:,1] += wfr.get_phase(slice_number = 0) #% np.pi*2
    

    #bl = Beamline()
    #bl.append(sp, propParams(1,1,1,1))
    #bl.append(Drift(.1), propParams(1,1,1,1, mode = 'normal'))
    
    #bl.propagate(wfr)

    II[:,:,1] =  wfr.get_intensity(slice_number = 0)  


    II[np.where(II == 0)] = 1e-10 


    for i in range(N):
       
        if i+1 != N:
            
            for val in [10,20,30,40,50]:
                b = II[:,:,i+1]
                results = processOneProjection(II[:,:,0],II[:,:,1])
                phi = (results['phi3'].real)
                phi = (results['phi3'].real+ np.pi) % (2 * np.pi) - np.pi
                print("Design Phase Diff: ",pm[val,val]*2)
                PH = (PH + np.pi) % (2 * np.pi) - np.pi
                pd = phi[val,val]
                pa = PH[val,val,i]-PH[val,val,i+1]
                print("Phase Diff: {}".format((pd)))
                print("Phase Diff Actual: {}".format(pa))
                #print(pd/pa)
        else:
            pass
        
        plt.imshow(PH[:,:,0]-PH[:,:,1])
        plt.imshow(phi)
            
        
        plt.imshow(wfr.get_phase(slice_number = 0))
        plt.show()
        plt.imshow(PH[:,:,0])
        plt.show()

        plt.imshow(PH[:,:,0]-PH[:,:,1])
    #if i+1 != N:
        #print("M Diff: {}".format((PH[val,val,i] - PH[val,val, i+1])/((phi[val,val] + np.pi) % (2 * np.pi) - np.pi)))
# =============================================================================
#         plt.imshow(norm(phi))
#         plt.show()
# =============================================================================



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




