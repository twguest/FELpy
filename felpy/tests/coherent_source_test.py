# -*- coding: utf-8 -*-
from tqdm import tqdm

from felpy.model.src.coherent import construct_SA1_pulse

from felpy.model.src.coherent import analytical_pulse_divergence, analytical_pulse_duration, analytical_pulse_energy, analytical_pulse_width
import numpy as np

 
from wpg.wpg_uti_wf import calc_pulse_energy

def get_properties(ekev, q, DEBUG = False):
    
    wfr = construct_SA1_pulse(1024, 1024, 4, ekev, q)
    
    
    if DEBUG:
        print("DEBUG")
        print("1) model")
        print("2) analytical result")
        
    fwhm = wfr.get_fwhm()[0]
    print("fwhm")
    print(fwhm)
    print(analytical_pulse_width(ekev))
    
    
    wfr.set_electric_field_representation('t')
    energy = calc_pulse_energy(wfr)[0]
    wfr.set_electric_field_representation('f')
    print("energy")
    print(energy)
    print(analytical_pulse_energy(q, ekev))
    
    duration = wfr.get_pulse_duration()
    print("duration")
    print(duration)
    print(analytical_pulse_duration(q))
    
    divergence = wfr.get_divergence()[0]
    print("divergence")
    print(divergence)
    print(analytical_pulse_divergence(q, ekev))
    

    return fwhm, energy, duration, divergence
 



if __name__ == '__main__':
    
     get_properties(11,0.5)
# =============================================================================
#     ekev = np.arange(5,10)
#     q = np.array([0.1,0.25,0.5,1.0])
#     
#     
#     dataset = np.zeros([ekev.shape[0], q.shape[0], 4]) #energy x charge x params
#     
#     for i in tqdm(range(0, ekev.shape[0])):
#         for k in range(q.shape[0]):
#             #print(get_properties(ekev[i], q[k]))
#             dataset[i, k, :] = get_properties(ekev[i], q[k])
#            
#     from matplotlib import pyplot as plt     
#     fig, ax = plt.subplots()
#     
#     #for i in range(dataset.shape[1]):
#     ax.scatter(ekev,dataset[:,0,0]*1e6)
#     ax.plot(ekev, analytical_pulse_width(ekev)*1e6)
#     
#     fig, ax = plt.subplots()
# 
#     for i in range(dataset.shape[1]):
#         ax.scatter(q[i],dataset[0,i,-1]*1e6)
#         ax.plot(q, analytical_pulse_divergence(q, ekev[0])*1e6)
# =============================================================================
