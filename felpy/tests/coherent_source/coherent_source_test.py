# -*- coding: utf-8 -*-
from tqdm import tqdm

from felpy.model.src.coherent import construct_SA1_pulse

from felpy.model.src.coherent import analytical_pulse_divergence, analytical_pulse_duration, analytical_pulse_energy, analytical_pulse_width
import numpy as np
from scipy.signal import peak_widths
from matplotlib import pyplot as plt
from wpg.wpg_uti_wf import calc_pulse_energy
from felpy.utils.vis_utils import scatter_plot
import seaborn as sns

fwhm2rms = np.sqrt(8*np.log(2)) ### FWHM = sqrt(8ln(2))*sigma

def get_pulse_energy(wfr):
    wfr.set_electric_field_representation('t')
    energy = calc_pulse_energy(wfr)[0]
    wfr.set_electric_field_representation('f')
    return energy 

def get_all_properties(ekev, q, DEBUG = False, VERBOSE = True):
    
    wfr = construct_SA1_pulse(1024, 1024, 6, ekev, q)
    
    
    if DEBUG:
        print("DEBUG")
        print("1) model")
        print("2) analytical result")
        
    fwhm = wfr.get_beam_size()[0]*2
 
    if VERBOSE:
        
        print("fwhm")
        print(wfr.get_beam_size()*2)
        print(analytical_pulse_width(ekev))
    
    
    wfr.set_electric_field_representation('t')
    energy = calc_pulse_energy(wfr)[0]
    wfr.set_electric_field_representation('f')
    
    if VERBOSE:    
        print("energy")
        print(energy)
        print(analytical_pulse_energy(q, ekev))
    
    duration = wfr.get_pulse_duration()
    if VERBOSE:
        print("duration")
        print(duration)
        print(analytical_pulse_duration(q))
    
    divergence = wfr.get_divergence()[0]
    if VERBOSE:
        print("divergence")
        print(divergence)
        print(analytical_pulse_divergence(q, ekev))
        

    return fwhm, energy, duration, divergence
 

def run_divergence_test():
 
 
    EKEV = np.arange(3,15)
    Q =  [0.1, 0.25, 0.5, 1.0]
    q = 0.1
    
    fig, ax = plt.subplots()
    
    for q in Q:
        divs = []
        for ekev in EKEV:
            divs.append(construct_SA1_pulse(2000, 2000, 2, ekev, q).get_divergence()[0]*1e6)
        
        ax = scatter_plot(EKEV, divs, return_axes = True, parse_axes = ax,
                          xlabel = "Energy (keV)",
                          ylabel = "Beam Divergence ($\mu$rad)",
                          marker = 'o',
                          marker_size = 6,
                          legend = True,
                          label = '{}'.format(q),
                          legend_title = "Beam Charge (nC)")
        
        ax.plot(EKEV, analytical_pulse_divergence(q, EKEV)*1e6)

    #ax.set_xlim(2, 15)


def run_width_test(efraction = 0.760968108550488):
 
 
    EKEV = np.arange(3,15)
 
    q = 0.5
    
    #sns.set_style('white')
    fig, ax = plt.subplots()
    

    divs = []
    errs = []
    for ekev in EKEV:
        
        wfr = construct_SA1_pulse(1024, 1024, 2, ekev, q)
        res, err = wfr.get_beam_size(efraction, threshold = 0.0001)
        res = res[0]* 1e6 *2
        
        err = err/(efraction/res)
        err *= wfr.get_spatial_resolution()[0]*1e6*2
        
        divs.append(res)
        errs.append(err)
    
    ax = scatter_plot(EKEV, divs, return_axes = True, parse_axes = ax,
                      xlabel = "Energy (keV)",
                      ylabel = "Beam Width ($\mu$m)",
                      marker = 'o',
                      marker_size = 6)
    
    ax.errorbar(EKEV, divs, yerr = errs, fmt = 'none')
    
    ax.plot(EKEV, analytical_pulse_width(EKEV)*1e6*np.sqrt(2*np.log(2)))

    ax.set_xlim(2, 15)
    

def run_duration_test():
 
 
    Q = np.arange(0.1, 1.0, (1.0-0.1)/20)
 
    ekev = 9.0
    
    #sns.set_style('white')
    fig, ax = plt.subplots()
    

    divs = []
    for q in Q:
        divs.append(construct_SA1_pulse(512, 512, 2, ekev, q).get_pulse_duration()*1e15)
    
    ax = scatter_plot(Q, divs, return_axes = True, parse_axes = ax,
                      xlabel = "Beam Charge (nC)",
                      ylabel = "Pulse Duration (fs)",
                      marker = 'o',
                      marker_size = 6)
    
    ax.plot(Q, analytical_pulse_duration(Q)*1e15)

    
def run_energy_test():
 
 
    EKEV = np.arange(3,15)
    Q =  [0.1, 0.25, 0.5, 1.0]
    q = 0.1
    
    #sns.set_style('white')
    fig, ax = plt.subplots()
     
    for q in Q:
        divs = []
        for ekev in EKEV:
            divs.append(get_pulse_energy(construct_SA1_pulse(512, 512, 5, ekev, q))*1e3*0.865)
        
        ax = scatter_plot(EKEV, divs, return_axes = True, parse_axes = ax,
                          xlabel = "Energy (keV)",
                          ylabel = "Pulse Energy (mJ)",
                          marker = 'o',
                          marker_size = 6,
                          legend = True,
                          label = '{}'.format(q),
                          legend_title = "Beam Charge (nC)")
        
        ax.plot(EKEV, analytical_pulse_energy(q,EKEV)*1e3)

    #ax.set_xlim(2, 15)

def run_all():
# =============================================================================
#     run_divergence_test()
#     
# =============================================================================
    run_width_test()
# =============================================================================
#     run_duration_test()
#     run_energy_test()
# =============================================================================

if __name__ == '__main__':
    DEBUG = False
    
    if DEBUG:
        get_all_properties(11,0.5)
    
    run_all()