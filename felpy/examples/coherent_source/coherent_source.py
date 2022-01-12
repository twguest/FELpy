import numpy as np

from felpy.model.src.coherent import construct_SA1_wavefront, construct_SA1_pulse
from felpy.model.src.coherent import analytical_pulse_divergence,analytical_pulse_duration,analytical_pulse_energy,analytical_pulse_width


from matplotlib import pyplot as plt

from wpg.wpg_uti_wf import calculate_fwhm, look_at_q_space

import wpg.srwlpy
from wpg import srwlib


import seaborn as sns

def test_coherent_pulse_fwhm(E, sdir = None):
    """
    Test that the FWHM of the source matches the analytical prediction
    
    :param E: list/array of energies [keV]
    :param sdir: output directory [str]
    """
    
    print("Testing Source Size vs. Energy and Resolution")
    
    sns.set()

    
    fig = plt.figure(figsize = [12,8])
    ax = fig.add_subplot(111)
    ax.set_title("Source Size Radiation Dependence", fontsize = 22)
    ax.set_ylabel("FWHM [$\mu$m]", fontsize = 22)
    ax.set_xlabel("Radiation Energy [keV]", fontsize = 22)
    ax.set_xlim([2.5,17.5])
    ax.set_ylim([0,55])
    
    E
    analytical_data = np.array([analytical_pulse_width(a)*1e6 for a in E])

    fwhms = []

    for energy in E:
        wfr = construct_SA1_wavefront(nx = 1024, ny = 1024, ekev = energy, q = 0.25)
        fwhms.append(calculate_fwhm(wfr)['fwhm_x']*1e6)

    simulation, = ax.plot(E, fwhms, 'r')
    analytical, = ax.plot(E, analytical_data, '--b')

    leg1 = ax.legend([simulation,analytical], ["Simulated Data", "Analytical Model"], loc = 'lower left', fontsize = 22)
    
    ax.add_artist(leg1)
    
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='both', which='minor', labelsize=14)
        

    
    if sdir is not None:
        fig.savefig(sdir + "coherent_source_test_fwhm.eps")
            
    plt.show()    
    return wfr


def test_coherent_pulse_divergence(E, Q, sdir):
    """
    Test that the FWHM of the source matches the analytical prediction
    
    :param E: list/array of energies [keV]
    :param q: list/array of beam-charges [nC]
    :param outdir: output directory [str]
    """
    
    sns.set()
    
    fig = plt.figure(figsize = [12,8])
    ax = fig.add_subplot(111)
    ax.set_title("Source Divergence FWHM", fontsize = 22)
    ax.set_ylabel(r"d$\theta$ [$\mu$rad]", fontsize = 16)
    ax.set_xlabel("Radiation Energy [keV]", fontsize = 16)
    ax.set_xlim([2.5,17.5])
    
    
    
    nx, ny = 1024,1024
  
    for q in Q:
        
        divergences = []
        
        for energy in E:
            wfr = construct_SA1_pulse(512, 512, 2, energy, q)
            divergences.append(wfr.get_divergence()[0])
            
        ax.scatter(E, np.array(divergences)*1e6, marker = 'o')
        
        analytical_data = [analytical_pulse_divergence(q, energy) for energy in E]
        ax.plot(E, np.array(analytical_data)*1e6, '--')
    
    ax.legend(["{} nC".format(q) for q in Q])
    
    if sdir == None:
        pass
    else:
        fig.savefig(sdir + "coherent_source_test_divergence.eps")




def core(E = np.linspace(3, 16, 15), Q = [0.25], sdir = None):
    
    test_coherent_pulse_fwhm(E, sdir)
    #test_coherent_pulse_divergence(E, Q, sdir)

if __name__ == '__main__':
    core(sdir = "")