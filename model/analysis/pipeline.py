#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script houses the analysis pipepline for a set of pulses in a directory.
It primarily concerns data extraction and cleaning; the structure of the main()
function of the script is as follows:
    
    * - yields plot
    ** - yields data of merit
    
    - Data handling methods:
        -Saving methods (ie building memmaps and assigning a location to each
                         pulse)

    - Pulse handling methods:
        - Extraction of intensity *
        - Construction of multi-pulse profiles etc. **
        
    ## PULSE STATISTICS
    
    - Basic statistics [individual pulse]:
        - Pulse Fluence *
        - System Efficiency *
        - Pulse Centroid Intensity and Location **Centroid Density Plot (3D - X,Y,INTENSITY)
            see(https://python-graph-gallery.com/371-surface-plot/)
        - Pulse Energy and Number of Photons. *
    
    - Coherence analysis [individual pulse]
        - Coherence Length Measurement *
        - Coherence Time Measurement *
        - Transverse Degree of Coherence *
    
    - Beam size analysis [Individual Pulse]:
        - Enclosed energy of a single pulse *Comparison w/ Legacy Values

    ## INTRA-PULSE STATISTICS:
            
    - Basic statistics [Intra-Pulse]:
        - Slice Centroid Intensity and Location vs. Time **(4D Scatter - X,Y,Centroid Size, Time)

    - Beam size analysis [Intra-Pulse]:
        - Intra-Pulse Enclosed Energy **PLOT BEAM SIZE v. Pulse Time

    ## PULSE TRAIN STATISTICS:
    
    - Coherence analysis [multi pulse]
        - Coherence Length Measurement ** (Coherence Length v. Number of Pulses)
        - Coherence Time Measurement ** (Coherence Time v. Number of Pulses)
        - Transverse Degree of Coherence** (TDOC v. Number of Pulses)
    
    - Beam size analysis [Multi-Pulse]:
        - Enclosed energy of a integrated pulse** (Beam Size v. Number of Pulses)
        
    
@author: twguest
@version: 0.0.0
"""


import sys
import shutil
import os

sys.path.append("../../")

import numpy as np


from model.tools import mkdir_p

from matplotlib import pyplot as plt ## for testing
from model.tools import constructPulse ## for testing

from model.analysis.coherence import run as coherence
from model.analysis.enclosedEnergy import run as beamSize

from wpg.wavefront import Wavefront


def memoryMap(savedir, fname, shape):
    """
    construct a memory map
    """
    memmap = np.memmap(savedir + fname, mode = "w+",
                       shape = shape, dtype = 'float64')
    return memmap

def readMap(mapdir,shape):
    """
    read a map from mapDir
    """
    
    mp = np.memmap(mapdir, dtype = 'float64', mode = 'r+', shape = shape)
    
    return mp

def setup():
    """
    hardcoded
    construct a dictionary containing all the relevant parameters and directories
    """
    
    params = {}
    
    params['indir'] = "../../data/tmpTest/" ## input directory
    params['global'] = "../../data/procTemp/" ## global data directory
    
    params['bsi'] = "bsi" ## integrated beam size memmap filename
    params['bss'] = "bss" ## pulsed beam size memmap filename
    params['coh'] = 'coh'
    
    nSlice = 6#2250 ### number of slices per pulse
    
    
    ### these should be moved pre-launch 
    
    mkdir_p(params['global'])
    
    bsi = memoryMap(params['global'], params['bsi'],
              shape = len(os.listdir(params['indir']))) ## integrated beam size
    
 
    bss = memoryMap(params['global'], params['bss'], 
              shape=(nSlice, 3, len(os.listdir(params['indir']))))
    
    coh = memoryMap(params['global'], params['coh'], 
              shape=(len(os.listdir(params['indir'])), 3))
    
    
    params["bsi_shape"] = bsi.shape
    params["bss_shape"] = bss.shape
    params["coh_shape"] = coh.shape
    


    print("\nInput Directory: {}".format(params['indir']))
    print("Global Data Directory: {}".format(params['global']))
    print("\nMemmap Directories:")
    print("Integrated Beam Size: {}".format(params['bsi']))
    print("Pulsed Beam Size: {}".format(params['bss']))
    print("Coherence Measurements: {}".format(params['coh']))    

    return params

def launch(multi):
    """ 
    job scheduler function
    """
    pass


def generateTestPulses(savedir, N = 5):
    """
    generate a set of test pulses
    
    :param savedir: directory for test pulses to be saved
    :param N: number of pulses to generate
    """
    print("Constructing {} Test Pulses".format(N))
    
    for n in range(N):
        
        wfr = constructPulse(1024, 1024, nz = 6, tau = 1e-12)
        
        wfr.data.arrEhor*= np.random.uniform(0.75, 1, size = wfr.data.arrEhor.shape)
        
        wfr.store_hdf5(savedir + "testWavefront_{}.h5".format(n))
        print("Storing Wavefront @" + savedir + "testWavefront_{}.h5".format(n))
    

def testBeamSizeIntegrated():
    
    print("making temporary directory")
    tmpdir = "../../data/tmpTest/"
    procdir = "../../data/procTemp/"
    
    mkdir_p(tmpdir)
    mkdir_p(procdir)
   
    generateTestPulses(tmpdir, N = 5)
    
    ### make memmaps, should be moved to launch
    
    beamRad = memoryMap(procdir, "integratedBeamSize", shape = len(os.listdir(tmpdir)))
    
    
    
    for itr in range(len(os.listdir(tmpdir))):
        wfr = Wavefront()
        print("opening "+ tmpdir + os.listdir(tmpdir)[itr])
        wfr.load_hdf5(tmpdir + os.listdir(tmpdir)[itr])
        
        bs = beamSize(wfr, mode = 'integrated')
        beamRad[itr] = bs

        
    np.save(procdir + "integratedBeamSize", beamRad)
    
    print("removing temporary directory")
    shutil.rmtree("../../data/tmpTest/")


def BeamSize(wfr, mode, memMap, ID, VERBOSE = False):
        

    if mode == 'integrated':
        memMap[ID] = beamSize(wfr, mode = mode, VERBOSE = VERBOSE)
    
    if mode == 'pulse':
        memMap[:,:,ID] = beamSize(wfr, mode = mode, VERBOSE = VERBOSE)
       
    
def Coherence(wfr, memMap, ID, VERBOSE):
    
    memMap[ID,:] = coherence(wfr)
    

def testUsage():
    pass


    

if __name__ == '__main__':
    
    params = setup()


    bsi = readMap(params['global'] + params['bsi'], shape = params['bss_shape'])
    bss = readMap(params['global'] + params['bss'], shape = params['bss_shape'])
    coh = readMap(params['global'] + params['coh'], shape = params['coh_shape'])
    
    for i in range(len(os.listdir(params['indir']))):
        
        
        
        ### SETUP PULSES FOR PIPELINE TEST
        fname = os.listdir(params['indir'])[i]
        ID = i
        wfr = Wavefront()
        wfr.load_hdf5(params['indir']+fname)
        
        
        
        ### integrated beam profile
        #print("Calculating Integrated Beam Profiles")
        #BeamSize(wfr, mode = "integrated", memMap = bsi, ID = ID)
        
        ### pulsed beam profile
        #print("calculating pulsed beam profiles")
        #BeamSize(wfr, mode = "pulse", memMap = bss, ID = ID, VERBOSE = False)
        
        ### single-pulse coherence analysis
        Coherence(wfr, coh, ID, VERBOSE = True)
    
    del bsi, bss
    
    bsi = readMap(params['global'] + params['bsi'], shape = params['bss_shape'])
    bss = readMap(params['global'] + params['bss'], shape = params['bss_shape'])
    coh = readMap(params['global'] + params['coh'], shape = params['coh_shape'])
    
