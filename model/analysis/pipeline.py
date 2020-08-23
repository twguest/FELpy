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

sys.path.append("../../")

import numpy as np

from matplotlib import pyplot as plt ## for testing
from model.tools import constructPulse ## for testing

from model.analysis.coherence import run as coherence


def launch(multi):
    """ 
    job scheduler function
    """
    pass

def run(multi = False):
    
    # do some stuff: ie setup memmap etc,
    # then launch() -- inter and intra-pulse
    # then launch multipulse simulations once all jobs have run
    # schedule multipulse stuff for after jobs have run
    
def generateTestPulses(savedir, N = 5):
    """
    generate a set of test pulses
    
    :param savedir: directory for test pulses to be saved
    :param N: number of pulses to generate
    """
    
    for n in range(N):
        
        wfr = constructPulse(500,500,5)
        
        wfr.data.arrEhor*= np.random.uniform(0.75, 1,
                                             size = wfr.data.arrEhor.shape)
        
        wfr.store_hdf5(savedir + "n.h5")
        
def testUsage():
    pass

if __name__ == '__main__':
    testUsage()
    
    





