#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

import os 
import numpy as np
from matplotlib import pyplot as plt
from felpy.utils.np_utils import get_mesh
from felpy.exp.shimadzu.preprocess import preprocess_shimadzu_data, shimadzu_test_data
from felpy.utils.vis_utils import extract_animation
from labwork.about import dCache, logs
from felpy.utils.job_utils import JobScheduler
from felpy.exp.shimadzu.correlation_analysis import correlation_analysis
from felpy.analysis.optics.scalar.enclosed_energy import get_enclosed_energy
import sys
def extract_whitefield_animation(ii, mesh, sdir, run, modes = ['all']):
       
    for m in modes:
        if os.path.exists(sdir + run + "_{}.gif".format(m)):
            pass
        else:
            print("extracting all images for .gif")
            print("saving to directory: {}{}".format(sdir, run + "_{}.gif".format(m)))
            
            extract_animation(ii, mesh, fname = run + "_{}".format(m), sdir = sdir,
                              mode = m)


def full_correlation_analysis(ii, mesh, sdir, run):
    
    c = correlation_analysis(arr = ii, mpi = True, VERBOSE = True)
    
    corr = c.inter_train_correlation()
    np.save(sdir + "inter_train_correlation_{}".format(run), corr)

    corr = c.intra_train_correlation()
    np.save(sdir + "intra_train_correlation_{}".format(run), corr)    
    
    corr = c.positional_pulse_correlation()
    np.save(sdir + "positional_train_correlation_{}".format(run), corr) 
    
    corr = c.sequential_pulse_correlation()
    np.save(sdir + "sequential_train_correlation_{}".format(run), corr) 
  
    

def beam_size_analysis(ii, mesh, sdir, run):
    
    ntrains = ii.shape[-1]
    npulses = ii.shape[-2]
    
    area = np.zeros([npulses, ntrains])
    
    px = py = 28.89e-06
    
    
    for train in range(ntrains):
        for pulse in range(npulses):
            print("Analysis of: Pulse {} // Train {}".format(pulse, train))
            rx, ry = get_enclosed_energy(ii[:,:,pulse,train], px, py)
            area[pulse,train] = np.pi*rx*ry

    
    np.save(sdir + "beam_area_{}".format(run), area)
    
def launch(methods):
    """
    This part launches the jobs that run in main 
    """
    
    cwd = os.getcwd()
    script = "routines.py"
    
    js = JobScheduler(cwd + "/" + script, logDir = logs,
                      jobName = "r0046_", partition = 'exfel', nodes = 4, jobType = 'array',
                      jobArray = methods)
        
    js.run(test = False)
    
    
if __name__ == '__main__':
    method = globals()[sys.argv[1]]
    sdir = dCache + "whitefield_data/"
    
    run = 'r0046'
    ci = np.load(sdir + "cropped_intensity_{}.npy".format(run))
    cmesh = np.load(sdir + "cropped_mesh_{}.npy".format(run))

    ii = np.load(sdir + "intensity_{}.npy".format(run))
    mesh = np.load(sdir + "mesh_{}.npy".format(run))
    method(ii, mesh, sdir, run)