# -*- coding: utf-8 -*-
"""
some routines that run successively for shimadzu whitefield data analysis
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

def burger_w_the_lot():
    pass


def extract_whitefield_animation(ii, mesh, sdir, run, modes = ['all']):
       
    for m in modes:
        if os.path.exists(sdir + run + "_{}.gif".format(m)):
            pass
        else:
            print("extracting all images for .gif")
            print("saving to directory: {}{}".format(sdir, run + "_{}.gif".format(m)))
            
            extract_animation(ii, mesh, fname = run + "_{}".format(m), sdir = sdir,
                              mode = m)


def full_correlation_analysis(ii, sdir, run):
    
    c = correlation_analysis(arr = ii, mpi = True, VERBOSE = True)
    
    corr = c.inter_train_correlation()
    np.save(sdir + "inter_train_correlation_{}".format(run), corr)

    corr = c.intra_train_correlation()
    np.save(sdir + "intra_train_correlation_{}".format(run), corr)    
    
    corr = c.positional_pulse_correlation()
    np.save(sdir + "positional_train_correlation_{}".format(run), corr) 
    
    corr = c.sequential_pulse_correlation()
    np.save(sdir + "sequential_train_correlation_{}".format(run), corr) 
  
    
def launch():
    """
    This part launches the jobs that run in main 
    """
    
    cwd = os.getcwd()
    script = "routines.py"
    
    js = JobScheduler(cwd + "/" + script, logDir = logs,
                      jobName = "r0046_correlation_analysis", partition = 'exfel', nodes = 4, jobType = 'single')
        
    js.run(test = False)
    
if __name__ == '__main__':
    sdir = dCache + "whitefield_data/"
    run = 'r0046'
    ii = np.load(sdir + "cropped_intensity_{}.npy".format(run))
    full_correlation_analysis(ii, sdir, run)