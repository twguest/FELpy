# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
from felpy.analysis.statistics.correlation import norm_difference
import multiprocessing as mp
from felpy.utils.np_utils import memory_map, readMap
#from felpy.exp.shimadzu.preprocess import shimadzu_test_data
from felpy.utils.os_utils import mkdir_p
import shutil
from functools import partial
#### correlation methods
from felpy.analysis.statistics.correlation import norm
import os
from felpy.utils.job_utils import JobScheduler
from labwork.about import logs, dCache
import sys
from felpy.utils.os_utils import mkdir_p

"""
to run on MAXWELL

in terminal

source activate optics
cd /gpfs/exfel/data/user/guestt/FELpy/felpy/analysis/optics/scalar/
python -c "from intensity_analysis import get_intensity_autocorrelation_ensemble as run;run()"
...
wait
"""
def get_intensity_autocorrelation_ensemble(array_dir = "/gpfs/exfel/data/user/guestt/labwork/dCache/whitefield_data/cropped_intensity_r0046.npy"):

    """
    This part launches the jobs that run in main 
    """
    
    arr = np.load(array_dir)
    
    cwd = os.getcwd()
    script = "intensity_correlations.py"
    
    js = JobScheduler(cwd + "/" + script, logDir = logs,
                      jobName = "autocorrelation_analysis_", partition = 'exfel', nodes = 1, jobType = 'array',
                      jobArray = range(arr.shape[-1]))
        
    js.run(test = False)
    
    
    

def get_intensity_autocorrelation_train(train_no, array, map_loc, sdir = None, mpi = True):
    
    array = array[:,:,:,train_no]
 
    map_loc = map_loc + "autocorrelation_map_{}".format(train_no)
    mmap = memory_map(map_loc, shape = (*array.shape, array.shape[-1]))
    
    pool = mp.Pool(mp.cpu_count()//2)
    
    func = partial(get_intensity_autocorrelation_pulse, array = array,
                   map_loc = map_loc)
    
    pool.map(func, range(array.shape[-1]))
    
    g = readMap(map_loc, (*array.shape, array.shape[-1]))       
    
    del mmap
    
    
    if sdir is not None:
        np.save(sdir + "autocorrelation_{}".format(train_no), g)
    
    del g
    
    print("Completed Autocorrelation Calculations for Train: {}".format(train_no))
    
    
def get_intensity_autocorrelation_pulse(pulse_no, array, map_loc):
    """
    computes the time-lagged intensity-intensity autocorrelation for a single
    slice in time (as denoted by slice_no)
    """
    ### create a memory map for this pulse-train
    
    ### send each pulse to its own processor
    
    mmap = memory_map(map_loc, shape = (*array.shape, array.shape[-1]))
    
    pulse = np.repeat(array[:,:,pulse_no][:, :, np.newaxis], array.shape[-1], axis=-1)

    mmap[:,:, pulse_no, :] = pulse*array 

if __name__ == '__main__':
    array_dir = "/gpfs/exfel/data/user/guestt/labwork/dCache/whitefield_data/cropped_intensity_r0046.npy"
    train_no = sys.argv[1]
    array = np.load(array_dir)
    map_loc = dCache + "/tmp/"
    sdir = dCache + "/whitefield_data/r0046/"
    
    mkdir_p(map_loc)
    mkdir_p(sdir)
    
    get_intensity_autocorrelation_train(train_no, array, map_loc, sdir)
    