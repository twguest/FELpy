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
from felpy.analysis.statistics.correlation import norm
"""
to run on MAXWELL

in terminal

source activate optics
cd /gpfs/exfel/data/user/guestt/FELpy/felpy/analysis/optics/scalar/
python -c "from intensity_correlations import get_intensity_autocorrelation_ensemble as run;run()"
...
wait
...
once run

python -c "from intensity_correlations import compile_data; compile_data()"
...
wait
...
"""

def get_normalised_difference(arr1, arr2):
    """
    |arr1-arr2|
    """
    return 1 - np.abs(norm(arr1)-norm(arr2))
    

def get_second_order_doc(arr1, arr2):
    """
    does not include ensemble averaging 
     
    [I(r_1, t_1)-I^{bar}(r_1)][I(r_1, t_2)-I^{bar}(r_1)]
    """
    mean = np.repeat(arr2.mean(axis = -1)[:, :, np.newaxis], arr2.shape[-1],
                     axis=-1)
    
    return (arr1 - mean)*(arr2 - mean)

def compile_data(indir = dCache + "whitefield_data/r0047/",
                 array_dir = dCache + "whitefield_data/cropped_intensity_r0047.npy",
                methods = [get_normalised_difference, get_second_order_doc]):
    
    nTrains = np.load(array_dir).shape[-1]
    
    for m in methods:
        
        for itr in range(nTrains):
            if itr == 0:
                corr = np.load(indir + "{}/autocorrelation_{}.npy".format(m.__name__, itr))
            else:
                corr += np.load(indir + "{}/autocorrelation_{}.npy".format(m.__name__, itr))
        
        corr /= float(nTrains)
        
        np.save(indir + m.__name__, corr)
        
    
    

def get_intensity_autocorrelation_ensemble(array_dir = dCache + "whitefield_data/cropped_intensity_r0047.npy",
                                           method = get_normalised_difference):

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
    
    
    

def get_intensity_autocorrelation_train(train_no, array, map_loc, method = get_normalised_difference, sdir = None, mpi = True):
    
    array = array[:,:,3:103,train_no]
 
    map_loc = map_loc + "autocorrelation_map_{}".format(train_no)
    mmap = memory_map(map_loc, shape = (*array.shape, array.shape[-1]))
    
    pool = mp.Pool(mp.cpu_count()//2)
    
    func = partial(get_intensity_autocorrelation_pulse, array = array,
                   map_loc = map_loc, method = method)
    
    pool.map(func, range(array.shape[-1]))
    
    g = readMap(map_loc, (*array.shape, array.shape[-1]))       
    
    del mmap
    
    
    if sdir is not None:
        np.save(sdir + "autocorrelation_{}".format(train_no), g)
    
    del g
    
    print("Completed Autocorrelation Calculations for Train: {}".format(train_no))
    
    

def get_intensity_autocorrelation_pulse(pulse_no, array, map_loc, method = get_normalised_difference):
    """
    computes the time-lagged intensity-intensity autocorrelation for a single
    slice in time (as denoted by slice_no)
    """
    ### create a memory map for this pulse-train
    
    ### send each pulse to its own processor
    
    mmap = memory_map(map_loc, shape = (*array.shape, array.shape[-1]))
    
    pulse = np.repeat(array[:,:,pulse_no][:, :, np.newaxis], array.shape[-1], axis=-1)
    
    
    mmap[:,:, pulse_no, :] = method(pulse, array)


def softmax_normalisation(arr):
    """
    softmax normalisation, ie., let sum of each 2D slice in an array = 1
    """
    
    if arr.ndim == 3:
        
        for itr in range(arr.shape[-1]):
            
            arr[:,:,itr] /= arr[:,:,itr].sum()
    
    return arr

if __name__ == '__main__':
    
    array_dir = dCache + "/whitefield_data/cropped_intensity_r0047.npy"
    

    train_no = int(sys.argv[1])
    
   
    method = [get_normalised_difference, get_second_order_doc] 
    
    
    if type(method) == list:
        
        for m in method:
        
            array = np.load(array_dir).astype('float64')[:,:,3:103:2,:]
            array = np.delete(array, 24, axis = -2)
            map_loc = dCache + "/tmp/"
            sdir = dCache + "/whitefield_data/r0047/"
            mkdir_p(sdir)
            
            
            sdir += "{}/".format(m.__name__)
            
            mkdir_p(map_loc)
            mkdir_p(sdir)
 
            get_intensity_autocorrelation_train(train_no, array, map_loc = map_loc, sdir = sdir,
                                                method = m)