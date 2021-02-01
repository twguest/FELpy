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

# =============================================================================
# def get_intensity_autocorrelation1(arr):
#     """
#     get the correlation between two intensity arrays
#     
#     :param arr: intensity array one (numpy array) [float]
#         this array should be of the form [nx, ny, npulses, ntrains]
#     :returns cor: intensity correlation (numpy array) [float]
#     """
#     return np.mean(arr1*arr2, axis = -1)/(np.mean(arr1, axis = -1)*np.mean(arr2, axis = -1))
# 
# =============================================================================



    
    

def intra_train_correlation(train, arr, method = norm_difference):

    tmp = np.zeros([arr.shape[0],
                    arr.shape[1],
                    arr.shape[-2],
                    arr.shape[-2]])

    for p1 in range(arr.shape[-2]):
        for p2 in range(arr.shape[-2]):
            
            tmp[:,:,p1,p2] = norm_difference(arr[:,:,p1,train],
                                             arr[:,:,p2, train],
                                  plot = False)
        
    return tmp.mean(axis = -1).mean(axis = -1)

def inter_train_correlation(train, arr):
    
    tmp = np.zeros([arr.shape[0],
                    arr.shape[1],
                    arr.shape[-1]])
    
    for t2 in range(arr.shape[-1]):
        
        a = norm_difference(arr[:,:,:,train],
                            arr[:,:,:,t2],
                            plot = False)
        tmp[:,:, t2] = a.mean(axis = -1)

    return tmp

def positional_pulse_correlation (position, arr):
 
            
    tmp = np.zeros([arr.shape[0],
                    arr.shape[1],
                    arr.shape[-1],
                    arr.shape[-1]])
    
    
    for t1 in range(arr.shape[-1]):
        
        for t2 in range(arr.shape[-1]):
            

            tmp[:,:, t1, t2] = norm_difference(arr[:,:,position,t1],
                                               arr[:,:,position,t2],
                                               plot = False)
        
    return tmp.mean(axis = -1)          
    
def sequential_pulse_correlation(position, arr):
    

    tmp = np.zeros([arr.shape[0], arr.shape[1],
                   arr.shape[2], arr.shape[-1]])
    
    for t in range(arr.shape[-1]):
        
        for p in range(arr.shape[-2]):
            tmp[:,:,p,t] = norm_difference(arr[:,:, position, t],
                                           arr[:,:, p,t],
                                           plot = False)
            
    tmp = tmp.mean(axis = -1)
    
    return tmp
            
        
        
        
            
            
        
def quick_plot(arr):
    
    if arr.ndim == 3:
        
        for i in range(arr.shape[-1]):
            
            plt.imshow(arr[:,:,i])
            plt.show()

if __name__ == '__main__':
    pass
# =============================================================================
#     from time import time
#     ii = shimadzu_test_data(250,250,20,10)
#      
#     start = time()
#     C = correlation_analysis(ii, mpi = True)
#     corr_a = C.inter_train_correlation( )
#     fin = time()
#     print("MPI: {:.4f} s".format(fin-start))
#     #quick_plot(corr_a)
# 
#     start = time()
#     C = correlation_analysis(ii, mpi = True)
#     corr_b = C.sequential_pulse_correlation()
#     fin = time()
#     print("MPI: {:.4f} s".format(fin-start))
#     quick_plot(abs(corr_b-corr_a))
# =============================================================================
