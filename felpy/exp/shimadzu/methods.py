import numpy as np
from matplotlib import pyplot as plt
from felpy.analysis.statistics.correlation import norm_difference
from multiprocessing import Pool, cpu_count
from felpy.model.tools import memoryMap
from felpy.exp.shimadzu.preprocess import shimadzu_test_data
from felpy.utils.os_utils import mkdir_p
import shutil
from functools import partial
#### correlation methods

def intra_train_correlation(train, arr, mode = 'lite'):
    
    if mode == 'lite':
        
        marr = arr.mean(axis = -2) ## mean array
        
        tmp = np.zeros([arr.shape[0],
                        arr.shape[1],
                        arr.shape[-1]])
    
        for pulse in range(arr.shape[-1]):
            
            tmp[:,:,pulse] = norm_difference(arr[:,:,pulse,train],
                                  marr[:,:,train],
                                  plot = False)
            
            return tmp.mean(axis = -1)

def inter_train_correlation(train, arr, mode = 'lite'):
    
    if mode == 'lite':
        
        tmp = np.zeros([arr.shape[0],
                    arr.shape[1],
                    arr.shape[-1]])
       
        marr = arr.mean(axis = -1).mean(axis = -1) ## mean array

        tmp = norm_difference(arr[:,:,:,train].mean(axis = -1), marr,
                                         plot = False)
        
    return tmp

def pulse_position_correlation(position, arr, mode = 'lite'):
    
    if mode == 'lite':
            
        tmp = np.zeros([arr.shape[0],
                        arr.shape[1],
                        arr.shape[-1]])
        
        marr = arr.mean(axis = -1) ## mean array

        for train in range(arr.shape[-1]):
            tmp[:,:, train] = norm_difference(arr[:,:,position,train],
                                              marr[:,:,position],
                                              plot = False)
            
        return tmp.mean(axis = -1)          
    
def sequential_pulse_correlation(position, arr, mode = 'lite'):
    
    if mode == 'lite':
        
        tmp = np.zeros([arr.shape[0], arr.shape[1],
                        arr.shape[2]])

        arr = arr.mean(axis = -1)
        
        
        for p in range(arr.shape[-1]):

            tmp[:,:,p] = norm_difference(arr[:,:,position],
                                         arr[:,:,p],
                                         plot = False)
            
    if mode == 'full':
        
        tmp = np.zeros([arr.shape[0], arr.shape[1],
                       arr.shape[2], arr.shape[-1]])
        
        for t in range(arr.shape[-1]):
            
            for p in range(arr.shape[-2]):
                tmp[:,:,p,t] = norm_difference(arr[:,:, position, t],
                                               arr[:,:, p,t],
                                               plot = False)
                
        tmp = tmp.mean(axis = -1)
        
    return tmp
            
        
        
        
            
            
            


class correlation_analysis():
    
    def __init__(self, arr, mpi = False, VERBOSE = True):
        
        if VERBOSE:
            print("Only MPI Supported.")
            
        self.arr = arr
        self.mpi = mpi
        
        self.output_shape = None
        self.tmp = None
        
        self.processes = cpu_count()//2
        
        if VERBOSE:    
            if self.mpi: 
                print("MPI Enabled: {} Processes".format(self.processes))
 
        self.p = Pool(processes=self.processes)
        


    def intra_train_correlation(self, mode = 'lite'):
        """
        intra-train correlation. 
        
        approaches the question: how well correlated is each pulse in a train
        w/ the average pulse intensity for that train.
        """
                

        
        
        if self.mpi:


            corr = self.p.map(partial(intra_train_correlation, arr = self.arr),
                         range(self.arr.shape[-1]))
 
            corr = np.dstack(corr)
 
        else:
            pass
               
        return corr


    def inter_train_correlation(self, mode = 'lite'):
        

        if self.mpi:


            corr = self.p.map(partial(inter_train_correlation, arr = self.arr),
                         range(self.arr.shape[-1]))
            
            corr = np.dstack(corr)
 
            
        
        else:
            pass

        return corr
   
    def pulse_position_correlation(self, mode = 'lite'):
        

    
        if self.mpi:


            corr = self.p.map(partial(pulse_position_correlation,
                                      arr = self.arr),
                              range(self.arr.shape[-2]))
            
            corr = np.dstack(corr)
 
            
        
        else:
            pass

        return corr
    
    
    def sequential_pulse_correlation(self, mode = 'lite'):
             
    
        if self.mpi:


            corr = self.p.map(partial(sequential_pulse_correlation,
                                      arr = self.arr,
                                      mode = mode),
                              range(self.arr.shape[-2]))
            
            corr = np.moveaxis(np.stack(corr), 0, -1) 
 
            
        
        else:
            pass
        
        return corr
        
def quick_plot(arr):
    
    if arr.ndim == 3:
        
        for i in range(arr.shape[-1]):
            
            plt.imshow(arr[:,:,i])
            plt.show()

if __name__ == '__main__':
    
    from time import time
    ii = shimadzu_test_data(250,250,20,10)
     
    start = time()
    C = correlation_analysis(ii, mpi = True)
    corr_a = C.sequential_pulse_correlation(mode = 'lite')
    fin = time()
    print("MPI: {:.4f} s".format(fin-start))
    #quick_plot(corr_a)

    start = time()
    C = correlation_analysis(ii, mpi = True)
    corr_b = C.sequential_pulse_correlation(mode = 'full')
    fin = time()
    print("MPI: {:.4f} s".format(fin-start))
    quick_plot(abs(corr_b-corr_a))