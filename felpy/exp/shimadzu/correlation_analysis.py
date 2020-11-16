import numpy as np
from matplotlib import pyplot as plt
from felpy.analysis.statistics.correlation import norm_difference
from multiprocessing import Pool, cpu_count
from felpy.model.tools import memoryMap
from felpy.exp.shimadzu.preprocess import shimadzu_test_data
from felpy.utils.os_utils import mkdir_p
import shutil
from functools import partial
from felpy.exp.shimadzu.methods import positional_pulse_correlation, inter_train_correlation, intra_train_correlation, sequential_pulse_correlation

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
                print("MPI Enabled: {} Proc".format(self.processes))
 
        self.p = Pool(processes=self.processes)
        


    def intra_train_correlation(self):
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


    def inter_train_correlation(self):
        

        if self.mpi:


            corr = self.p.map(partial(inter_train_correlation, arr = self.arr),
                         range(self.arr.shape[-1]))
            
            corr = np.stack(corr)
 
            
        
        else:
            pass

        return corr
   
    def positional_pulse_correlation(self):
        

    
        if self.mpi:


            corr = self.p.map(partial(positional_pulse_correlation,
                                      arr = self.arr),
                              range(self.arr.shape[-2]))
            
            corr = np.moveaxis(np.stack(corr), 0, -2) 
 
            
        
        else:
            pass

        return corr
    
    
    def sequential_pulse_correlation(self):
             
    
        if self.mpi:


            corr = self.p.map(partial(sequential_pulse_correlation,
                                      arr = self.arr),
                              range(self.arr.shape[-2]))
            
            corr = np.moveaxis(np.stack(corr), 0, -1) 
 
        else:
            pass
        
        return corr
    
    
if __name__ == '__main__':
    
    from time import time
    ii = shimadzu_test_data(250,250,20,10)
     
    start = time()
    C = correlation_analysis(ii, mpi = True)
    corr_a = C.inter_train_correlation( )
    fin = time()
    print("MPI: {:.4f} s".format(fin-start))
    #quick_plot(corr_a)

    start = time()
    C = correlation_analysis(ii, mpi = True)
    corr_b = C.sequential_pulse_correlation()
    fin = time()
    print("MPI: {:.4f} s".format(fin-start))
    quick_plot(abs(corr_b-corr_a))