import numpy as np
from felpy.analysis.statistics.correlation import norm_difference
from multiprocessing import Pool, cpu_count
from felpy.model.tools import memoryMap#### correlation methods

def correlation_method_a(arr, mpi = False):
    ### this is equivalent to 'inter' analysis
    
    marr = arr.mean(axis = -1).mean(axis = -1) ## mean array
    
    
    corr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[-1]])
    
    ### get mean correlation between all sequential pulses
    for train in range(arr.shape[-1]):
        
        tmp = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2]])

        for pulse in range(arr.shape[-2]):
        
            tmp[:,:, pulse] = norm_difference(arr[:,:,pulse,train],
                                               marr,
                                               plot = False)
            
        corr[:,:,train] = tmp.mean(axis = -1)
        
        return corr
    


def get_correlation(arr, mode = 'sequential'):
    
    if mode == 'sequential':
        
        tcorr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2]-1, arr.shape[-1]])

        ### get mean correlation between all sequential pulses
        for train in range(arr.shape[-1]):
            

            for pulse in range(tcorr.shape[-2]):
            
                tcorr[:,:, pulse, train] = norm_difference(arr[:,:,pulse,train],
                                                   arr[:,:,pulse+1,train],
                                                   plot = False)         
    
    elif mode == 'inter':
        ### get mean correlation between each train and the mean wavefield
        
        marr = arr.mean(axis = -1).mean(axis = -1) ## mean array
        
        tcorr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[-1]])
        
        ### get mean correlation between all sequential pulses
        for train in range(arr.shape[-1]):
            
            tmp = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2]])

            for pulse in range(arr.shape[-2]):
            
                tmp[:,:, pulse] = norm_difference(arr[:,:,pulse,train],
                                                   marr,
                                                   plot = False)
                
            tcorr[:,:,train] = tmp.mean(axis = -1)
    
    
    elif mode == 'intra':
        ### get the mean correlation between all pulses in a train
        ### looks to be heavy
        
        tcorr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[-1]])
             
        ### get mean correlation between all sequential pulses
        for train in range(arr.shape[-1]):
            
            tmp = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2]])
            marr = arr[:,:,:,train].mean(axis = -1)
            
            for pulse in range(arr.shape[-2]):
    
                tmp[:,:, pulse] = norm_difference(arr[:,:,pulse,train],
                                                   marr,
                                                   plot = False)
                
            tcorr[:,:,train] = tmp.mean(axis = -1)
            
    elif mode == 'position':
        ### get the correlation between pulses position wise

        
        tcorr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[-2]])
             
        ### get mean correlation between all sequential pulses
        for pulse in range(arr.shape[-2]):
            
            tmp = np.zeros([arr.shape[0], arr.shape[1], arr.shape[-1]])
            marr = arr[:,:,pulse,:].mean(axis = -1)
            
            for train in range(arr.shape[-1]):
    
                tmp[:,:, train] = norm_difference(arr[:,:,pulse,train],
                                                  marr,
                                                  plot = False)
                
            tcorr[:,:,train] = tmp.mean(axis = -1)
            
    return tcorr

    
    

def generate_correlation_data(arr, mode, mesh, sdir, fname, delay = 0.06):

    
    if mode == 'sequential':
        
        tcorr = get_correlation(arr, mode = 'sequential')
        
        for train in range(arr.shape[-1]):
               
            tdir = sdir + "/tmp/"
            mkdir_p(tdir)
            
            for pulse in range(tcorr.shape[-2]):

                correlation_plot(tcorr[:,:,pulse,train], mesh,
                                 sdir = tdir + "train_{:04d}_pulse_{:04d}.png".format(train, pulse),
                                 label = "train_{}".format(train))
                
        
            animate(indir = tdir, outdir = sdir,
                    fname = fname + "_sequential_correlation",
                    delay = 0.06,
                    rmdir = True)
            
        np.save(sdir + "/sequential_correlation", tcorr)
            
    if mode == 'inter':
        
        tcorr = get_correlation(arr, mode = 'inter')
        
        for train in range(arr.shape[-1]):
               
            tdir = sdir + "/tmp/"
            mkdir_p(tdir)
            
       
            correlation_plot(tcorr[:,:,train], mesh,
                             sdir = tdir + "train_{:04d}.png".format(train),
                             label = "train_{}".format(train))
            
    
        animate(indir = tdir, outdir = sdir,
                fname = fname + "_interpulse_correlation",
                delay = 0.06,
                rmdir = True)
        
        np.save(sdir + "/interpulse_correlation", tcorr)
        
    if mode == 'intra':
        
        tcorr = get_correlation(arr, mode = 'intra')
        
        for train in range(arr.shape[-1]):
               
            tdir = sdir + "/tmp/"
            mkdir_p(tdir)
 
            correlation_plot(tcorr[:,:,train], mesh,
                             sdir = tdir + "train_{:04d}.png".format(train),
                             label = "train_{}".format(train))
            

        animate(indir = tdir, outdir = sdir,
                fname = fname + "_intrapulse_correlation",
                delay = 0.06,
                    rmdir = True)

        np.save(sdir + "/intrapulse_correlation", tcorr)

    if mode == 'position':
        
        tcorr = get_correlation(arr, mode = 'position')

               
        tdir = sdir + "/tmp/"
        mkdir_p(tdir)
        
        for pulse in range(arr.shape[-2]):

            correlation_plot(tcorr[:,:,pulse], mesh,
                             sdir = tdir + "pulse_{:04d}.png".format(pulse),
                             label = "pulse position {}".format(pulse))
            
        
        animate(indir = tdir, outdir = sdir,
                fname = fname + "_positional_correlation",
                delay = 0.06,
                rmdir = True)
        np.save(sdir + "/positional_correlation", tcorr)