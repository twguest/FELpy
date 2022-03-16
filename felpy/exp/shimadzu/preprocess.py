# -*- coding: utf-8 -*-
"""
This is a script used exclusively to process the shimadzu data.  
"""
import os 
import shutil

import numpy as np
from labwork.about import dCache
from felpy.utils.os_utils import mkdir_p

from felpy.utils.np_utils import get_mesh, complex_gaussian_2d
from copy import copy

from felpy.utils.daq_utils import load_data, shimadzu_reshape


def shimadzu_test_data(nx, ny, npulses, ntrains, weight = 0.01,
                       processed = False):
    """
    generate data matching the structure of 'treated' shimadzu test data
    
    :param nx: number of horizontal pixels
    :param ny: number of vertical pixels
    :param npulses: number of pulses per train
    :param ntrains: number of trains
    :param weight: weighting of random modulation (1 leads to unccorelated fields)
     
    :returns data: numpy array of shape [nx, ny, npulses, ntrains]
    """
    
    data = np.zeros([nx,ny,npulses,ntrains]).astype('complex64')
    
    for train in range(ntrains):
        for pulse in range(npulses):
            data[:,:,pulse,train] = complex_gaussian_2d(nx,ny, sigma = nx/10) + np.random.rand(nx, ny)*weight
    
    return data


    
def preprocess_shimadzu_data(proposal, exp, run, px):
    
        sdir = dCache + "/NKB_sensing/whitefield/" + run + "/"
        print("Saving results to dCache: {}".format(sdir))
        mkdir_p(sdir)
        print("Directory exists: {}".format(os.path.exists(sdir)))
        
   
        data = load_data(run, proposal, exp)
        
        ii = data.get_array('SPB_EHD_HPVX2_1/CAM/CAMERA:daqOutput',
                            'data.image.pixels')
        
 
        ii = shimadzu_reshape(ii)

        mesh = get_mesh(ii, px, px)
        print("Mesh Shape: {}".format(mesh.shape))
        
        print("")
        print("Shimadzu Data Pre-Processed")
        return ii, mesh

    


if __name__ == '__main__':
# =============================================================================
#         
#     gifDir = sdir + "gifs/"
#     mkdir_p(gifDir)
# =============================================================================
        
# =============================================================================
#     
#     proposal = 1
#     run = 2 
#     exp = 3
#     px = 4
#     
#     sdir = dCache + 1
#     
# =============================================================================

    pass