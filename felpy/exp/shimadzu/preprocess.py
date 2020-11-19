# -*- coding: utf-8 -*-
"""
This is a script used exclusively to process the shimadzu data.  
"""
import os 
import shutil
#from karabo_data import RunDirectory
#from felpy.utils.daq_utils import load_data, shimadzu_reshape

import numpy as np
from labwork.about import dCache
from felpy.utils.os_utils import mkdir_p
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from felpy.model.tools import create_circular_mask
from felpy.analysis.statistics.univariate import mean_intensity
from felpy.utils.vis_utils import extract_animation, basic_plot, animate
from felpy.utils.np_utils import get_mesh, gaussian_2d
from felpy.analysis.statistics.correlation import norm_difference, correlation_plot
from copy import copy
import seaborn as sns
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import gaussian_kde
from felpy.analysis.optics.scalar.centroid import get_centroid
from felpy.analysis.optics.scalar.enclosed_energy import get_enclosed_energy

from felpy.analysis.statistics.kde import get_kde_kernel, plot_kde, get_kde
from felpy.analysis.statistics.cdf import get_cdf, plot_cdf_data
from felpy.utils.daq_utils import load_data, shimadzu_reshape


def shimadzu_test_data(nx, ny, npulses, ntrains, weight = 0.1,
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
    
    data = np.zeros([nx,ny,npulses,ntrains])
    
    for train in range(ntrains):
        for pulse in range(npulses):
            data[:,:,pulse,train] = gaussian_2d(nx,ny) + np.random.rand(nx, ny)*weight
    
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