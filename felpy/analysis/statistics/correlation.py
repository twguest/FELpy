#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.0"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

import numpy as np

from sklearn.preprocessing import minmax_scale
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
def norm(arr, lim = (0,1)):
    """
    Normalise each slice of a 3D array
    Repackaging of sklearn.preprocessing.minmax_scale for a 3D array
    
    :param arr: array to be normalised
    :param lim: normalisation range
    
    :returns arr: normalised array [np array]
    """
    
    if arr.ndim == 3:
        for itr in range(arr.shape[-1]):
            
            arr[:,:,itr] = minmax_scale(arr[:,:,itr], feature_range = lim)
    else:
        arr = minmax_scale(arr, feature_range = lim)
        
    return arr

def norm_difference(arr1, arr2, plot = True, sdir = None, VERBOSE = False):
    """
    get the normalised difference of a pair of arrays
    if arr1 and arr2 are 3D, the normalised difference will be determined
    slice-wise. ie., norm(arr)
    :param arr1, arr2: arrays for comparison [numpy arrays]
    
    :returns arr3: normalised difference map
    """
    
    if VERBOSE:
        print("Getting Normalised Difference")
    
    arr1, arr2 = norm(arr1), norm(arr2)
    arr3 = arr1 - arr2
    
    if plot:

        if arr3.ndim == 3:   
            
            for slc in range(arr3.shape[-1]):
                fig = plt.figure()
                ax1 = fig.add_subplot()
                fig.suptitle("Normalised Difference Map: Slice {}".format(slc))
                
                img = ax1.imshow(arr3[:,:,slc], cmap = 'jet', vmin = -1, vmax = 1)
                ax1.set_xticks([])
                ax1.set_yticks([])
                fig.colorbar(img)
                plt.show()
                if sdir != None:
                    fig.savefig(sdir + "_slice{}.png".format(slc))
                else:
                    pass

        else:
            fig = plt.figure()
            ax1 = fig.add_subplot()
            fig.suptitle("Normalised Difference Map")
            img = ax1.imshow(arr3, cmap = 'jet', vmin = -1, vmax = 1)
            fig.colorbar(img)
            ax1.set_xticks([])
            ax1.set_yticks([])
             
            if sdir != None:
                fig.savefig(sdir + "_slice{}.png".format(slc))
            else:
                pass
            
    return abs(arr3)  


def euclidian_distance(arr1, arr2, VERBOSE = False):
    """
    calculate the normalised euclidian distance between two arrays 
    
    :param arr1, arr2: arrays for comparison  [numpy arrays]
    
    :returns: normalised euclidian distance (0-1) [float]
    """
        
    if VERBOSE:
        print("Getting Normalised Euclidian Distance")
    
    
    return 0.5*(np.std(arr1-arr2)**2) / (np.std(arr1)**2+np.std(arr2)**2)

 
def correlation_plot(corr, mesh, label = "", title = "", sdir = None, cmap = 'jet'):
    """ 
    plot the correlation function of a pair of 2D arrays (x,y)
    
    :param corr: 2D correlation array (via get_correlation)
    :param mesh: coordinate mesh [np array]
    :param sdir: save directory for output .png
    :param label: figure label
    :param title: figure title
    :param cmap: figure color map
    """
    
    fig, ax1 = plt.subplots()
    
    img = ax1.imshow(corr, cmap = cmap,
                     extent = [np.min(mesh[1])*1e6, np.max(mesh[1])*1e6,
                               np.min(mesh[0])*1e6, np.max(mesh[0])*1e6],
                     vmin = 0, vmax = 0.5)
    
    fig.suptitle = title

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='7.5%', pad=0.05)

    cbar = fig.colorbar(img, cax)
    cbar.set_label("Normalised Difference")
    
    ax1.set_xlabel("x [$\mu$m]")
    ax1.set_ylabel("y [$\mu$m]")
    
    ax1.annotate(label, horizontalalignment = 'left',   
                    verticalalignment = 'bottom',
                    xy = (0,1),
                    c = 'white')
    if sdir is None:
        fig.show()
    else:
        fig.savefig(sdir)
        plt.show()
        
        

if __name__ == '__main__':
    
    mode = input("wpg/exp: ")
    
    if mode == 'exp':
        
        proposal = input("Proposal ID: ")
        exp = input("Exp ID: ")
        run = input("run ID: ")