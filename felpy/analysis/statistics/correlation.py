#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 13:45:57 2020

@author: twguest
"""
import numpy as np

from sklearn.preprocessing import minmax_scale
from matplotlib import pyplot as plt

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
            
    return arr3  


def euclidian_distance(arr1, arr2, VERBOSE = False):
    """
    calculate the normalised euclidian distance between two arrays 
    
    :param arr1, arr2: arrays for comparison  [numpy arrays]
    
    :returns: normalised euclidian distance (0-1) [float]
    """
        
    if VERBOSE:
        print("Getting Normalised Euclidian Distance")
    
    
    return 0.5*(np.std(arr1-arr2)**2) / (np.std(arr1)**2+np.std(arr2)**2)

if __name__ == '__main__':
    
    a = np.random.rand(100,100,4)
    b = np.random.rand(100,100,4)

    c = norm_difference(a, b, plot = True)
    ed = euclidian_distance(a, b)
    print(ed)