# -*- coding: utf-8 -*-

import numpy as np
import seaborn as sns
import matplotlib as mpl
from multiprocessing import pool, cpu_count
from felpy.analysis.statistics.correlation import norm as normalise
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from os import listdir
from scipy.ndimage import center_of_mass

import pandas as pd


def get_com(arr):
    
    """ 
    return the center of mass of a 2D array. If the input array is three-dimensional,
    return the slice wise list of center. If four-dimensional, likewise, iterate
    over the third and fourth dimensions.
    
    assumes dimensions [nx, ny, nz, nt] or likewise
    
    :param arr: input array, ie., beam intensity [np array]
    
    :returns centroid: numpy array containing index of beam centroids.
    """
    
    if arr.ndim == 2:
        
        ### centroid is a tuple
        centroid = center_of_mass(arr)
 
    elif arr.ndim == 3:
        
        ### centroid shape [nz, 2]
        centroid = [center_of_mass(arr[:,:,i]) for i in range(arr.shape[-1])]
    
    elif arr.ndim == 4:
        
        centroid = np.zeros([arr.shape[-2], 2, arr.shape[-1]])
        
        for itr in range(arr.shape[-1]):
        
            c = [center_of_mass(arr[:,:,i,itr]) for i in range(arr.shape[-2])]
            centroid[:,:,itr] = np.asarray(c)
    centroid = np.asarray(centroid)

    return centroid

def com_to_h5(image, outdir, px = 1, py = 1):
    """ 
    wrapper function to write center-of-mass outputs to a h5 file
    
    :param image: image to be analysed
    :param outdir: directory of h5 file
    :param px: horizontal pixel size (m)
    :paray py: vertical pixel size (m)
    """
    
    dict_com_x = {}
    dict_com_y = {}

    com = get_com(image)


    for itr in range(image.shape[-1]):

        dict_com_x['Train {}'.format(itr)] = com[:,0,itr]*px
        dict_com_y['Train {}'.format(itr)] = com[:,1,itr]*py


    df_x = pd.DataFrame.from_dict(dict_com_x)
    df_y = pd.DataFrame.from_dict(dict_com_y)

    df_x.to_hdf(outdir, key = "com_x")
    df_y.to_hdf(outdir, key = "com_y")


if __name__ == '__main__':
    
    arr = np.random.rand(100,100, 128, 5)*50
    
    cnt = get_com(arr)
 
    #arg = np.random.rand(arr.shape[-1], arr.shape[-2])
    #plot_centroid(cnt, arg,clabel = "Energy")