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




def plot_centroid(centroid, arg = None, clabel = None):
    
    """
    utility for plotting distribution of beam centroids w/ respect to
    i) the train they come from, ii) their position in the train.
    
    If centroid is 2 dimensional, all pulses will be assumed to be of the same
    train and will be colored with respect to their position in the train. 
    
    If centroid is 3 dimensional, pulses will be colored depending on their position
    within their train.
    
    If arg is supplied, items will be colored depending on arg (ie., time, energy,
    freckle dist. on my body...) and trains will be denoted by marker shape
    
    :param centroid: beam centroids [2D or 3D np array]
    :param arg: optional argument [np array of shape (centroid.shape[0], 1, centroid.shape[-1])]
    """
    
    sns.set_style("white")
    
    
    fig, ax1 = plt.subplots()

   
    if centroid.ndim == 2:
       
        cmap = mpl.cm.jet
        
        if arg is None:
            
            norm = mpl.colors.Normalize(vmin = 0, vmax = centroid.shape[0])
            m = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
            cmap = m.get_cmap()
            
            arg = np.linspace(0, centroid.shape[0], centroid.shape[0])
            
            im1 = ax1.scatter(centroid[:, 0]*1e6, centroid[:,2], d[:,1]*1e6,
                        c = arg,
                        s = 16)
            
            ax1.set_title("Centroid Position", fontsize = 16)
            
            ax1.set_xlabel("x [$\mu$m]", fontsize = 14)
            ax1.set_ylabel("y [$\mu$m]", fontsize = 14)
            
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes('right', size='7.5%', pad=0.05)
            
            cbar = fig.colorbar(im1, cax = cax)
            cbar.set_label("Position in Train", fontsize = 14)

        else:
            norm = mpl.colors.Normalize(vmin = 0, vmax = np.max(arg))
            m = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
            cmap = m.get_cmap()
            
            im1 = ax1.scatter(centroid[:, 0], centroid[:,1],
                  c = arg,
                  cmap = cmap,
                  s = 16)    
            
            
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes('right', size='7.5%', pad=0.05)

            cbar = fig.colorbar(im1, cax = cax)
            
            ax1.set_title("Centroid Position", fontsize = 16)
            
            if clabel != None:
                cbar.set_label(clabel, fontsize = 16)
                
        plt.show()
        
        
            
    elif centroid.ndim == 3:
        
        cmap = mpl.cm.jet    
        
        if arg is None:
            
            norm = mpl.colors.Normalize(vmin = 0, vmax = 256)
            m = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
            cmap = m.get_cmap()
            
            t = normalise(np.linspace(0, centroid.shape[-1],centroid.shape[-1]), lim = (0, 1))
            
            for train in range(centroid.shape[-1]):
                
 
                im1 = ax1.scatter(centroid[:, 0, train]*1e6, centroid[:,1, train]*1e6,
                                  c = cmap(np.linspace(t[train],t[train], centroid.shape[0])),
                                  cmap = cmap,
                                  s = 16)
            
            ax1.set_title("Centroid Position", fontsize = 16)
            
            ax1.set_xlabel("x [$\mu$m]", fontsize = 14)
            ax1.set_ylabel("y [$\mu$m]", fontsize = 14)
            
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes('right', size='7.5%', pad=0.05)
            cbar = fig.colorbar(im1, cax = cax)
            cbar.mappable.set_clim([0,centroid.shape[-1]])
            cbar.mappable.set_cmap(cmap)
            cbar.set_label("Train Number", fontsize = 14)
            
            
        else:
            norm = mpl.colors.Normalize(vmin = 0, vmax = 256)
            m = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
            cmap = m.get_cmap()
            
            narg = normalise(arg, (10,500))
            
            t = normalise(np.linspace(0, centroid.shape[-1],centroid.shape[-1]), lim = (0, 1))

            for train in range(centroid.shape[-1]):
                
 
                im1 = ax1.scatter(centroid[:, 0, train]*1e6, centroid[:,1, train]*1e6,
                                  c = cmap(np.linspace(t[train],t[train], centroid.shape[0])),
                                  cmap = cmap,
                                  s = narg[train,:],
                                  alpha = 0.60)
            
            ax1.set_title("Centroid Position", fontsize = 16)
            1
            ax1.set_xlabel("x [$\mu$m]", fontsize = 14)
            ax1.set_ylabel("y [$\mu$m]", fontsize = 14)
            
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes('right', size='7.5%', pad=0.05)
            cbar = fig.colorbar(im1, cax = cax)
            cbar.mappable.set_clim([0,centroid.shape[-1]])
            cbar.mappable.set_cmap(cmap)
            cbar.set_label("Train Number", fontsize = 14)
                
        plt.show()
        


if __name__ == '__main__':
    
    arr = np.random.rand(100,100, 128, 5)*50
    
    cnt = get_com(arr)
 
    #arg = np.random.rand(arr.shape[-1], arr.shape[-2])
    #plot_centroid(cnt, arg,clabel = "Energy")