#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 20:06:28 2020

@author: twguest
"""
 
import numpy as np

from matplotlib import pyplot as plt
import matplotlib as mpl

import os
from moviepy.editor import ImageSequenceClip
from sklearn.preprocessing import minmax_scale as norm


def triple(ii, title = None, xlabel = None, ylabel = None, clabel = None,
           extent = None, cmap = 'hot', vmin = None, vmax = None, savedir = None,
           cticks = None, clabels = None, resolution = 100):
    
    
    
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize = [12, 8], dpi = resolution)
    
    for itr in range(ii.shape[-1]):
        
        
        im = axes[itr].imshow(ii[:,:,itr], cmap = cmap, vmin = vmin, vmax = vmax)
        
        axes[itr].set_xticks([])
        axes[itr].set_yticks([])
    
        axes[itr].set_xlabel(xlabel)
        axes[itr].set_xlabel(ylabel)

    fig.subplots_adjust(right=0.75)
    fig.suptitle(title)    
    
    
    cbar_ax = fig.add_axes([0.78, 0.365, 0.05, 0.277])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(clabel, fontsize = 16)
    
    if cticks != None and clabels != None:
        cbar.set_ticks(cticks)
        cbar.set_ticklabels(clabels)
    
    if savedir is not None:
        fig.savefig(savedir)
    else:
        plt.show()
        

def arr2gif(fname, array, fps=10, scale=1.0):
    """Creates a gif given a stack of images using moviepy
    see: https://gist.github.com/nirum/d4224ad3cd0d71bfef6eba8f3d6ffd59
    
    :params fname: save directory (str)
    :params array: array to be made into gif (along -1 axis)
    :params fps: frames per second (int)
    :params scale: rescaling factor of image
    """

    # ensure that the file has the .gif extension
    fname, _ = os.path.splitext(fname)
    filename = fname + '.gif'

    # normalise 
    for itr in range(array.shape[-1]):
        array[:,:,itr] = norm(array[:,:,itr]*255)    

    # copy into the color dimension if the images are black and white
    if array.ndim == 3:
        array = array[..., np.newaxis] * np.ones(3)
        
    # make the moviepy clip
    clip = ImageSequenceClip(list(array), fps=fps).resize(scale)
    clip.write_gif(filename, fps=fps)
    print("gif saved @: {}".format(fname))