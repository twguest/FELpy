#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 20:06:28 2020

@author: twguest
"""
import shutil 
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
    
def animate(indir, outdir, fname, delay = 0.1, rmdir = False):
    """
    create gif from pngs in directory
    """        
    os.system("convert -delay {} {}/*.png {}{}.gif".format(delay, indir, outdir, fname) )
    
    if rmdir == True:
        shutil.rmtree(indir)
        

def extract_animation(ii, mesh, fname, sdir, mode = 'train'):
    
    tdir = sdir + "/tmp/"
    mkdir_p(tdir)
    
    if mode == 'all':
        for train in range(ii.shape[-1]):
            for pulse in range(ii.shape[-2]):
                
                basic_plot(ii[:,:,pulse,train], mesh,
                           sdir = tdir + "train_{}_pulse_{}.png".format(train, pulse),
                           label = "train/pulse: {}/{}".format(train+1,pulse+1))
    if mode == 'train':
        
        ii = mean_intensity(ii, mode = 'train')
        
        for train in range(ii.shape[-1]):
            
            basic_plot(ii[:,:,train], mesh,
                       sdir = tdir + "train_{}.png".format(train),
                       label = 'train: {}'.format(train+1))
        
    animate(indir = tdir, outdir = sdir, fname = fname, delay = 0.06)


def basic_plot(ii, mesh, sdir = None,
               label = None,
               title = None,
               cmap = 'bone'):
    """ 
    a simple plot of some two-dimensional intensity array
    
    :param ii: 2D intensity array
    :param mesh: coordinate mesh [np array]
    """
    
    fig, ax1 = plt.subplots()
    
    img = ax1.imshow(ii, cmap = cmap,
                     extent = [np.min(mesh[1])*1e6, np.max(mesh[1])*1e6,
                               np.min(mesh[0])*1e6, np.max(mesh[0])*1e6])
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='7.5%', pad=0.05)
    
    cbar = fig.colorbar(img, cax)
    cbar.set_label("Intensity (a.u.)")
    
    ax1.set_xlabel("x [$\mu$m]")
    ax1.set_ylabel("y [$\mu$m]")
    
    if title != None:
        ax1.set_title = title
    if label != None:
        ax1.annotate(label, horizontalalignment = 'left',
                     verticalalignment = 'bottom',
                     xy = (0,1),
                     c = 'white')
    if sdir is None:
        plt.show()
    else:
        fig.savefig(sdir)
        plt.show()
    