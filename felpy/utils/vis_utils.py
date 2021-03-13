#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

import shutil 
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from moviepy.editor import ImageSequenceClip
from sklearn.preprocessing import minmax_scale as norm
from felpy.utils.os_utils import mkdir_p
from felpy.analysis.statistics.univariate import mean_intensity
from mpl_toolkits.mplot3d import Axes3D
from felpy.utils.np_utils import extent_from_mesh

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
               xlims = None,
               ylims = None,
               crop = None,
               label = None,
               title = None,
               cmap = 'bone',
               scale = 1e6):
    
    """ 
    a simple plot of some two-dimensional intensity array
    
    :param ii: 2D intensity array
    :param mesh: coordinate mesh [np array]
    """
    
    if scale == 1e6:
        mode = '[$\mu$m]'
        
    elif scale == 1e3:
        mode = '[mm]'
        
    elif scale == 1:
        mode = '[m]'
        
    else:
        mode = "[" + str(scale) + "m]"
        
    fig, ax1 = plt.subplots()
    
    if xlims is not None and ylims is not None:
        roi = extent_from_mesh(mesh, xlims, ylims)
        ii = ii[roi[0]:roi[1], roi[2]:roi[3]]
            
        img = ax1.imshow(ii, cmap = cmap,
                         extent = [mesh[0,0,roi[0]]*scale, mesh[0,0,roi[1]]*scale,
                                   mesh[1,roi[2],0]*scale, mesh[1,roi[3],0]*scale])
    else:
        img = ax1.imshow(ii, cmap = cmap,
                         extent = [np.min(mesh[1])*scale, np.max(mesh[1])*scale,
                                   np.min(mesh[0])*scale, np.max(mesh[0])*scale])
  
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='7.5%', pad=0.05)
    
    cbar = fig.colorbar(img, cax)
    cbar.set_label("Intensity (a.u.)")
    
    ax1.set_xlabel("x {}".format(mode))
    ax1.set_ylabel("y {}".format(mode))
    
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

    

def scatter_3D(x, y, z,
               xlabel = "", ylabel = "", zlabel = "",
               azim = 60, elev = 10, get_2D = False,
               sdir = None):
  
    fig = plt.figure()

    ax1 = fig.add_subplot(projection = '3d')
    ax1.view_init(elev=elev, azim=azim)
    
    if x.ndim == 2 and y.ndim == 2:
        for itr in range(x.shape[-1]):
            ax1.scatter(x[:,itr], z, y[:,itr])
    else:
        for pos in range(len(x)):
            ax1.scatter(x[pos],z[pos],y[pos])
        
    ax1.set_xlabel(xlabel)
    ax1.set_zlabel(ylabel)
    ax1.set_ylabel(zlabel)
    
    if get_2D is True:
        fig2, ax2 = plt.subplots()
        fig3, ax3 = plt.subplots()
        if x.ndim == 2: 
            ax2.scatter(x.mean(axis = -1),z)
            ax2.set_xlabel(xlabel)
            ax2.set_ylabel(zlabel)
        
            
            
            ax3.scatter(y.mean(axis = -1),z)
            ax3.set_xlabel(ylabel)
            ax3.set_ylabel(zlabel)
            
        else:
            ax2.scatter(x,z)
            ax2.set_xlabel(xlabel)
            ax2.set_ylabel(zlabel)
            plt.show()
            
            
            ax3.scatter(y,z)
            ax3.set_xlabel(ylabel)
            ax3.set_ylabel(zlabel)
            plt.show()
    
    if sdir is not None:
        fig.savefig(sdir)
        
        if get_2D is True:
            fig2.savefig(sdir + "2D_xz.png")
            fig3.savefig(sdir + "2D_yz.png")
    
    plt.show()
    
def generate_3D_animation(x, y, z,
                          sdir,
                          title = "",
                          npts = 360,
                          xlabel = "",
                          ylabel = "",
                          zlabel = "",
                          mode = 'default', 
                          get_2D = True):
    
    tdir = sdir + "/tmp/"
    mkdir_p(tdir)
    
    if mode == 'default':
        alim, elim = 360, 360
    
    for ang in range(0, 360, int(360/npts)):
        scatter_3D(x, y, z, 
                   sdir = tdir + "{:5d}.png".format(ang),
                   azim = ang % alim,
                   elev = ang % alim,
                   xlabel = xlabel,
                   ylabel = ylabel, zlabel = zlabel)
    
    animate(tdir, sdir, title + "_3D_animation",
            rmdir = True)


def plot_fill_between(data, title = "",
                      xlabel = "",
                      ylabel = "",
                      alpha = 0.25,
                      context = 'notebook',
                      plot_max = True,
                      xlim = None,
                      ylim = None,
                      color = 'blue'):
    """
    plot the mean of a 2D data set with the 95% confidence interval filled
    """
    #sns.set()
    sns.set_style("dark")
    sns.set_context(context)
    
    fig, ax1 = plt.subplots()
    
    ax1.plot(data.mean(-1), color = color)
    ax1.fill_between(np.arange(data.shape[0]),
                     data.mean(-1) - np.std(data),
                     data.mean(-1) + np.std(data),
                     alpha = alpha, color = color,
                     edgecolor= None, facecolor = None,
                     linewidth=0.0)
    
    if xlim is not None:
        ax1.set_xlim(xlim)
    if ylim is not None:
        ax1.set_ylim(ylim)

    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    

    if plot_max:
        ax1.plot(np.max(data, -1), c = 'blue', linestyle = 'dashdot')
        ax1.plot(np.min(data, -1), c = 'blue', linestyle = 'dotted')
        
    plt.show()
    

def simple_line_plot(data, title = "",
                      xlabel = "",
                      ylabel = "",
                      context = 'notebook',
                      xlim = None,
                      ylim = None,
                      color = 'blue'):
    """
    plot the mean of a 2D data set with the 95% confidence interval filled
    """
    #sns.set()
    sns.set_style("dark")
    sns.set_context(context)
    
    fig, ax1 = plt.subplots()
    
    ax1.plot(data, color = color)

    if xlim is not None:
        ax1.set_xlim(xlim)
    if ylim is not None:
        ax1.set_ylim(ylim)

    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    
    plt.show()
    
    
def colorbar_plot(arr,
                  mesh = None,
                  label = None,
                  title = None,
                  xlabel = None,
                  ylabel = None,
                  clabel = "",
                  context = 'notebook',
                  sdir = None,
                  cmap = 'bone',
                  normalise = False,
                  vmin = 0, vmax = 1,
                  scale = 1e6,
                  aspect = 'auto'):
    
    """ 
    plot a 2D array with a colorbar (x,y)
    
    :param corr: 2D correlation array (via get_correlation)
    :param mesh: coordinate mesh [np array]
    :param sdir: save directory for output .png
    :param label: figure label
    :param title: figure title
    :param cmap: figure color map
    """
    
    if normalise:
        arr = norm(arr)
        
    sns.set_context(context)
    sns.set_style('dark')
    
    fig, ax1 = plt.subplots()
    
    img = ax1.imshow(arr, cmap = cmap,
                     extent = [np.min(mesh[1])*scale, np.max(mesh[1])*scale,
                               np.min(mesh[0])*scale, np.max(mesh[0])*scale],
                     vmin = vmin, vmax = vmax,
                     aspect = aspect)
    
    ax1.set_title(title)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='7.5%', pad=0.05)

    cbar = fig.colorbar(img, cax)
    cbar.set_label(clabel)
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    
    ax1.annotate(label, horizontalalignment = 'left',   
                    verticalalignment = 'bottom',
                    xy = (0,1),
                    c = 'white')
    
    if sdir is None:
        fig.show()
    else:
        fig.savefig(sdir + ".png")
        plt.show()
        
def signal_plot(xdata, ydata,
        xlabel,
        ylabel,
        title,
        context = 'notebook',
        return_axes = False):

    sns.set_context(context)
    sns.set_style('dark')
    fig, ax1 = plt.subplots()
    
    ax1.plot(xdata, ydata)
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    
    if return_axes:
        return ax1
    else:
        plt.show()
        
def scatter_plot(xdata, ydata = None,
        xlabel = "",
        ylabel = "",
        title = "",
        context = 'notebook',
        return_axes = False):

    if ydata is None:
        ydata = np.arange(0, len(xdata), len(xdata))
                          
    sns.set_context(context)
    sns.set_style('dark')
    fig, ax1 = plt.subplots()
    
    ax1.scatter(xdata, ydata)
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    
    if return_axes:
        return ax1
    else:
        plt.show()
        