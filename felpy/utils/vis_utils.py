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
import matplotlib 
import matplotlib.pyplot as plt

from matplotlib.colors import Normalize
import matplotlib.cm as cm

#!python numbers=disable
fig_width_pt = 400.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
fig_size = [fig_width,fig_width/2]

plt.style.use(['science','ieee'])
mpl.rcParams['axes.linewidth'] = 1 #set the value globally


 


class Grids:
    
    def __init__(self, fig_width = 400., global_aspect = 1, scale = 1, context = 'paper'):
        
        sns.set_context(context)
        
        fig_width_pt = 400.0 * scale  # Get this from LaTeX using \showthe\columnwidth
        inches_per_pt = 1.0/72.27               # Convert pt to inches
        fig_width = fig_width_pt*inches_per_pt  # width in inches
        self.fig_size = [fig_width,fig_width/global_aspect]


    def label(axis):
        pass
    
    def pad(self, pad):
        self.fig.tight_layout(pad = pad)

    def add_global_colorbar(self, clabel, cmap = 'bone', tick_values = None, tick_labels = None, fontsize = 12):
                
        cmap=cm.get_cmap(cmap)
        normalizer=Normalize(np.min(ii),np.max(ii))
        im=cm.ScalarMappable(norm=normalizer, cmap = cmap)
        
        if self.n*self.m > 1:
            cbar = self.fig.colorbar(im, ax=self.axes.ravel().tolist(), pad = 0.025)
        else:
            cbar = self.fig.colorbar(im, ax=self.axes, pad = 0.025)

        cbar.set_label(clabel, fontsize = fontsize)
        
        if tick_values != None and tick_labels != None:
            cbar.set_ticks(tick_values)
            cbar.set_ticklabels(tick_labels)
            
        
    def create_grid(self, n = 1, m = 1, title = None, xlabel = None, ylabel = None,
                            resolution = 100, fontsize = 12, sharex = True, sharey = True):
    
        self.n = n
        self.m = m
        
        fig, axes = plt.subplots(nrows=n, ncols=m,
                                 figsize = self.fig_size,
                                 dpi = resolution,
                                 sharex = sharex, sharey = sharey,  gridspec_kw={'width_ratios': np.ones(m)})
 
        
        if n*m > 1:
            
            for itr, axis in enumerate(axes.ravel()):
                                

                axis.set_aspect('auto')
        else:                            
                axes.set_aspect('auto')
                
 
        
        
        fig.text(-0.015, 0.5, ylabel, va='center', rotation='vertical', fontsize= fontsize)
        fig.text(0.475, 0, xlabel, ha='center', fontsize = fontsize)
       
        fig.subplots_adjust(right=0.5)
        
        fig.suptitle(title)    
        fig.tight_layout()

        self.fig = fig
        self.axes = axes

def colorbar_plot(dataset,
                  mesh = None,
                  label = None,
                  title = None,
                  xlabel = None,
                  ylabel = None,
                  clabel = "",
                  context = 'paper',
                  cmap = 'bone',
                  normalise = True,
                  scale = 1,
                  aspect = 'auto',
                  sdir = None, 
                  lognorm = False):
    
    """ 
    plot a 2D datasetay with a colorbar (x,y)
    
    :param corr: 2D correlation datasetay (via get_correlation)
    :param mesh: coordinate mesh [np datasetay]
    :param sdir: save directory for output .png
    :param label: figure label
    :param title: figure title
    :param cmap: figure color map
    """
    
    sns.set_context(context)
    
    if normalise:
        dataset = norm(dataset)
        vmin, vmax = 0,1
    else: 
        vmin, vmax = np.min(dataset), np.max(dataset)
    
    if mesh is not None:
        extent = [np.min(mesh[1])*scale, np.max(mesh[1])*scale,
                  np.min(mesh[0])*scale, np.max(mesh[0])*scale]
    else:
        extent = None
    
    if lognorm == True:
        lognorm = matplotlib.colors.LogNorm()
        vmin += 1e-100
    else: 
        lognorm = None
        
    fig, ax1 = plt.subplots(figsize = fig_size)
    
    img = ax1.imshow(dataset,
                     cmap = cmap,
                     extent = extent,
                     vmin = vmin, vmax = vmax,
                     aspect = aspect, 
                     norm = lognorm)
    

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='7.5%', pad=0.05)

    cbar = fig.colorbar(img, cax)

    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    cbar.set_label(clabel)
    ax1.annotate(label, horizontalalignment = 'left',   
                    verticalalignment = 'bottom',
                    xy = (0,1),
                    c = 'white')
   
    if sdir is None:
        fig.show()
    else:
        fig.savefig(sdir + ".png")

    plt.show()




       

def triple_colorbar(ii, title = None, xlabel = None, ylabel = None, clabel = None,
           extent = None, cmap = 'hot', vmin = None, vmax = None, savedir = None,
           cticks = None, clabels = None, resolution = 100):
    
    #sns.set_style('dark')
    sns.set_context('paper')
    
 
    
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize = fig_size, dpi = resolution, sharex = True, sharey = True)
    
    for itr in range(ii.shape[-1]):
        
        
        im = axes[itr].imshow(ii[:,:,itr], cmap = cmap, vmin = vmin, vmax = vmax)

        axes[itr].set_xlabel(xlabel)
        axes[itr].set_xlabel(ylabel)
        axes[itr].set_aspect(1.5)
    
    
    fig.text(-0.015, 0.5, 'common Y', va='center', rotation='vertical', fontsize= 12)
    fig.text(0.55, 0.06, 'common X', ha='center', fontsize = 12)
    fig.subplots_adjust(right=0.75)
    fig.suptitle(title)    
    fig.tight_layout()
    
    cbar_ax = fig.add_axes([.99, 0.158, 0.04, 0.687])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(clabel, fontsize = 16)
    
    if cticks != None and clabels != None:
        cbar.set_ticks(cticks)
        cbar.set_ticklabels(clabels)
    
    if savedir is not None:
        fig.savefig(savedir)
    else:
        plt.show()

def small_colorbar_grid(ii, title = None, xlabel = None, ylabel = None, clabel = None,
           extent = None, cmap = 'hot', vmin = None, vmax = None, savedir = None,
           cticks = None, clabels = None, resolution = 100):
    
    #sns.set_style('dark')
    sns.set_context('paper')
    
 
    
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize = fig_size, dpi = resolution, sharex = True, sharey = True)
    
    for itr in range(ii.shape[-1]):
        
        
        im = axes[itr].imshow(ii[:,:,itr], cmap = cmap, vmin = vmin, vmax = vmax)

        axes[itr].set_xlabel(xlabel)
        axes[itr].set_xlabel(ylabel)
        axes[itr].set_aspect(1.5)
    
    
    fig.text(-0.015, 0.5, 'common Y', va='center', rotation='vertical', fontsize= 12)
    fig.text(0.55, 0.06, 'common X', ha='center', fontsize = 12)
    fig.subplots_adjust(right=0.75)
    fig.suptitle(title)    
    fig.tight_layout()
    
    cbar_ax = fig.add_axes([.99, 0.158, 0.04, 0.687])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(clabel, fontsize = 16)
    
    if cticks != None and clabels != None:
        cbar.set_ticks(cticks)
        cbar.set_ticklabels(clabels)
    
    if savedir is not None:
        fig.savefig(savedir)
    else:
        plt.show()

def dataset2gif(fname, datasetay, fps=10, scale=1.0):
    """Creates a gif given a stack of images using moviepy
    see: https://gist.github.com/nirum/d4224ad3cd0d71bfef6eba8f3d6ffd59
    
    :params fname: save directory (str)
    :params datasetay: datasetay to be made into gif (along -1 axis)
    :params fps: frames per second (int)
    :params scale: rescaling factor of image
    """

    # ensure that the file has the .gif extension
    fname, _ = os.path.splitext(fname)
    filename = fname + '.gif'

    # normalise 
    for itr in range(datasetay.shape[-1]):
        datasetay[:,:,itr] = norm(datasetay[:,:,itr]*255)    

    # copy into the color dimension if the images are black and white
    if datasetay.ndim == 3:
        datasetay = datasetay[..., np.newaxis] * np.ones(3)
        
    # make the moviepy clip
    clip = ImageSequenceClip(list(datasetay), fps=fps).resize(scale)
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
               scale = 1e6,
               context = 'notebook'):
    
    """ 
    a simple plot of some two-dimensional intensity datasetay
    
    :param ii: 2D intensity datasetay
    :param mesh: coordinate mesh [np datasetay]
    """
    
    if scale == 1e6:
        mode = '[$\mu$m]'
        
    elif scale == 1e3:
        mode = '[mm]'
        
    elif scale == 1:
        mode = '[m]'
        
    else:
        mode = "[" + str(scale) + "m]"
        
    plt.style.use(['science','ieee'])
    sns.set_context(context)
    
    fig, ax1 = plt.subplots()
    
    if xlims is not None and ylims is not None:
        roi = extent_from_mesh(mesh, xlims, ylims)
        ii = ii[roi[0]:roi[1], roi[2]:roi[3]]
            
        img = ax1.imshow(ii, cmap = cmap,
                         extent = [mesh[1,roi[2],0]*scale, mesh[1,roi[3],0]*scale,
                                     mesh[0,0,roi[0]]*scale, mesh[0,0,roi[1]]*scale],
                        aspect = 'equal')
    else:
        img = ax1.imshow(ii, cmap = cmap,
                         extent = [np.min(mesh[1])*scale, np.max(mesh[1])*scale,
                                   np.min(mesh[0])*scale, np.max(mesh[0])*scale],
                        aspect = 'auto')
  
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
    

def simple_line_plot(x, y = None,
                     title = "",
                      xlabel = "",
                      ylabel = "",
                      context = 'notebook',
                      xlim = None,
                      ylim = None,
                      color = 'blue',
                      parse_axes = None,
                      return_axes = False,
                      label = ""):
    """
    plot the mean of a 2D data set with the 95% confidence interval filled
    """
    #sns.set()
    plt.style.use(['science','ieee'])
    sns.set_context(context)
    
    if parse_axes is None:    
        fig, ax1 = plt.subplots()
    else: 
        ax1 = parse_axes
        
    if y is not None:
        ax1.plot(x, y, color = color, label = label)
    else:    
        ax1.plot(x, color = color, label = label)

    if xlim is not None:
        ax1.set_xlim(xlim)
    if ylim is not None:
        ax1.set_ylim(ylim)

    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    
    if return_axes:
        return ax1
    else:
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
        parse_axes = None,
        return_axes = False,
        marker= 'o',
        marker_size = 12,
        color = 'b',
        legend = False,
        legend_title = "",
        label = None):
    
    plt.style.use(['science','ieee'])

    if ydata is None:
        ydata = np.arange(0, len(xdata), len(xdata))
                  

    if parse_axes is not None:
        ax1 = parse_axes
    else:
        fig, ax1 = plt.subplots()
    
    ax1.scatter(xdata, ydata, label = label, marker = marker, s = marker_size)
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    
    if legend:
        ax1.legend(title=legend_title)
    #ax1.autoscale(tight=True) 
    #sns.set_context(context)
    if return_axes:
        return ax1
    else:
        plt.show()
        
def double_colorbar_plot(dataset1, dataset2,
                         extent1 = None,
                         extent2 = None,
                         xlabel1 = "",
                         xlabel2 = "",
                         ylabel1 = "",
                         ylabel2 = "",
                         clabel1 = "",
                         clabel2 = "",
                         title1 = "",
                         title2 = "",
                         cmap1 = 'bone',
                         cmap2 = 'hsv',
                         vmin1 = None,
                         vmax1 = None,
                         vmin2 = None,
                         vmax2 = None,
                         aspect = 'auto',
                         context = 'talk',
                         sdir = None):
    
    sns.set_style('dark')
    sns.set_context(context)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    
    im1 = ax1.imshow(dataset1, interpolation='None', extent = extent1,
                     cmap = cmap1, vmin = vmin1, vmax = vmax1)
    
    ax1.set_xlabel(xlabel1)
    ax1.set_ylabel(ylabel1)
    ax1.set_title(title1)
    
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='7.5%', pad=0.05)
    
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
    cbar1.set_label(clabel1)

    ax2 = fig.add_subplot(122)
    im2 = ax2.imshow(dataset2, interpolation='None', extent = extent2,
                     cmap = cmap2, vmin = vmin2, vmax = vmax2)
    ax2.set_xlabel(xlabel2)
    ax2.set_ylabel(ylabel2)
    ax2.set_title(title2)
    
    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes('right', size='7.5%', pad=0.05)
    
    cbar2 = fig.colorbar(im2, cax=cax2, orientation='vertical')
    cbar2.set_label(clabel2)
    #fig.tight_layout()
    #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=None)

    if sdir is not None:
        fig.savefig(sdir)
        
def contour_plot(x,y,z,
                 xlabel = "",
                 ylabel = "",
                 title = "",
                 clabel = "",
                 aspect = 'auto',
                 return_axes = False,
                 cbar_scale = 1,
                 context = 'notebook',
                 cmap = 'viridis'):
    
    fig, ax1 = plt.subplots()
    plt.style.use(['science','ieee'])
    sns.set_context(context)

        
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(title)
    
    im1 = ax1.contourf(x,y,z*cbar_scale,
                       cmap = cmap)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='7.5%', pad=0.05)
    
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
    cbar1.set_label(clabel, fontsize = 12)
    
    ax1.set_aspect(aspect)
    
    if return_axes:
        return ax1


if __name__ == '__main__':
    
    from scipy.ndimage import gaussian_filter
    #matplotlib.rcParams.update({'font.size': 22})
# =============================================================================
#     ii = gaussian_filter(np.random.rand(250,250,3),8)
#     triple_colorbar(ii, cmap = 'bone')
#     
#     ii = np.random.rand(250,250,2)
#     double_colorbar(ii)
# 
#     
#     ii = np.random.rand(250,250,4)
#     small_colorbar_grid(ii)
# =============================================================================
    CMAP = 'hsv'
    ii = gaussian_filter(np.random.rand(2000,780,3),8)
    plots = Grids(fig_width = 400, global_aspect = 1.5)
    plots.create_grid(n = 1, m = 3, xlabel = "x", ylabel = "y")
    plots.axes[0].imshow(ii[:,:,0], cmap = CMAP)
    plots.axes[1].imshow(ii[:,:,1], cmap = CMAP)
    plots.axes[2].imshow(ii[:,:,2], cmap = CMAP)
    
    plots.add_global_colorbar("treys colorbar", cmap = CMAP)
    
    
    n = 2
    m = 2
    
    CMAP = 'hsv'
    plots = Grids(fig_width = 400, global_aspect = 1.5)
    plots.create_grid(n = n, m = m, xlabel = "x", ylabel = "y",)
    
    idx = np.unravel_index(np.arange(n*m), [n,m])
    
    plots.axes[idx[0][0], idx[0][1]].plot(np.arange(20))
    
    plots.axes[idx[0][0], idx[0][1]].set_xlabel("trey")
    plots.pad(3)
    