#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.1"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""

import shutil 
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import string
from sklearn.preprocessing import minmax_scale as norm
from felpy.utils.os_utils import mkdir_p


from mpl_toolkits.mplot3d import Axes3D
from felpy.utils.np_utils import extent_from_mesh
import matplotlib 
import matplotlib.pyplot as plt
from PIL import ImageColor
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

exfel_c1 = [i/256 for i in ImageColor.getrgb("#0D1546")]
exfel_c2 = [i/256 for i in ImageColor.getrgb("#F39200")]
exfel_c3 = [i/256 for i in ImageColor.getrgb("#559DBB")]
exfel_c4 = [i/256 for i in ImageColor.getrgb("#B2B2B2")]

#!python numbers=disable
fig_width_pt = 400.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
fig_size = [fig_width,fig_width/2]

plt.style.use(['science','ieee'])
mpl.rcParams['axes.linewidth'] = 1 #set the value globally


from mpl_toolkits.axes_grid1.inset_locator import InsetPosition,inset_axes
import matplotlib 
from matplotlib import pyplot as plt
import numpy as np

def complex_array_to_rgb(X, theme='dark', rmax=None):
    '''
    via: https://stackoverflow.com/questions/15207255/is-there-any-way-to-use-bivariate-colormaps-in-matplotlib

    Takes an array of complex number and converts it to an array of [r, g, b],
    where phase gives hue and saturaton/value are given by the absolute value.
    Especially for use with imshow for complex plots.
    '''
    absmax = rmax or np.abs(X).max()
    Y = np.zeros(X.shape + (3,), dtype='float')
    Y[..., 0] = np.angle(X) / (2 * np.pi) % 1
    if theme == 'light':
        Y[..., 1] = np.clip(np.abs(X) / absmax, 0, 1)
        Y[..., 2] = 1
    elif theme == 'dark':
        Y[..., 1] = 1
        Y[..., 2] = np.clip(np.abs(X) / absmax, 0, 1)
    Y = matplotlib.colors.hsv_to_rgb(Y)
    return Y


def grad_to_rgb(angle, absolute):
    """Get the rgb value for the given `angle` and the `absolute` value

    Parameters
    ----------
    angle : float
        The angle in radians
    absolute : float
        The absolute value of the gradient
    
    Returns
    -------
    array_like
        The rgb value as a tuple with values [0..1]
    """
    global max_abs

    # normalize angle
    angle = angle % (2 * np.pi)

    if angle < 0:
        angle += 2 * np.pi

    return matplotlib.colors.hsv_to_rgb((angle / 2 / np.pi, 
                                         absolute / max_abs, 
                                         absolute / max_abs))



def complex_plot(wfr, ax, scale = 1, rmax = None):
    """ 
    function to generate a complex-valued plot of a wpg.wavefront
    
    :param wfr: felpy.model.wavefront type object to be plotted
    :param ax: axis to be plotted on
    :param scale: scale of extent (i.e. 1 corresponds to m, 1e3 to mm, 1e-03 to km)
    :param rmax: maximum value of intensity colormap
    """
    ii = wfr.get_intensity().sum(-1)
    ph = wfr.get_imag_part().sum(-1)

    rmax = rmax or np.max(ii)
    
    data = ii + 1j*ph
    ax.imshow(complex_array_to_rgb(data, theme = 'dark', rmax = rmax), extent = [scale*l for l in wfr.extent])


def add_inset_polar_map(grid, rmax = 1, inset_x = 0.8, inset_y = 0.65, inset_size = 0.5):
    """ 
    add a polar colorbar to a currently existing axis
    saves axes to felpy grid utility 

    :param grid: grid to add colorbar to
    :param rmax: value corresponding to maximum brightness value
    """
    ax = grid.fig.add_subplot(1, 1, 1, projection='polar')
    
    n = 200  # the number of secants for the mesh
    t = np.linspace(0, 2 * np.pi, n)
    r = np.linspace(0, rmax, n)
    rg, tg = np.meshgrid(r, t)
    c = np.array(list(map(grad_to_rgb, tg.T.flatten(), rg.T.flatten())))
    cv = c.reshape((n, n, 3))

    m = ax.pcolormesh(t, r, cv[:,:,1], color=c, shading='auto')

    m.set_array(None)
    ax.set_yticklabels([])

    ip = InsetPosition(ax, [inset_x - inset_size / 2,
                            inset_y - inset_size / 2,
                            inset_size,
                            inset_size])

    ax.set_axes_locator(ip)
    ax.tick_params(colors='white', which='both')
    ax.set_xticklabels(["0,2$\pi$","$\pi/4$","$\pi$/2", "3$\pi$/4", "$\pi$", "5$\pi$/4", "$3\pi/2$","7$\pi$/4"], fontsize = 10)
    
    grid.axes['cb'] = ax

def exfel_cmap():

    N = 256
    c1 = np.ones([N, 4])

    c1[:, 0] = np.linspace(exfel_c1[0], 1, N)  
    c1[:, 1] = np.linspace(exfel_c1[1], 1, N)  
    c1[:, 2] = np.linspace(exfel_c1[2], 1, N)  

    cmap_1 = ListedColormap(c1)

    c2 = np.ones([N, 4])

    c2[:, 0] = np.linspace(exfel_c2[0], 1, N)  
    c2[:, 1] = np.linspace(exfel_c2[1], 1, N)  
    c2[:, 2] = np.linspace(exfel_c2[2], 1, N)  

    cmap_2 = ListedColormap(c2)

    newcolors2 = np.vstack((cmap_1(np.linspace(0, 1, 128)),
                            cmap_2(np.linspace(1, 0, 128))))

    cmap = ListedColormap(newcolors2, name='exfel')

    return cmap  

def add_colorbar(mappable, ax, fig,
                 orientation = 'vertical',
                 cmap = 'bone',
                 vmin = 0,
                 vmax = 1,
                 size = '5%',
                 pad = 0.25,
                 clabel = "",
                 fontsize = 16,
                 labelpad = 0.05):
    """
    add a colorbar to an axis
    """
    divider = make_axes_locatable(ax)
    
    if orientation == 'horizontal':
        cax = divider.new_vertical(size = size,
                                    pad = pad,
                                    pack_start = True)
    if orientation == 'vertical':
        cax = divider.new_horizontal(size = size,
                                     pad = pad,
                                     pack_start = False)
    fig.add_axes(cax)
    cbar = fig.colorbar(mappable, cax = cax, orientation = orientation, cmap = cmap)
    cbar.set_label(clabel, fontsize = fontsize, labelpad = labelpad)
    return cbar

class Grids:
    
    def __init__(self, fig_width = 400., global_aspect = 1, scale = 1, context = 'paper'):
        
        sns.set_context(context)
        
        fig_width_pt = 400.0 * scale  # Get this from LaTeX using \showthe\columnwidth
        inches_per_pt = 1.0/72.27               # Convert pt to inches
        fig_width = fig_width_pt*inches_per_pt  # width in inches
        self.fig_size = [fig_width,fig_width/global_aspect]
        self.fontsize = 16

    def label(axis):
        pass
    
    def pad(self, pad):
        self.fig.tight_layout(pad = pad)

    def add_global_colorbar(self, clabel, vmin = 0, vmax = 1, cmap = 'bone', tick_values = None, tick_labels = None, fontsize = 12, orientation = 'vertical', pad = 0.025):
                
        cmap=cm.get_cmap(cmap)
        normalizer=Normalize(vmin, vmax)
        im=cm.ScalarMappable(norm=normalizer, cmap = cmap)
        
        if type(self.axes) == np.ndarray:
            if len(self.axes.flatten()) > 1:
                self.cbar = self.fig.colorbar(im, ax=self.axes.ravel().tolist(), pad = pad, orientation = orientation)
            else:
                self.cbar = self.fig.colorbar(im, ax=self.axes, pad = pad, orientation = orientation)
    


            
        elif type(self.axes) == dict:
            aa = [self.axes[key] for key in self.axes.keys()]
            self.cbar = self.fig.colorbar(im, ax=aa, pad = pad, orientation = orientation)
        else:
            self.cbar = self.fig.colorbar(im, ax=self.axes, pad = pad, orientation = orientation)

        self.cbar.set_label(clabel, fontsize = fontsize)
        self.cbar.set_label(clabel, fontsize = fontsize)

        if tick_values != None and tick_labels != None:
            self.cbar.set_ticks(tick_values)
            self.cbar.set_ticklabels(tick_labels)
            
            
            
    def get_axes(self):
        
        if type(self.axes) == np.ndarray:
            a = self.axes.flatten()
        elif type(self.axes) == dict:
            a = self.axes
        return a
    
    def set_fontsize(self, fontsize = 16):
        
        self.fontsize = fontsize 
        
        if type(self.axes) == np.ndarray:
            if len([*self.axes]) == 1: 
                self.axes.tick_params(axis='both', which='major', labelsize=fontsize)
                self.axes.xaxis.label.set_size(fontsize)
                self.axes.yaxis.label.set_size(fontsize)
            else:
                for ax in self.axes.flatten():
                    ax.tick_params(axis='both', which='major', labelsize=fontsize)
                    ax.xaxis.label.set_size(fontsize)
                    ax.yaxis.label.set_size(fontsize)
                    
        elif type(self.axes) == dict:
            for key in self.axes.keys():
                self.axes[key].tick_params(axis='both', which='major', labelsize=fontsize)
                self.axes[key].xaxis.label.set_size(fontsize)
                self.axes[key].yaxis.label.set_size(fontsize)
    
        else:
            self.axes.tick_params(axis='both', which='major', labelsize=fontsize)
            self.axes.xaxis.label.set_size(fontsize)
            self.axes.yaxis.label.set_size(fontsize)
        
        if hasattr(self, "cbar"):
            self.cbar.ax.tick_params(labelsize=fontsize)
            
    def savefig(self, sdir):
        self.fig.savefig(sdir, dpi = 600)
        
        
    def create_grid(self, n = 1, m = 1, title = None, xlabel = None, ylabel = None,
                            resolution = 100, fontsize = 12, sharex = True, sharey = True, **kwargs):
    
        self.n = n
        self.m = m
        
        if 'width_ratios' in kwargs:
            width_ratios = kwargs['width_ratios']
        else:
            width_ratios = np.ones(m)
        
        if 'height_ratios' in kwargs:
            height_ratios = kwargs['height_ratios']
        else:
            height_ratios = np.ones(n)
            
        fig, axes = plt.subplots(nrows=n, ncols=m,
                                 figsize = self.fig_size,
                                 dpi = resolution,
                                 sharex = sharex, sharey = sharey,  gridspec_kw={'width_ratios': width_ratios,
                                                                                'height_ratios': height_ratios})
 
        
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
        
        
    def create_mosaic(self, mosaic, title = None, xlabel = None, ylabel = None,
                            resolution = 100, fontsize = 12, sharex = False, sharey = False, **kwargs):
    
        self.mosaic = mosaic
        
        if 'width_ratios' in kwargs:
            width_ratios = kwargs['width_ratios']
        else:
            width_ratios = None
        
        if 'height_ratios' in kwargs:
            height_ratios = kwargs['height_ratios']
        else:
            height_ratios = None
        
        fig, axes = plt.subplot_mosaic(mosaic,
                                 figsize = self.fig_size,
                                 dpi = resolution,
                                 sharex = sharex, sharey = sharey, gridspec_kw={'width_ratios': width_ratios,
                                                                                'height_ratios': height_ratios})
    

            
        fig.text(-0.015, 0.5, ylabel, va='center', rotation='vertical', fontsize= fontsize)
        fig.text(0.475, 0, xlabel, ha='center', fontsize = fontsize)
       
        fig.subplots_adjust(right=0.5)
        
        fig.suptitle(title)    
        fig.tight_layout()

        self.fig = fig
        self.axes = axes
        
    def annotate(self, characters = string.ascii_lowercase, color = 'black', ha = 'left', va = 'bottom', bg = None, fontsize = 16):
        """
        in progress
        
        annotation function for sets of figures
        """
        
        if type(self.get_axes()) == np.ndarray:
            for itr in range(self.get_axes().shape[0]):
                ax = self.get_axes()[itr]
                ax.annotate(characters[itr],
                            horizontalalignment = ha,   
                            verticalalignment = va,
                            xy = (ax.get_xlim()[0],ax.get_ylim()[1]),
                            c = color,
                            fontsize = fontsize)
                
        elif type(self.axes) == dict:
            
            for itr, key in enumerate([*self.axes.keys()]):
                ax = self.get_axes()[key]

                ax.annotate(key,
                            horizontalalignment = ha,   
                            verticalalignment = va,
                            xy = (ax.get_xlim()[0],ax.get_ylim()[1]),
                            c = color,
                            fontsize = fontsize)


import ipywidgets as widgets
from ipywidgets import interactive, FloatSlider, IntSlider, HBox, Label

def single_plot(pulse,train, scale = 2, global_aspect = 1.25):
    """ 
    for use with viewer and mhz viewer class
    """
    grid = Grids(scale = scale, global_aspect = global_aspect)
    grid.create_mosaic([['a']])
    grid.axes['a'].imshow(image[:,:,pulse,train],vmin=0, cmap= 'afmhot', extent = [np.min(mesh[1])*1e3, np.max(mesh[1])*1e3,
                                   np.min(mesh[0])*1e3, np.max(mesh[0])*1e3])
    
    grid.axes['a'].set_xlabel("x (mm)")
    grid.axes['a'].set_ylabel("y (mm)")
    
    grid.add_global_colorbar(cmap = 'afmhot', clabel = "Intensity", fontsize = 22)
    grid.set_fontsize(22)


def mhz_viewer(image, mesh, scale = 2, global_aspect = 1.25):

    a = widgets.IntSlider(min=0, max= image.shape[-2]-1, step=1, value = 10, description='Pulse #:', continous_update = True, orientation='horizontal')
    b = widgets.IntSlider(min=0, max=image.shape[-1]-1, step=1, value = 10, description='Train #:', continous_update = True, orientation='horizontal')
    a.style.handle_color="white"
    b.style.handle_color="white"
    widgetsDisplay = widgets.HBox([a,b])

    interactive_plot = widgets.interactive_output(single_plot, {'pulse': a,'train': b})

    display(widgetsDisplay, interactive_plot)
    



def animate(indir, outdir, fname, delay = 0.1, rmdir = False):
    """
    create gif from pngs in directory
    """        
    os.system("convert -delay {} {}/*.png {}{}.gif".format(delay, indir, outdir, fname) )
    
    if rmdir == True:
        shutil.rmtree(indir)


if __name__ == '__main__':
    
    grid = Grids()
    grid.create_mosaic(mosaic = [['a','a','a'],['b','c','d']], sharex = True)
    grid.set_fontsize(15)
    grid.add_global_colorbar(clabel = 'trey', fontsize = 15)
    grid.annotate()
    
    