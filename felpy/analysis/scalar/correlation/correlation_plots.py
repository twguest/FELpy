# -*- coding: utf-8 -*-
import numpy as np
import seaborn as sns
from labwork.about import dCache
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from felpy.utils.vis_utils import animate
from felpy.utils.os_utils import mkdir_p 


def correlation_plot(corr, mesh, label = "", title = "", sdir = None, context = 'talk', 
                     cmap = 'jet'):
    """ 
    plot the correlation function of a pair of 2D arrays (x,y)
    
    :param corr: 2D correlation array (via get_correlation)
    :param mesh: coordinate mesh [np array]
    :param sdir: save directory for output .png
    :param label: figure label
    :param title: figure title
    :param cmap: figure color map
    """
    #sns.set({'axes.grid' : False})
    sns.set_context(context)
    fig, ax1 = plt.subplots()
    
    img = ax1.imshow(corr, cmap = cmap,
                     extent = [np.min(mesh[1])*1e6, np.max(mesh[1])*1e6,
                               np.min(mesh[0])*1e6, np.max(mesh[0])*1e6],
                    vmin = 0.5, vmax = 1.0)
    
    fig.suptitle = title
    
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='7.5%', pad=0.05)

    cbar = fig.colorbar(img, cax)
    cbar.set_label("$g^{(N)}_T$")
    
    ax1.set_xlabel("x [$\mu$m]")
    ax1.set_ylabel("y [$\mu$m]")
    
    ax1.annotate(label, horizontalalignment = 'left',   
                    verticalalignment = 'bottom',
                    xy = (0,1),
                    c = 'white')
    
    plt.rcParams["axes.grid"] = False
    
    if sdir is None:
        fig.show()
    else:
        fig.savefig(sdir)
        plt.show()
        
 
    
def correlation_lag_plot(corr, mesh, lag = 1,
                        title = "", xlabel = "",
                        ylabel = "",
                        label = "",
                        context = 'talk',
                        cmap = 'viridis',
                        sdir = None):
    
   
    fig, ax1 = plt.subplots(figsize = [5,5])

    ax1.annotate(label, horizontalalignment = 'left',   
                verticalalignment = 'bottom',
                xy = (0,1),
                c = 'black')
    

    rng = corr.shape[-1]-lag
    
    for itr in range(0,rng):
        if itr == 0:
            lc = corr[:,:,itr,itr + lag]
        else:
            lc += corr[:,:,itr,itr + lag]
    lc /= float(rng)
    print(np.min(lc))
    img = ax1.imshow(lc,cmap = cmap,
           extent = [np.min(mesh[1])*1e6, np.max(mesh[1])*1e6,
           np.min(mesh[0])*1e6, np.max(mesh[0])*1e6],
                    vmin = 0.75, vmax = 1.0)
    

    
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='7.5%', pad=0.05)

    cbar = fig.colorbar(img, cax)
    cbar.set_label("$g^{(N)}_T$")
    
       
    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    
    sns.set_context(context)
    
    if sdir is not None:
        fig.savefig(sdir + "{}.png".format(label))
    
def correlogram(corr,
                context = 'talk',
                title = "",
                xlabel = "",
                ylabel = ""):
    
    sns.set_context(context)
    
    fig, ax1 = plt.subplots()
    
    ax1.set_title(title)    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    
    img = ax1.imshow(corr)
    
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='7.5%', pad=0.05)

    cbar = fig.colorbar(img, cax)
    cbar.set_label("$g^{(N)}_T$")
    