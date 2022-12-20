import os 
import shutil
from karabo_data import RunDirectory
from felpy.utils.daq_utils import load_data, shimadzu_reshape
from extra_data import open_run, RunDirectory, H5File
import numpy as np
from felpy.utils.os_utils import mkdir_p
from felpy.utils.vis_utils import  plot_fill_between, simple_line_plot, colorbar_plot
from felpy.analysis.optics.enclosed_energy import get_enclosed_energy
from felpy.analysis.optics.centroid import get_com
from labwork.about import dCache
from felpy.utils.np_utils import get_mesh
from felpy.utils.vis_utils import animate
from felpy.utils.vis_utils import plot_fill_between
from matplotlib import pyplot as plt
from felpy.analysis.optics.enclosed_energy import get_enclosed_energy
from tqdm import tqdm


def create_gif(ii, mesh, sdir, fname, title = "", xlabel = "", ylabel = ""):
    plt.style.use(['science','ieee'])
    mkdir_p(sdir)
    
    tmp_dir = sdir + "/tmp/"
    mkdir_p(tmp_dir)
        
    for i in range(ii.shape[-2]):
        colorbar_plot(ii[:,:,i,:].mean(-1),
                     mesh, 
                     title = title,
                     xlabel = xlabel,
                     ylabel = ylabel,
                     clabel = "Intensity (a.u.)",
                     normalise = False,
                     context = 'talk',
                     cmap = 'jet',
                     label = "pulse: {}".format(i+1),
                     sdir = tmp_dir + "pulse_{:0>4d}".format(i))

    animate(tmp_dir, sdir, fname, delay = 0.1, rmdir = True)
    
def plot_com(ii, mesh, xlabel = "", ylabel = "", 
             parse_axes = None, return_axes = True,
             color_x = 'blue', color_y = 'orange',
             label = "", linestyle = '-', sdir = None):
    
    plt.style.use(['science','ieee'])

    px = mesh[0][0,1]-mesh[0][0,0]
    py = mesh[1][1,0]-mesh[1][0,0]

    com = get_com(ii)    
    
    std_x = com[:,0,:].std(axis = -1)*px
    std_y = com[:,1,:].std(axis = -1)*py

    mean_x = com[:,0,:].mean(axis = -1)*px - com[:,0,:].mean()*px
    mean_y = com[:,1,:].mean(axis = -1)*py - com[:,1,:].mean()*py
   
    if parse_axes is not None:
        ax1 = parse_axes[0]
    else:
        fig1, ax1 = plt.subplots()

    
    
    ax1.plot(np.arange(com.shape[0]),mean_x, color = color_x, label = label, linestyle = linestyle)
    
    ax1.fill_between(np.arange(com.shape[0]),
                     mean_x - std_x,
                     mean_x + std_x,
                     alpha = 0.5,
                     color = color_x,
                     edgecolor= None, facecolor = None,
                     linewidth=0.0)
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("x (pixels)")
    
    ax1.legend()
 
    if sdir is not None:
        plt.savefig(sdir + "_a")

    if parse_axes is not None:
        ax2 = parse_axes[1]
    else:
        fig2, ax2 = plt.subplots()


    ax2.plot(np.arange(com.shape[0]), mean_y, color = color_y, label = label,  linestyle = linestyle)
    
    ax2.fill_between(np.arange(com.shape[0]),
                     mean_y - std_y,
                     mean_y + std_y,
                     alpha = 0.5,
                     color = color_y,
                     edgecolor= None, facecolor = None,
                     linewidth=0.0)
    
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("y (pixels)")
    
    ax2.legend()
    
    if sdir is not None:
        plt.savefig(sdir + "_b")
         
    if return_axes:
        return ax1, ax2
    
def plot_enclosed_energy(ii, mesh, xlabel = "", ylabel = "", 
             parse_axes = None, return_axes = True, linestyle = 'solid',
             color_x = 'blue', color_y = 'orange',
             label = "", sdir = None):
    plt.style.use(['science','ieee'])

        
    px = mesh[0][0,1]-mesh[0][0,0]
    py = mesh[1][1,0]-mesh[1][0,0]
    
    ntrains = ii.shape[-1]
    npulses = ii.shape[-2]
    
    data = np.ones([ntrains, npulses, 3])
    
    for train in tqdm(range(ntrains)):
        for pulse in range(npulses):
            results = get_enclosed_energy(ii[:,:,pulse,train], px, py)
            data[train,pulse,0] = results[0][0]
            data[train,pulse,1] = results[0][1]
            data[train,pulse,2] = results[1]
    
    x_area = data[:,:,0]
    y_area = data[:,:,1]
    err_area = data[:,:,2]

    mean_x = x_area.mean(0)
    std_x = x_area.std(0)
    
    mean_y = y_area.mean(0)
    std_y = y_area.std(0)
    
    mean_err = err_area.mean(0)

    if parse_axes is not None:
        ax1 = parse_axes[0]
    else:
        fig, ax1 = plt.subplots()

    
    ax1.errorbar(np.arange(mean_x.shape[0]),mean_x, color = color_x, yerr = mean_err,
                 label = label, linestyle = linestyle,
                 elinewidth = .1, ecolor = 'black')
    
    ax1.fill_between(np.arange(mean_x.shape[0]),
                     mean_x - std_x,
                     mean_x + std_x,
                     alpha = 0.5,
                     color = color_x,
                     edgecolor= None, facecolor = None,
                     linewidth=0.0)
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("x (pixels)")
    
    plt.legend()
 
    if sdir is not None:
        plt.savefig(sdir + "_a")
    

    
    if parse_axes is not None:
        ax2 = parse_axes[1]
    else:
        fig, ax2 = plt.subplots()
    

    
    ax2.errorbar(np.arange(mean_y.shape[0]), mean_y, yerr = mean_err,
                 color = color_y, label = label,  linestyle = linestyle,
                 elinewidth = .1, ecolor = 'black')
    
    ax2.fill_between(np.arange(mean_y.shape[0]),
                     mean_y - std_y,
                     mean_y + std_y,
                     alpha = 0.5,
                     color = color_y,
                     edgecolor= None, facecolor = None,
                     linewidth=0.0)
    
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("y (pixels)")
    
    plt.legend()
    
    
    if sdir is not None:
        plt.savefig(sdir + "_b")
        
    if return_axes:
        return ax1, ax2


def plot_integrated_beam_intensity(ii, mesh, xlabel = "", ylabel = "", 
             parse_axes = None, return_axes = True, linestyle = 'solid',
             color = 'blue',
             label = "", sdir = None):
    
    plt.style.use(['science','ieee'])

    if parse_axes is None:
        fig,ax1 = plt.subplots()
    else:
        ax1 = parse_axes[0]
        
    px = mesh[0][0,1]-mesh[0][0,0]
    py = mesh[1][1,0]-mesh[1][0,0]
    
    ntrains = ii.shape[-1]
    npulses = ii.shape[-2]
    
    data = np.ones([ntrains, npulses, 3])
    
    intensity = ii.sum(0).sum(0).mean(-1)
    std_dev = ii.sum(0).sum(0).std(-1)
    
    ax1.plot(np.arange(intensity.shape[0]), intensity, color = color,
                 label = label, linestyle = linestyle)
    
    ax1.fill_between(np.arange(intensity.shape[0]),
                     intensity-std_dev,
                     intensity + std_dev,
                     alpha = 0.5,
                     color = color,
                     edgecolor= None, facecolor = None,
                     linewidth=0.0)
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("$\sum$I (a.u.)")
    
    plt.legend()
 
    if sdir is not None:
        plt.savefig(sdir)
    

    if return_axes:
        return ax1


def plot_xgm(data, xlabel = "", 
             parse_axes = None, return_axes = True, linestyle = 'solid',
             color = 'blue',
             label = "", sdir = None):
    
    plt.style.use(['science','ieee'])

    if parse_axes is None:
        fig,ax1 = plt.subplots()
    else:
        ax1 = parse_axes[0]
 
    
     
    intensity = data.mean(0)
    std_dev = data.std(0)
    
    ax1.plot(np.arange(intensity.shape[0]), intensity, color = color,
                 label = label, linestyle = linestyle)
    
    ax1.fill_between(np.arange(intensity.shape[0]),
                     intensity-std_dev,
                     intensity + std_dev,
                     alpha = 0.5,
                     color = color,
                     edgecolor= None, facecolor = None,
                     linewidth=0.0)
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("XGM Intensity (a.u.)")
    
    plt.legend()
 
    if sdir is not None:
        plt.savefig(sdir)
    

    if return_axes:
        return ax1
