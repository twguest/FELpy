import os 
import shutil
from karabo_data import RunDirectory
from felpy.utils.daq_utils import load_data, shimadzu_reshape
from extra_data import open_run, RunDirectory, H5File
import numpy as np
from felpy.utils.os_utils import mkdir_p

from felpy.analysis.enclosed_energy import get_enclosed_energy
from felpy.analysis.centroid import get_com
from labwork.about import dCache
from felpy.utils.np_utils import get_mesh
from felpy.utils.vis_utils import animate
from matplotlib import pyplot as plt
from felpy.analysis.enclosed_energy import get_enclosed_energy
from tqdm import tqdm
import seaborn as sns
from felpy.utils.vis_utils import Grids

def fluence(ii, w, E0):
    return
                
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
             return_grid = True,
             color_x = 'blue', color_y = 'orange',
             label = "", linestyle = '-',
             sdir = None, legend = True, context = 'notebook',
             avg = 'pulse', global_aspect = 2.5, report = True,
            parse_grid = None, xlim = None, ylim = None):
    
    plt.style.use(['science','ieee'])
    sns.set_context(context)
    
    if parse_grid is not None:
        plots = parse_grid
    else:
        plots = Grids(global_aspect = global_aspect, scale = 2)
        plots.create_grid(n = 1, m = 2, sharex = False, sharey = True, fontsize = 30)
        
    ax1,ax2 = plots.axes

    px = mesh[0][0,1]-mesh[0][0,0]
    py = mesh[1][1,0]-mesh[1][0,0]

    com = get_com(ii) ### because plot scale is in mm    
    com -= np.mean(com)
    
    if avg == 'pulse':
 
        std_x = com[:,0,:].std(axis = -1)*px
        std_y = com[:,1,:].std(axis = -1)*py

        mean_x = com[:,0,:].mean(axis = -1)*px - com[:,0,:].mean()*px
        mean_y = com[:,1,:].mean(axis = -1)*py - com[:,1,:].mean()*py
    

        ax1.plot(np.arange(com[:,0,:].mean(-1).shape[0]),mean_x, color = color_x, label = label, linestyle = linestyle)

        ax1.fill_between(np.arange(com[:,0,:].mean(-1).shape[0]),
                         mean_x - std_x,
                         mean_x + std_x,
                         alpha = 0.25,
                         color = color_x,
                         edgecolor= None, facecolor = None,
                         linewidth=0.0)
        
        ax2.plot(np.arange(com[:,0,:].mean(-1).shape[0]), -mean_y, color = color_y, label = label,  linestyle = linestyle)
    
        ax2.fill_between(np.arange(com[:,0,:].mean(-1).shape[0]),
                         -(mean_y - std_y),
                         -(mean_y + std_y),
                         alpha = 0.25,
                         color = color_y,
                         edgecolor= None, facecolor = None,
                         linewidth=0.0)
    
    if avg == 'train':
        
        std_x = com[:,0].std()*px
        std_y = com[:,1].std()*py

        mean_x = com[:,0]*px - com[:,0].mean()*px
        mean_y = com[:,1]*py - com[:,1].mean()*py

        ax1.plot(np.arange(com.shape[0]),mean_x, color = color_x, label = label, linestyle = linestyle)

        ax2.plot(np.arange(com.shape[0]), mean_y, color = color_y, label = label,  linestyle = linestyle)

    if legend:
        ax1.legend()

    if report:
        print("Horizontal C.O.M")
        #print("Mean Beam COM: {} mm".format(mean_x.mean()))
        print("Min/Max Beam Displacement: {:.2f}/{:.2f} mm".format(mean_x.min(), mean_x.max()))
        print("Variance in Beam Displacement: {} mm".format(mean_x.var()))
        print()
 
    

    
    if report:
        print("Vertical C.O.M")
        print("Min/Max Beam Displacement: {:.2f}/{:.2f} mm".format(mean_y.min(), mean_y.max()))
        print("STd in Beam Displacement: {} mm".format(mean_y.std()))
        print()
        
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("$\Delta$x (mm)")
    
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("$\Delta$ y (mm)")
    

 
    if legend:
        ax2.legend()
    
    if sdir is not None:
        plots.fig.savefig(sdir)

    if return_grid:
        return plots
    
    
def plot_enclosed_energy(ii, mesh, xlabel = "", ylabel = "", 
             parse_grid = None, linestyle = 'solid',
             color_x = 'blue', label = "", sdir = None, legend = True, context = 'talk', avg = 'train', global_aspect = 2.5,
                        return_grid = True, report = True, xlim = None, ylim = None, marker = None):
    
    plt.style.use(['science','ieee'])
    sns.set_context(context)
        

    if parse_grid is not None:
        plots = parse_grid
    else:
        plots = Grids(global_aspect = global_aspect, scale = 2)
        plots.create_grid(n = 1, m = 1, sharex = False, sharey = True, fontsize = 30)
    
    
    context = sns.set_context(context)
    ax1  = plots.axes

    
    px = mesh[0][0,1]-mesh[0][0,0]
    py = mesh[1][1,0]-mesh[1][0,0]
    
    if avg == 'pulse':
        
        ntrains = ii.shape[-1]
        npulses = ii.shape[-2]
    
        data = np.ones([ntrains, npulses, 3])
        
        for pulse in tqdm(range(npulses)):
            for train in range(ntrains):
                    results = get_enclosed_energy(ii[:,:,pulse,train], px, py)
                    data[train,pulse,0] = results[0][0]
                    data[train,pulse,1] = results[0][1]
                    data[train,pulse,2] = results[1]

        beam_area = data[:,:,0]
        err_area = data[:,:,2]

        mean_x = beam_area.mean(0)
        std_x = beam_area.std(0)

        mean_err = err_area.mean(0)
        
        ax1.fill_between(np.arange(npulses),
                     mean_x - std_x,
                     mean_x + std_x,
                     alpha = 0.25,
                     color = color_x,
                     edgecolor= None, facecolor = None,
                     linewidth=0.0)  
    

        ax1.errorbar(np.arange(npulses),mean_x, marker = marker, color = color_x, yerr = mean_err,
                     label = label, linestyle = linestyle,
                     elinewidth = .5, ecolor = color_x)

    elif avg == 'train':
        
        ntrains = ii.shape[-1]
    
        data = np.ones([ntrains, 3])

        for train in tqdm(range(ntrains)):
                results = get_enclosed_energy(ii[:,:,train], px, py)
                data[train,0] = results[0][0]
                data[train,1] = results[0][1]
                data[train,2] = results[1]

        beam_area = data[:,0]
        err_area = data[:,2]

        mean_x = beam_area-np.mean(beam_area)
        


        ax1.errorbar(np.arange(ntrains),beam_area, marker = marker, color = color_x, yerr = err_area,
                     label = label, linestyle = linestyle,
                     elinewidth = .5, ecolor = color_x)


    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("r ($mm$)")
    
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    
    if sdir is not None:
        plots.fig.savefig(sdir)
    
    if report:
        print("Mean Beam Area: {} mm".format(mean_x.mean()))
        print("Min/Max Beam Area: {:.2f}/{:.2f} mm".format(mean_x.min(), mean_x.max()))
        print("Error in Beam Area: {} mm".format(mean_x.std()))
    
    
    if legend:
        plt.legend()
     
    if return_grid:
        return plots
   



def plot_integrated_beam_intensity(ii, mesh, xlabel = "", ylabel = "", 
             parse_grid = None, return_grid = True, linestyle = 'solid',
             color = 'blue',
             label = "", sdir = None):
    
    plt.style.use(['science','ieee'])

    if parse_grid is None:
        fig,ax1 = plt.subplots()
    else:
        ax1 = parse_grid[0]
        
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
    
    if return_grid:
        return plots


def plot_xgm(data, xlabel = "", 
             parse_grid = None, return_axes = True, linestyle = 'solid',
             color = 'blue',
             label = "", sdir = None, global_aspect = 2.5,
            return_grid = True):
    
    plt.style.use(['science','ieee'])

    if parse_grid is not None:
        ax1 = parse_grid
    else:
        plots = Grids(global_aspect = global_aspect, scale = 2)
        plots.create_grid(n = 1, m = 1, sharex = False, sharey = True, fontsize = 20)
        ax1  = plots.axes
     
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

    if return_grid:
        return plots