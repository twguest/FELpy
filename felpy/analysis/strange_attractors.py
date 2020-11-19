# -*- coding: utf-8 -*-
import numpy as np

from matplotlib import pyplot as plt
from labwork.about import dCache
from felpy.utils.os_utils import mkdir_p
from felpy.utils.vis_utils import animate
import seaborn as sns


def generate_attractor_animation(data,
                                 sdir,
                                 xlabel = "",
                                 ylabel = "",
                                 title = ""):
    
    
    sns.set()
    
    tdir = sdir + "/tmp/"
    mkdir_p(tdir)
    
    if data.ndim == 1:
        
        for n in range(data.shape[-1] - 1):
            
            fig, ax1 = plt.subplots(1,1)     
            
            ax1.plot(data[:n], np.roll(data, 1)[:n])
            ax1.set_xlabel(xlabel)
            ax1.set_ylabel(ylabel)
            ax1.set_title(title)
            

            
            fig.savefig(tdir + "{:04d}".format(n))
        
            
            
    elif data.ndim == 2 :
        
        for m in range(data.shape[-1]):
            
            for n in range(data.shape[1] - 1):
                fig, ax1 = plt.subplots()     
                ax1.plot(data[:n, m], np.roll(data, 1)[:n, m])
                ax1.set_xlabel(xlabel)
                ax1.set_ylabel(ylabel)
                ax1.set_title(title)
                
                fig.savefig(tdir + "{:4d}{:4d}".format(m,n))
            
                
        
    animate(tdir, sdir, title + "_attractor", rmdir = True)            

if __name__ == '__main__':
    
    fname = '/home/twguest/Downloads/r0033/centroid_data.npy'
    
    cdata = np.load('/home/twguest/Downloads/r0033/centroid_data.npy')
    
    generate_attractor_animation(cdata[:,0,0], sdir = 
                                 "/opt/labwork/scratch/")