#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 12:09:30 2020

@author: twguest
"""

###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################

import numpy as np
from scipy.ndimage import gaussian_filter
from matplotlib import pyplot as plt 

def genMirrorSurface(nx, ny, mirDim, outdir, mode = 'Flat', plot = False, mirrorName = None):
    """
    Generate a plane mirror surface
    
    :param nx: number of horizontal pixels [int]
    :param ny: number of vertical pixels [int]
    :param mirDim: list of mirror dimensions [dx,dy] [m]
    :param outdir: save directory
    :param mode: type of mirror surface to be generated
    """
    
    mirLen = mirDim[0]
    mirWid = mirDim[1]
    
    if mode == 'flat':
        surface = np.zeros((nx,ny))

    if mode == 'random':
        surface = np.random.normal(size = [nx,ny])*1e-09
        surface = gaussian_filter(surface, 1)
        
        if plot == True:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            
            if mirrorName is not None:
                ax.set_title(mirrorName + " Surface")
            elif mirrorName is None:
                ax.set_title("Mirror Surface")
            
            img = ax.imshow(surface*1e9,
                            extent = [-mirLen/2*1e6, mirLen/2*1e6, -mirWid/2*1e6, mirWid/2*1e6],
                            aspect = 'auto')
            
            ax.set_xlabel("x ($\mu$m)")
            ax.set_ylabel("y ($\mu$m)")
            
            cb = plt.colorbar(img, ax = ax)
            cb.ax.get_yaxis().labelpad = 15
            cb.ax.set_ylabel("Height Error (nm)", rotation = 270)
            
            fig.savefig(outdir + "mir_"+mode+".png")
            
    surface[1:, 0] = np.linspace(-mirLen/2, mirLen/2, nx-1) 
    surface[0, 1:] = np.linspace(-mirWid/2, mirWid/2, ny-1) 
    
    np.savetxt(outdir+"mir_"+ mode +".dat", surface, delimiter='\t')


            
if __name__ == '__main__':
    genMirrorSurface(100, 100, [1e-06, 1e-06], "../../tmp/", mode = 'random', plot = True)
        