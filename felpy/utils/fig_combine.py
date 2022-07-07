# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from felpy.utils.vis_utils import Grids

def combine_figures(F, orientation = 'v', dpi = 600, title = None):
    """
    combines a set of pre-defined matplotlib figures
    
    :param F: list of matplotlib.pyplot.fig objects
    :param orientation: define the orientation in which the figures will be stacked
    :param dpi: defines image resolution
    :param title: global title for the image
    
    :returns fig: concatenated figures
    
    """
    backend = mpl.get_backend()
    mpl.use('agg')

    canvas = []
    for f in F:
        canvas.append(f.canvas) 
        
    for c in canvas:
        c.draw()
    
    if orientation == 'h':
        a = np.hstack([np.array(c.buffer_rgba()) for c in canvas])
    if orientation == 'v':
        a = np.vstack([np.array(c.buffer_rgba()) for c in canvas])
        
    
    mpl.use(backend)
    fig,ax = plt.subplots(dpi=dpi)
    fig.suptitle(title)

    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)
    
    ax.matshow(a)
    return fig

if __name__ == '__main__':
    


    
    
    grid = Grids(global_aspect = 2.5)
    grid.create_grid(1,2)
    ax1,ax2 = grid.get_axes()
    ax1.plot(np.random.rand(100), color = 'b')
    ax2.plot(np.random.rand(100), color = 'g')
    f1 = grid.fig 
    
    grid = Grids(scale = 1)
    grid.create_grid(1,1)
    ax3 = grid.get_axes()
    ax3.plot(np.random.rand(100), color = 'r')
    f2 = grid.fig 
    
    combine_figures([f1,f2], orientation = 'v')
    
