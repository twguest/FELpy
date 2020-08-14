#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 13:41:39 2020

@author: twguest
"""

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

def scatter3(x, y, a, xlabel = "", ylabel = "", alabel = ""):
    """
    scatter 3 data sets (two spatial, and a third quantity q1)
    
    :param x: horizontal spatial array
    :param y: vertical spatial array
    :param a: tertiary function, reflected in marker size
    """
    sns.set()
    
    fig = plt.figure()

    ax1 = fig.add_subplot()
    
    ax1.scatter(x,y, s = a, label = 'Pulse Energy')
    
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    
    plt.legend()
    
    plt.show()
    
    

def scatter4(x, y, a, b, title = "", xlabel = "", ylabel = "", alabel = "", blabel = ""):
    """
    scatter 3 data sets (two spatial, and a third quantity q1)
    
    :param x: horizontal spatial array
    :param y: vertical spatial array
    :param a: tertiary function, reflected in marker size
    """
    sns.set()
    
    fig = plt.figure()

    ax1 = fig.add_subplot()
    
    im1 = ax1.scatter(x,y, s = a, c = b, label = 'Pulse Energy', cmap='jet')
    
    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    
    
    plt.legend()
    
    
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax = cax)
    cax.set_ylabel(blabel)
    plt.show()


if __name__ == '__main__':
    
    x = np.random.rand(1, 25)*10
    y = np.random.rand(1, 25)*10
    a = np.random.rand(1, 25)*100
    b = np.random.rand(1, 25)*10
    #scatter3(x,y,a, xlabel = "x", ylabel = "y")
    
    scatter4(x,y,a,b, blabel = "doodle")