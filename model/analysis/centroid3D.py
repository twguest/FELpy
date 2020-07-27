#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:45:25 2020

@author: twguest
"""


import numpy as np

from matplotlib import pyplot as plt

from mpl_toolkits.mplot3d import Axes3D  

def plot(positions,t):
    
    fig = plt.figure()
    ax1 = fig.add_subplot(projection = '3d')
    
    for itr in range(len(positions)):
        
        ax1.scatter(positions[itr][0]*1e6, t[itr]*1e15, positions[itr][1]*1e6)
        
    ax1.set_xlabel("x [$\mu m$]")
    ax1.set_zlabel("y [$\mu m$]")
    ax1.set_ylabel("time [fs]")
    plt.show()

def testPlot():
    
    t = np.linspace(-4.5e-15, 4.5e-15, 100)
    
    centroids = []
    
    for itr in range(100):
        centroids.append(np.random.rand(2))
    
    plot(centroids, t)    
    
if __name__ == '__main__':
    
    testPlot()