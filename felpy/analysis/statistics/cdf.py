# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def get_cdf(data, mean = None):
    
    if mean == None:
        mean = np.mean(data)
        
    data = abs(data-mean)
    
    prob = ECDF(data[:,0]).y
    
    cdfx = ECDF(data[:,0]).x
    cdfy = ECDF(data[:,1]).x
    
    cdf = np.array([prob, cdfx, cdfy])
    
    return cdf [:,1:]

def plot_cdf_data(cdf, mode = '', method = np.mean):
    
    if mode == 'single':
        
        sns.set()
        fig, ax1 = plt.subplots()
        
        xcdf = ax1.plot(cdf[0], cdf[1]) ## x cdf
        ycdf = ax1.plot(cdf[0], cdf[2]) ## y cdf
        
        ax1.legend(["x axis", "y axis"])
        
        ax1.set_xlabel("Distance from Mean Centroid")
        ax1.set_ylabel("Probability")
        
        ax1.set_title("Cumulative Density Function")
        
        
    else:
        
        sns.set()
        fig, ax1 = plt.subplots()
        
        xcdf = ax1.plot(cdf[1,:,0], method(cdf[0,:,], axis = -1)) ## x cdf
        ycdf = ax1.plot(cdf[2,:,0], method(cdf[0,:,], axis = -1)) ## y cdf
        
        ax1.legend(["x axis", "y axis"])
        
        ax1.set_xlabel("Distance from Mean Centroid")
        ax1.set_ylabel("Probability")
        
        ax1.set_title("Empirical Cumulative Density Estimate")
        
    return fig
