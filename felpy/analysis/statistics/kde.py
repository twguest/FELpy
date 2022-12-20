"""
DEPR
"""

# -*- coding: utf-8 -*-
import numpy as np
import seaborn as sns
from scipy.stats import gaussian_kde



def plot_kde(data, mesh, title = None, lims = None, weights = None):
    

    sns.set()
    fig, ax1 = plt.subplots()
 
    ax1 = sns.kdeplot(data[:,0], data[:,1], shade = True, cbar = True,
                      common_norm=True, common_grid = True, weights = weights)
    
    if lims != None:
        ax1.set_xlim(lims[0][0], lims[0][1])
        ax1.set_ylim(lims[1][0], lims[1][1])
        
    ax1.set_xlabel("x [$\mu$m]")
    ax1.set_ylabel("y [$\mu$m]")
    
    if title != None:
        ax1.set_title(title)
    
    return fig
 

def get_kde_kernel(data, mesh, weights = None):
    
    X, Y = mesh
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = data.T
    kernel = gaussian_kde(values, weights = weights)

    return kernel

def get_kde(data, mesh, weights = None):
    
    X, Y = mesh
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = data.T
    kernel = gaussian_kde(values, weights = weights)
    kde = np.reshape(kernel(positions).T, X.shape)
    return kde 


if __name__ == '__main__':
    
    pass
    