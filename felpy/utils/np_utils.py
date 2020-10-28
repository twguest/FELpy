# -*- coding: utf-8 -*-
import numpy as np

def get_mesh(ii, dx, dy):
    """
    returns a coordinate mesh ie, conversion of pixels to r-space
    
    :param ii: intensity array
    :param dx: horizontal pixel width
    :param dy: vertical pixel width
    
    :returns grid: pixel coordinate grid (2, nx, ny) [numpy array],
    where first dimension gives x- and y- position vector
    """
    
    xc = np.arange(ii.shape[0])*dx
    yc = np.arange(ii.shape[1])*dy
    
    grid = np.array(np.meshgrid(xc, yc))
    
    return grid

def get_wpg_mesh(wfr):
    
    dx, dy = wfr.pixelsize()
    
    nx, ny = wfr.params.Mesh.nx, wfr.params.Mesh.ny
    
    ii = np.zeros([nx, ny])
    
    mesh = get_mesh(ii, dx, dy)
    
    return mesh

