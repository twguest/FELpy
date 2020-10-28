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

def crop_array(arr, llim, ulim):
    
    coords = np.argwhere(np.logical_and(arr >= llim, arr <= ulim))
    x_min, y_min = coords.min(axis=0)
    x_max, y_max = coords.max(axis=0)
    cropped = arr[x_min:x_max+1, y_min:y_max+1]
    extent = [x_min, x_max, y_min, y_max]
    return cropped, extent

def extent_from_mesh(mesh, xlims, ylims):
    
    xm = mesh[0]
    ym = mesh[1]
    
    croppedx, extentx = crop_array(xm, *xlims)
    croppedy, extenty = crop_array(ym, *ylims)
    
    extent = [extentx[2], extentx[3], extenty[0], extenty[1]]
    
    return extent
    
    
if __name__ == '__main__':
    from matplotlib import pyplot as plt
    ii = np.random.rand(1000,1000)
    mesh = get_mesh(ii, 1, 1)
    
    arr = mesh[0]
    
    extent = extent_from_mesh(mesh, (10,1000), (49, 1003))
    plt.imshow(ii[extent[0]:extent[1], extent[2]:extent[3]])
    print(ii[extent[0]:extent[1], extent[2]:extent[3]].shape)
    