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
    
def gaussian_2d(nx, ny):
    """
    generate a 2d gaussian-like beam in array form
    
    :param nx: number of columns
    :param ny: number of rows
    """
    x, y = np.meshgrid(np.linspace(-1,1,nx), np.linspace(-1,1,ny))
    d = np.sqrt(x*x+y*y)
    sigma, mu = 1.0, 0.0
    g = np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
    return g

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    g = gaussian_2d(10,25)
    plt.imshow(g)