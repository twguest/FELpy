#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.0"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""

import os
import numpy as np
from copy import copy

from PIL import Image


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

    dx, dy = wfr.get_spatial_resolution()

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


def read_map(mapdir, shape, dtype='float64'):
    """
    read a map from mapDir
    """

    mp = np.memmap(mapdir, dtype=dtype, mode='r+', shape=shape)

    return mp


def memory_map(map_loc, shape, dtype='float64'):
    """
    construct a memory map
    """

    if os.path.exists(map_loc):
        memmap = read_map(map_loc, shape, dtype)
    else:
        memmap = np.memmap(map_loc, mode="w+",
                           shape=shape, dtype=dtype)
    return memmap


def equate_mesh(mesh1, mesh2):
    """
    crops mesh2 to the size of mesh1
    """

    xmin, xmax = np.min(mesh1[0, ]), np.max(mesh1[0, ])
    ymin, ymax = np.min(mesh1[1, ]), np.max(mesh1[1, ])
    print(xmin, xmax, ymin, ymax)

    mesh3 = mesh2[1,
                  find_nearest_above(mesh2[1, 0, ], xmin): find_nearest_above(mesh2[1, 0, ], xmax),
                  find_nearest_above(mesh2[1, 0, ], xmin): find_nearest_above(mesh2[1, 0, ], xmax)]

    mesh4 = mesh2[0,
                  find_nearest_above(mesh2[0, 0, ], ymin): find_nearest_above(mesh2[0, 0, ], ymax)+1,
                  find_nearest_above(mesh2[0, 0, ], ymin): find_nearest_above(mesh2[0, 0, ], ymax)+1]

    mesh3 = np.array([mesh3, mesh4])

    return mesh3

    #mesh3[0,] = mesh2[0, np.where(mesh2[0,] == xmin), np.where(mesh2[0,] == xmin)]
    #mesh3[1,] = mesh2[1, np.where(mesh2[1,] == xmin), np.where(mesh2[1,] == xmin)]


def find_nearest_above(my_array, target):
    diff = my_array - target
    mask = np.ma.less_equal(diff, 0)
    # We need to mask the negative differences and zero
    # since we are looking for values above
    if np.all(mask):
        return None  # returns None if target is greater than any value
    masked_diff = np.ma.masked_array(diff, mask)
    return masked_diff.argmin()


def test_equate_mesh():
    a = np.random.rand(50, 50)
    b = np.random.rand(75, 75)

    da = 1
    db = 1.25

    mesh1 = get_mesh(a, da, da)
    mesh2 = get_mesh(b, db, db)

    mesh3 = equate_mesh(mesh1, mesh2)
    return mesh3


def crop_to(arr1, arr2):

    if arr1.shape[0] < arr2.shape[0]:
        arr2 = arr2[arr2.shape[0]//2 - arr1.shape[0]
                    // 2: arr2.shape[0]//2 + arr1.shape[0]//2]
    elif arr2.shape[0] < arr1.shape[0]:
        arr1 = arr1[arr1.shape[0]//2 - arr2.shape[0]
                    // 2: arr1.shape[0]//2 + arr2.shape[0]//2]
    if arr1.shape[1] < arr2.shape[1]:
        arr2 = arr2[:, arr2.shape[1]//2 - arr1.shape[1]
                    // 2: arr2.shape[1]//2 + arr1.shape[1]//2]
    elif arr2.shape[1] < arr1.shape[1]:
        arr1 = arr1[:, arr1.shape[1]//2 - arr2.shape[1]
                    // 2: arr1.shape[1]//2 + arr2.shape[1]//2]

    return arr1, arr2


def load_tif(directory):
    """
    load a tif file as an np array
    """
    return (np.asarray(Image.open(directory)))


if __name__ == '__main__':

    arr1 = np.random.rand(530, 200)
    arr2 = np.random.rand(406, 410)

    arr1, arr2 = crop_to(arr1, arr2)
    print(arr1.shape, arr2.shape)
