# -*- coding: utf-8 -*-
"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.1"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""

import numpy as np

def wavefront_tilt(meshgrid, kx, ky):
    """ 
    define a two-dimensional complex tilted plane-wave, tilted along transverse pointing vectors kx and ky
    
    :param mesh: [2,nx,ny] array of coordinates, first-axis corresponds to x and y axis,
    see felpy.utils.np_utils.get_mesh for more
    :param kx: horizontal pointing vector [-pi, pi]
    :param ky: vertical pointing vector [-pi, pi]
    """
    
    return np.exp(1j*(kx*meshgrid[0]+ky*meshgrid[1]))