# -*- coding: utf-8 -*-

""" 
a script to test the com functions
"""

import numpy as np

from felpy.analysis.optics.scalar.enclosed_energy_DEPR import get_enclosed_energy as depr
from felpy.analysis.optics.scalar.enclosed_energy import get_enclosed_energy as get_enclosed_energy

from felpy.utils.os_utils import timing


@timing 
def old(ii, dx, dy):
    return depr(ii, dx, dy)[0]

@timing 
def new(ii, dx, dy):
    return get_enclosed_energy(ii, dx, dy)[0]
    
def speed_test(ii,dx,dy):
    
    for i in range(5):
        ans1 = old(ii,dx,dy)
        ans2 = new(ii,dx,dy)

        print(ans1)
        print(ans2)
        print("")
        
    
if __name__ == '__main__':
    
    from scipy.ndimage import gaussian_filter
    from matplotlib import pyplot as plt
    
    nx, ny = 500,500
    ii = np.zeros([nx,ny])
    ii[nx//2, ny//2] = 100
    
    ii = gaussian_filter(ii, 25)
    plt.imshow(ii)
    
    dx = 1e-06
    dy = 1e-06
    
    speed_test(ii, 1, 1)
    
