# -*- coding: utf-8 -*-

import numpy as np
from scipy.ndimage import gaussian_filter
from matplotlib import pyplot as plt
from felpy.utils.analysis_utils import get_windows
from felpy.analysis.statistics.correlation import norm_difference
from PIL import Image
def test_array(npts, sx, sy):
    
    arr1 = np.random.rand(npts, npts)
    

    
    
    arr2 = np.roll(arr1, sy, axis = 0)
    arr2 = np.roll(arr2, sx, axis = 1)
    
    #arr1 = gaussian_filter(arr1, 10)
    #arr2 = gaussian_filter(arr2, 10)
    return arr1, arr2
    

def odd_even(N):
    """
    check if a number is odd or even
    """
    
    if N % 2 == 0:
        ans = "even"
    elif N % 2 == 1:
        ans = 'odd'
    return ans
 

def get_window(arr, c, l):
    """
    return the window of analysis given a feature pixel centered at C
    
    :param c: coordinates (tuple of window center) (cx,cy)
    :param l: side length of square window
    """
    cx = c[1]
    cy = c[0]
    
    pad = l//2
    f = np.pad(arr, l)
    fl = l//2
    
    if odd_even(l) == 'odd':
        feature = f[(cx+pad-l):(cx+pad),
                    (cy+pad):(cy+pad+l)]
 
    elif odd_even(l) == 'even':
        feature = f[cx+pad-l//2:cx+pad+l//2,
                    cy+pad-l//2:cy+pad+l//2]
    return feature


def correlation(arr1, arr2):
    corr = 1 - norm_difference(arr1, arr2, plot = False).mean()
    return corr


def double_plot(arr1, arr2):
    
    vmin = min(np.min(arr1), np.min(arr2))
    vmax = min(np.max(arr1), np.max(arr2))
    
    fig, axs = plt.subplots(1,2)
    ax1, ax2 = axs
    ax1.imshow(arr1, vmin = 0, vmax = 255)
    ax1.set_title("reference")
    ax2.set_title("perturbed")
    ax2.imshow(arr2, vmin = 0, vmax = 255)
    plt.show()

def window_stack(a, stepsize=1, width=3):
    n = a.shape[0]
    return np.hstack( [a[i:1+n+i-width:stepsize] for i in range(0,width)] )

def rolling_window(a, shape):  # rolling window for 2D array
    s = (a.shape[0] - shape[0] + 1,) + (a.shape[1] - shape[1] + 1,) + shape
    strides = a.strides + a.strides
    return np.lib.stride_tricks.as_strided(a, shape=s, strides=strides)

    
if __name__ == '__main__':
    
    arr1, arr2 = test_array(100, sx = 4, sy = 0)
    #arr1 = np.asarray(Image.open('../../data/samples/test.jpg').convert('LA'))[:,:,0]
    #arr1 = np.arange(10000).reshape(100,100)
    
    print(arr1[25,25])
    double_plot(arr1,arr1)
    
 
    ### analysis window  
    #window = get_window(arr1, c = c, l = 24)
    c = [25,25]
    ### feature window
    

    f = rolling_window(np.pad(arr2, 2), shape = (5,5))
    fwindow = f[25,25,]
    
    w = rolling_window(arr1,(5,5))
    #double_plot(window, fwindow)
    ans = np.zeros([w.shape[0],w.shape[1]])
    
    for dx in range(w.shape[0]):
        for dy in range(w.shape[1]):
            ans[dx,dy] = correlation(w[dx,dy,], fwindow )
            
    ans = np.pad(ans, 2)
    print(np.where((ans == np.max(ans))))
# =============================================================================
#     corrs = np.zeros([window.shape[0],
#                       window.shape[1]])
#     
# 
# =============================================================================
