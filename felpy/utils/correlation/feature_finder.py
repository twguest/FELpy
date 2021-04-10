# -*- coding: utf-8 -*-

import numpy as np
from scipy.ndimage import gaussian_filter
from matplotlib import pyplot as plt
from felpy.utils.analysis_utils import get_windows
from felpy.analysis.statistics.correlation import norm_difference
from PIL import Image

def test_array(npts, sx, sy):
    
    arr1 = np.random.rand(npts, npts)
    arr2 = shift_image(arr1, sx, sy)
    #arr1 = gaussian_filter(arr1, 10)
    #arr2 = gaussian_filter(arr2, 10)
    return arr1, arr2
    
def shift_image(X, dx, dy):
    X = np.roll(X, dy, axis=0)
    X = np.roll(X, dx, axis=1)
    if dy>0:
        X[:dy, :] = 0
    elif dy<0:
        X[dy:, :] = 0
    if dx>0:
        X[:, :dx] = 0
    elif dx<0:
        X[:, dx:] = 0
    return X

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
    ax1.imshow(arr1, vmin = 0, vmax = 1)
    ax1.set_title("reference")
    ax2.set_title("perturbed")
    ax2.imshow(arr2, vmin = 0, vmax = 1)
    plt.show()

def window_stack(a, stepsize=1, width=3):
    n = a.shape[0]
    return np.hstack( [a[i:1+n+i-width:stepsize] for i in range(0,width)] )

def rolling_window(a, shape):  # rolling window for 2D array
    s = (a.shape[0] - shape[0] + 1,) + (a.shape[1] - shape[1] + 1,) + shape
    strides = a.strides + a.strides
    return np.lib.stride_tricks.as_strided(a, shape=s, strides=strides)

    
if __name__ == '__main__':
    
    arr1, arr2 = test_array(100, sx = 0, sy = 4)
    #arr1 = np.asarray(Image.open('../../data/samples/test.jpg').convert('LA'))[:,:,0]
    #arr1 = np.arange(10000).reshape(100,100)
    
    double_plot(arr1,arr1)
    
 
    ### analysis window  
    c = [25,25]
    window = get_window(arr1, c = c, l = 25)
    print(window.shape)
    ### feature window
    

    f = rolling_window(np.pad(arr2, 2), shape = (25,25))
    fwindow = f[25,25,]
    from skimage.registration import phase_cross_correlation
    print(phase_cross_correlation(fwindow, window,upsample_factor = 5))

# =============================================================================
#     corrs = np.zeros([window.shape[0],
#                       window.shape[1]])
#     
# 
# =============================================================================
